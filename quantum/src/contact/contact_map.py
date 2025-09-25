from __future__ import annotations

from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

from logger import get_logger
from utils.qubit_utils import build_full_identity, build_pauli_z_operator, convert_to_qubits

from typing import TYPE_CHECKING

if TYPE_CHECKING:  # Imports used only for type checking to reduce runtime dependencies here
    from qiskit.quantum_info import SparsePauliOp  # pyright: ignore[reportMissingTypeStubs]
    from protein import Protein
    from protein.bead import Bead

logger = get_logger()


MIN_DISTANCE_BETWEEN_CONTACTS = (
    5  # Minimum bonds between two beads to consider a contact
)


class ContactMap:
    def __init__(self, protein: Protein):
        self.main_main_contacts: dict[int, dict[int, SparsePauliOp]] = defaultdict(
            lambda: defaultdict()
        )

        """
        # self.main_side_contacts: dict[int, dict[int, SparsePauliOp]] = defaultdict(
        #     lambda: defaultdict()
        # )
        # self.side_main_contacts: dict[int, dict[int, SparsePauliOp]] = defaultdict(
        #     lambda: defaultdict()
        # )
        # self.side_side_contacts: dict[int, dict[int, SparsePauliOp]] = defaultdict(
        #     lambda: defaultdict()
        # )
        """

        self._contacts_detected: int = 0
        self._protein: Protein = protein

        self._num_contact_qubits: int = pow(len(self._protein.main_chain) - 1, 2)
        self._full_identity: SparsePauliOp = build_full_identity(
            num_qubits=self._num_contact_qubits
        )

        try:
            self.initialize_contact_map()
        except Exception:
            logger.exception("Error in initializing contact map")
            raise

        else:
            logger.debug(
                f"Contact map initialized successfully. Contacts detected: {self._contacts_detected}"
            )

        # Automatically dump a pretty-printed version of the main-main contact map
        try:
            self.dump_main_main_contacts()
        except Exception:  # pragma: no cover - logging side effect only
            logger.exception("Failed to pretty print main-main contact map")

    def initialize_contact_map(self):
        """Initializes all contact maps to empty dictionaries."""
        main_chain_length: int = len(self._protein.main_chain)

        for lower_bead_idx in range(main_chain_length - 2):
            for upper_bead_idx in range(lower_bead_idx + 2, main_chain_length):
                upper_bead: Bead = self._protein.main_chain[upper_bead_idx]
                lower_bead: Bead = self._protein.main_chain[lower_bead_idx]

                if upper_bead.sublattice == lower_bead.sublattice:
                    continue

                if (
                    abs(upper_bead.index - lower_bead.index)
                    < MIN_DISTANCE_BETWEEN_CONTACTS
                ):
                    continue

                contact_operator: SparsePauliOp = self._create_main_main_contact(
                    upper_bead=upper_bead, lower_bead=lower_bead
                )

                self.main_main_contacts[upper_bead.index][lower_bead.index] = (
                    contact_operator
                )
                self._contacts_detected += 1

    def _create_main_main_contact(
        self, upper_bead: Bead, lower_bead: Bead
    ) -> SparsePauliOp:
        """Creates a contact operator between two main chain beads."""
        z_op_index: int = (lower_bead.index) * (len(self._protein.main_chain) - 1) + (
            upper_bead.index
        )
        logger.debug(f"cwel z_op_index for {upper_bead.index}, {lower_bead.index}: {z_op_index}")

        contact_operator = build_pauli_z_operator(
            num_qubits=self._num_contact_qubits, pauli_z_indices={z_op_index}
        )
        logger.debug(f"Contact operator before conversion: {contact_operator}")

        logger.debug(
            f"Created main-main contact between beads {upper_bead.index} and {lower_bead.index} | Z index: {z_op_index} | Num qubits: {self._num_contact_qubits}"
        )
        return convert_to_qubits(contact_operator)

    # ---------------------------------------------------------------------
    # Pretty Printing / Export
    # ---------------------------------------------------------------------
    def dump_main_main_contacts(self, file_path: str | None = None) -> Path:
        """
        Pretty print the main-main contact map to a text file.

        The format is a matrix where rows and columns correspond to main chain bead
        indices (with their residue symbol). A contact is marked with 'X', absence
        with '.'. Only contacts that pass initialization constraints are marked.

    Additionally, an explicit list of contact triples and a raw (plain) dump of
    the nested main_main_contacts dictionary (with SparsePauliOp string
    representations) are appended at the end for debugging purposes.

        Parameters
        ----------
        file_path: str | None
            Optional path to the output file. If omitted, a file named
            'contact_map_main_main.txt' is created in the current working directory.

        Returns
        -------
        Path
            Path to the written file.

        """
        # Resolve output path
        output_path = Path(file_path) if file_path else Path("contact_map_main_main.txt")

        main_chain = self._protein.main_chain
        n = len(main_chain)

        # Build a quick lookup set for faster membership tests
        contact_pairs: set[tuple[int, int]] = set()
        for upper_idx, lower_dict in self.main_main_contacts.items():
            for lower_idx in lower_dict:
                # Store as (min,max) so lookup order agnostic if we ever change iteration
                contact_pairs.add((lower_idx, upper_idx))

        header_indices = [f"{i}" for i in range(n)]
        header_symbols = [bead.symbol for bead in main_chain]

        lines: list[str] = []
        lines.append("# Main-Main Contact Map")
        lines.append(
            f"# Generated: {datetime.now(tz=timezone.utc).isoformat(timespec='seconds')}"
        )
        lines.append(f"# Protein length: {n}")
        lines.append(f"# Contacts detected: {self._contacts_detected}")
        lines.append("# Legend: X = contact, . = no contact (matrix is symmetric; diagonal blank)")
        lines.append("")
        lines.append("Indices : " + " ".join(f"{idx:>3}" for idx in header_indices))
        lines.append("Symbols : " + " ".join(f"{sym:>3}" for sym in header_symbols))
        lines.append("")

        # Matrix body
        for i, bead_i in enumerate(main_chain):
            row_cells: list[str] = []
            for j, _ in enumerate(main_chain):
                if i == j:
                    row_cells.append("   ")
                    continue
                lower, upper = sorted((i, j))
                marker = " X " if (lower, upper) in contact_pairs else " . "
                row_cells.append(marker)
            lines.append(f"{i:>2}({bead_i.symbol}) |" + "".join(row_cells))

        # Detailed list (optional): provide explicit contacts with derived Z index
        lines.append("")
        lines.append("# Explicit contact pairs (upper_index, lower_index, z_op_index)")
        for upper_idx, lower_dict in sorted(self.main_main_contacts.items()):
            for lower_idx in sorted(lower_dict):
                z_index = lower_idx * (len(main_chain) - 1) + upper_idx
                lines.append(f"{upper_idx:>3}, {lower_idx:>3}, {z_index:>5}")

        # Raw / plain representation of the dictionary (may be verbose)
        lines.append("")
        lines.append("# Raw main_main_contacts dump (upper_index -> lower_index: SparsePauliOp)")
        for upper_idx in sorted(self.main_main_contacts):
            lines.append(f"upper {upper_idx}:")
            for lower_idx in sorted(self.main_main_contacts[upper_idx]):
                op = self.main_main_contacts[upper_idx][lower_idx]
                lines.append(f"  lower {lower_idx}: {op}")

        output_path.write_text("\n".join(lines), encoding="utf-8")
        logger.info(f"Main-main contact map written to {output_path.resolve()}")
        return output_path
