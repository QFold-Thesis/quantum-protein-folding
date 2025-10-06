from __future__ import annotations

from collections import defaultdict
from datetime import datetime, timezone
import sys
from pathlib import Path
from typing import TYPE_CHECKING

from constants import DIST_VECTOR_AXES, QUBITS_PER_TURN, EMPTY_SIDECHAIN_PLACEHOLDER
from logger import get_logger
from utils.qubit_utils import create_empty_sparse_pauli_op, fix_qubits

if TYPE_CHECKING:  # pragma: no cover
    from protein import Protein
    from protein.bead import Bead
    from qiskit.quantum_info import SparsePauliOp

logger = get_logger()


class DistanceMap:
    def __init__(self, protein: Protein):
        self._protein: Protein = protein
        self._main_chain_len: int = len(self._protein.main_chain)

        self._pauli_op_len: int = (self._main_chain_len - 1) * QUBITS_PER_TURN
        self._distance_map: defaultdict[Bead, defaultdict[Bead, SparsePauliOp]] = (
            defaultdict(
                lambda: defaultdict(
                    lambda: create_empty_sparse_pauli_op(self._pauli_op_len)
                )
            )
        )

        self._main_chain_distances_detected: int = 0

        try:
            self._calc_distances_main_chain()
            self._add_distances_side_chain()
        except Exception:
            logger.exception(
                "Error occurred while calculating distances for main_chain"
            )
            raise
        else:
            logger.debug(
                f"Distance map for main_chain initialized successfully with {self._main_chain_distances_detected} distances detected."
            )
            # Dump distance map to a file and exit (debug requirement)
            try:
                self._dump_distance_map_to_file()
            except Exception:
                logger.exception("Failed to dump distance map to file before exit.")
            # Exit program immediately after dump as per user request
            sys.exit(0)

    def _calc_distances_main_chain(self) -> None:
        for lower_bead_idx in range(self._main_chain_len):
            for upper_bead_idx in range(lower_bead_idx + 1, self._main_chain_len):
                lower_bead = self._protein.main_chain[lower_bead_idx]
                upper_bead = self._protein.main_chain[upper_bead_idx]
                axes_vector: list[SparsePauliOp] = [
                    create_empty_sparse_pauli_op(self._pauli_op_len)
                    for _ in range(DIST_VECTOR_AXES)
                ]

                for k in range(lower_bead_idx, upper_bead_idx):
                    indic_funcs = self._protein.main_chain[k].turn_funcs()
                    if indic_funcs is None:
                        logger.debug(
                            f"No turn functions for bead {k}, skipping calculating distance..."
                        )
                        continue

                    sub_lattice_sign: int = (-1) ** k

                    for axis_idx, indic_fun_x in enumerate(indic_funcs):
                        axes_vector[axis_idx] += sub_lattice_sign * indic_fun_x

                for axis_idx in range(len(axes_vector)):
                    axes_vector[axis_idx] = fix_qubits(axes_vector[axis_idx])
                    self._distance_map[lower_bead][upper_bead] += (
                        axes_vector[axis_idx] ** 2
                    )

                self._distance_map[lower_bead][upper_bead] = fix_qubits(
                    self._distance_map[lower_bead][upper_bead]
                )
                self._main_chain_distances_detected += 1

                logger.debug(
                    "Calculated distance for main_chain_%d -> main_chain_%d | Num qubits: %d",
                    lower_bead_idx,
                    upper_bead_idx,
                    self._distance_map[lower_bead][upper_bead].num_qubits,
                )

    def _add_distances_side_chain(self) -> None:
        """Add distance contributions that involve side-chain beads (if present)."""
        logger.debug("siemaaa")
        # Quick exit if no side beads are present anywhere.
        for lower_bead_idx in range(self._main_chain_len):
            logger.debug(f"lower_bead_idx = {lower_bead_idx}")
            lower_main_bead = self._protein.main_chain[lower_bead_idx]
            logger.debug(f"lower_main_bead.side_chain[0] = {lower_main_bead.side_chain}")
            lower_side_bead = None
            if lower_main_bead.side_chain[0].symbol is not EMPTY_SIDECHAIN_PLACEHOLDER:
                lower_side_bead = lower_main_bead.side_chain.beads[0]
            lower_indic_funcs = lower_main_bead.turn_funcs()

            for upper_bead_idx in range(lower_bead_idx + 1, self._main_chain_len):
                logger.debug(f"upper_bead_idx = {upper_bead_idx}")
                upper_main_bead = self._protein.main_chain[upper_bead_idx]
                upper_side_bead = None
                if upper_main_bead.side_chain[0].symbol is not EMPTY_SIDECHAIN_PLACEHOLDER:
                    upper_side_bead = upper_main_bead.side_chain.beads[0]
                upper_indic_funcs = upper_main_bead.turn_funcs()

                for axis_idx in range(DIST_VECTOR_AXES):
                    lower_indic = (
                        lower_indic_funcs[axis_idx] if lower_indic_funcs is not None else None
                    )
                    logger.debug(f"lower_indic = {lower_indic}")
                    upper_indic = (
                        upper_indic_funcs[axis_idx] if upper_indic_funcs is not None else None
                    )
                    logger.debug(f"upper_indic = {upper_indic}")

                    # side(lower) -> main(upper)
                    if lower_side_bead and lower_indic is not None:
                        logger.debug("side(lower) -> main(upper)")
                        distance_term = self._calc_distance_term(
                            self._distance_map,
                            lower_bead_idx,
                            lower_indic,
                            upper_bead_idx,
                            None,
                        )
                        self._distance_map[lower_side_bead][upper_main_bead] += distance_term

                    # main(lower) -> side(upper)
                    if upper_side_bead and upper_indic is not None:
                        logger.debug("main(lower) -> side(upper)")
                        distance_term = self._calc_distance_term(
                            self._distance_map,
                            lower_bead_idx,
                            None,
                            upper_bead_idx,
                            upper_indic,
                        )
                        self._distance_map[lower_main_bead][upper_side_bead] += distance_term

                    # side(lower) -> side(upper)
                    if (
                        lower_side_bead
                        and upper_side_bead
                        and lower_indic is not None
                        and upper_indic is not None
                    ):
                        logger.debug("side(lower) -> side(upper)")
                        distance_term = self._calc_distance_term(
                            self._distance_map,
                            lower_bead_idx,
                            lower_indic,
                            upper_bead_idx,
                            upper_indic,
                        )
                        self._distance_map[lower_side_bead][upper_side_bead] += distance_term

    def _calc_distance_term(
        self,
        distance_map_axis_x,
        lower_bead_idx,
        lower_indic_func,
        upper_bead_idx,
        upper_indic_func,
    ) -> 'SparsePauliOp':
        lower_bead_obj = self._protein.main_chain[lower_bead_idx]
        upper_bead_obj = self._protein.main_chain[upper_bead_idx] 
        result = distance_map_axis_x[lower_bead_obj][upper_bead_obj]
        if lower_indic_func is not None:
            result -= (-1) ** lower_bead_idx * lower_indic_func
        if upper_indic_func is not None:
            result -= (-1) ** upper_bead_idx * upper_indic_func
        return fix_qubits(result)

    def __getitem__(self, key: 'Bead') -> defaultdict['Bead', 'SparsePauliOp']:  # type: ignore[override]
        return self._distance_map[key]

    def __setitem__(self, key: 'Bead', value: defaultdict['Bead', 'SparsePauliOp']) -> None:
        self._distance_map[key] = value

    def __len__(self) -> int:
        return len(self._distance_map)

    # ------------------------------------------------------------------
    # Debug / dumping utilities
    # ------------------------------------------------------------------
    def _dump_distance_map_to_file(self) -> None:
        """
        Serialize the distance map to a human-readable text file.

        Each line lists: lower_bead(index,symbol) upper_bead(index,symbol) num_qubits term_summary.
        The term summary is the string representation of the SparsePauliOp (may be long).
        """
        timestamp = datetime.now(tz=timezone.utc).strftime("%Y%m%d_%H%M%S")
        out_path = Path(f"distance_map_{timestamp}.txt")
        lines: list[str] = []
        lines.append("# Distance Map Dump\n")
        lines.append(f"# Main chain length: {self._main_chain_len}\n")
        lines.append(f"# Distances detected (main-chain pairs): {self._main_chain_distances_detected}\n\n")

        for lower_bead, inner in self._distance_map.items():
            for upper_bead, op in inner.items():
                if op is None:
                    continue
                try:
                    op_str = str(op)
                except Exception as exc:  # pragma: no cover - defensive
                    op_str = f"<error obtaining op string: {exc}>"
                lines.append(
                    f"{lower_bead.index}:{lower_bead.symbol} "
                    f"{upper_bead.index}:{upper_bead.symbol} "
                    f"num_qubits={op.num_qubits} op={op_str}"
                )

        out_path.write_text("\n".join(lines), encoding="utf-8")
        logger.info("Distance map dumped to %s", out_path.resolve())
