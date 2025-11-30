"""Contact map utilities for protein folding simulations.

This module provides the `ContactMap` class that:

- Initializes contact operators between MainBeads of a Protein's MainChain,
- Ensures minimum bond separation between MainBeads to consider contacts,
- Represents contacts as qubit operators (SparsePauliOp) suitable for quantum simulations.
"""

from collections import defaultdict

from qiskit.quantum_info import SparsePauliOp

from constants import MIN_DISTANCE_BETWEEN_CONTACTS
from logger import get_logger
from protein import Protein
from protein.bead import Bead
from utils.qubit_utils import (
    build_identity_op,
    build_pauli_z_operator,
    convert_to_qubits,
)

logger = get_logger()


class ContactMap:
    """Represents a contact map for a protein's main chain.

    Stores pairwise contact operators between MainBeads, using Pauli operators
    to encode whether contacts are present, respecting minimum bond distances
    and sublattice constraints.

    Attributes:
        main_main_contacts (dict[int, dict[int, SparsePauliOp]]): Contact operators between main chain beads.
        contacts_detected (int): Total number of contacts detected in the map.

    """

    def __init__(self, protein: Protein) -> None:
        """Initializes the contact map for the given protein.

        Args:
            protein (Protein): The Protein object that includes all information about protein.

        Raises:
            Exception: If contact map initialization fails.

        """
        self.main_main_contacts: dict[int, dict[int, SparsePauliOp]] = defaultdict(
            lambda: defaultdict()
        )

        self.contacts_detected: int = 0
        self._protein: Protein = protein

        self._num_contact_qubits: int = pow(len(self._protein.main_chain) - 1, 2)
        self._full_identity: SparsePauliOp = build_identity_op(
            num_qubits=self._num_contact_qubits
        )

        try:
            logger.debug("Initializing ContactMap...")
            self._initialize_contact_map()
        except Exception:
            logger.exception("Error in initializing ContactMap")
            raise
        else:
            logger.info(
                "ContactMap initialized with %s contacts detected.",
                self.contacts_detected,
            )

    def _initialize_contact_map(self) -> None:
        """Initializes all contact maps to empty dictionaries.

        Note:
            The minimum distance between residues for forming a contact is set by the constant
            MIN_DISTANCE_BETWEEN_CONTACTS = 5. This ensures that contacts are only considered
            between residues separated by at least five positions in the sequence, which is
            consistent with the geometric constraints of the tetrahedral lattice representation.

        """
        main_main_contacts_count: int = 0
        main_chain_length: int = len(self._protein.main_chain)
        logger.debug("Initializing MainBead-MainBead contacts...")

        for lower_bead_idx in range(main_chain_length - 2):
            for upper_bead_idx in range(lower_bead_idx + 2, main_chain_length):
                upper_bead: Bead = self._protein.main_chain[upper_bead_idx]
                lower_bead: Bead = self._protein.main_chain[lower_bead_idx]
                logger.debug(
                    "Evaluating potential contact between MainBeads: %s (index: %s) and %s (index: %s)",
                    upper_bead.symbol,
                    upper_bead.index,
                    lower_bead.symbol,
                    lower_bead.index,
                )

                if upper_bead.sublattice == lower_bead.sublattice:
                    logger.debug(
                        "Skipping contact between MainBeads: %s (index: %s) and %s (index: %s) due to same sublattice \n",
                        upper_bead.symbol,
                        upper_bead.index,
                        lower_bead.symbol,
                        lower_bead.index,
                    )
                    continue

                if (
                    abs(upper_bead.index - lower_bead.index)
                    < MIN_DISTANCE_BETWEEN_CONTACTS
                ):
                    logger.debug(
                        "Skipping contact between MainBeads: %s (index: %s) and %s (index: %s) due to insufficient bond separation (min=%s, actual=%s)\n",
                        upper_bead.symbol,
                        upper_bead.index,
                        lower_bead.symbol,
                        lower_bead.index,
                        MIN_DISTANCE_BETWEEN_CONTACTS,
                        abs(upper_bead.index - lower_bead.index),
                    )
                    continue

                contact_operator: SparsePauliOp = self._create_main_main_contact(
                    upper_bead=upper_bead, lower_bead=lower_bead
                )

                self.main_main_contacts[lower_bead.index][upper_bead.index] = (
                    contact_operator
                )
                self.contacts_detected += 1
                main_main_contacts_count += 1
                logger.debug(" ")

        logger.info(
            "Calculated %s MainBead-MainBead contacts", main_main_contacts_count
        )

    def _create_main_main_contact(
        self, upper_bead: Bead, lower_bead: Bead
    ) -> SparsePauliOp:
        """Creates a contact operator between two main chain MainBeads.

        Args:
            lower_bead (Bead): The bead from the main chain at the lower index.
            upper_bead (Bead): The bead from the main chain at the upper index.

        Returns:
            SparsePauliOp: Pauli-Z operator for the contact between the two MainBeads.

        """
        z_op_index: int = (lower_bead.index) * (len(self._protein.main_chain) - 1) + (
            upper_bead.index
        )

        contact_operator = build_pauli_z_operator(
            num_qubits=self._num_contact_qubits, pauli_z_indices={z_op_index}
        )

        logger.debug(
            "Created contact operator between MainBeads: %s (index: %s) and %s (index: %s)",
            lower_bead.symbol,
            lower_bead.index,
            upper_bead.symbol,
            upper_bead.index,
        )
        return convert_to_qubits(contact_operator)
