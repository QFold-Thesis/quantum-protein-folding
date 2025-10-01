"""
Contact map utilities for protein folding simulations.

This module provides the `ContactMap` class that:

- Initializes contact operators between beads of a protein's main chain,
- Ensures minimum bond separation between beads to consider contacts,
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
    """
    Represents a contact map for a protein's main chain.

    Stores pairwise contact operators between beads, using Pauli operators
    to encode whether contacts are present, respecting minimum bond distances
    and sublattice constraints.
    """

    def __init__(self, protein: Protein):
        """Initializes the contact map for the given protein."""
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

        self.contacts_detected: int = 0
        self._protein: Protein = protein

        self._num_contact_qubits: int = pow(len(self._protein.main_chain) - 1, 2)
        self._full_identity: SparsePauliOp = build_identity_op(
            num_qubits=self._num_contact_qubits
        )

        try:
            self._initialize_contact_map()
        except Exception:
            logger.exception("Error in initializing contact map")
            raise

        else:
            logger.debug(
                f"Contact map initialized successfully. Contacts detected: {self.contacts_detected}"
            )

    def _initialize_contact_map(self):
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

                self.main_main_contacts[lower_bead.index][upper_bead.index] = (
                    contact_operator
                )
                self.contacts_detected += 1

    def _create_main_main_contact(
        self, upper_bead: Bead, lower_bead: Bead
    ) -> SparsePauliOp:
        """Creates a contact operator between two main chain beads."""
        z_op_index: int = (lower_bead.index) * (len(self._protein.main_chain) - 1) + (
            upper_bead.index
        )

        contact_operator = build_pauli_z_operator(
            num_qubits=self._num_contact_qubits, pauli_z_indices={z_op_index}
        )

        logger.debug(
            f"Created main-main contact between beads main_chain_{upper_bead.index} -> main_chain_{lower_bead.index} | Z index: {z_op_index} | Num qubits: {self._num_contact_qubits}"
        )
        return convert_to_qubits(contact_operator)
