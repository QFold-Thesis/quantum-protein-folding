"""
Distance map calculations for a protein's main chain.

This module provides the ``DistanceMap`` class, which:

- computes pairwise distances between beads in the main chain,
- stores distances in a vectorized form for multiple axes,
- applies qubit fixes to account for predefined bead states,
- and maintains a dictionary of squared distances for downstream
  quantum simulations of protein folding.
"""

from collections import defaultdict

from qiskit.quantum_info import SparsePauliOp

from constants import DIST_VECTOR_AXES, EMPTY_OP_COEFF, QUBITS_PER_TURN
from logger import get_logger
from protein import Protein
from protein.bead.bead import Bead
from utils.qubit_utils import (
    build_identity_op,
    fix_qubits,
)

logger = get_logger()


class DistanceMap:
    def __init__(self, protein: Protein):
        """
        Initializes the distance map for the given protein's main chain,
        setting up data structures to store distances along multiple axes
        and computing initial pairwise distances.

        Args:
            protein (Protein): The Protein object that includes all information about protein.

        Raises:
            Exception: If distance calculation for the main chain fails.

        """
        self._protein: Protein = protein
        self._main_chain_len: int = len(self._protein.main_chain)

        self._pauli_op_len: int = (self._main_chain_len - 1) * QUBITS_PER_TURN
        self._distance_map: defaultdict[int, defaultdict[int, SparsePauliOp]] = (
            defaultdict(
                lambda: defaultdict(
                    lambda: build_identity_op(self._pauli_op_len, EMPTY_OP_COEFF)
                )
            )
        )

        self._main_chain_distances_detected: int = 0

        try:
            logger.debug("Initializing DistanceMap for MainChain...")
            self._calc_distances_main_chain()
        except Exception:
            logger.exception(
                "Error in initializing DistanceMap for MainChain"
            )
            raise
        else:
            logger.info(
                f"DistanceMap for MainChain initialized with {self._main_chain_distances_detected} distances detected"
            )

    def _calc_distances_main_chain(self) -> None:
        """
        Calculates pairwise quantum distances between beads in the main chain.

        For each bead pair (lower_bead_idx, upper_bead_idx), computes a vector of quantum operators
        representing the squared distance along each axis, using turn functions and sublattice signs.
        Results are stored in the distance map for use in quantum Hamiltonian construction.
        """
        for lower_bead_idx in range(self._main_chain_len):
            for upper_bead_idx in range(lower_bead_idx + 1, self._main_chain_len):

                lower_bead: Bead = self._protein.main_chain[lower_bead_idx]
                upper_bead: Bead = self._protein.main_chain[upper_bead_idx]
            
                axes_vector: list[SparsePauliOp] = [
                    build_identity_op(self._pauli_op_len, EMPTY_OP_COEFF)
                    for _ in range(DIST_VECTOR_AXES)
                ]

                for k in range(lower_bead_idx, upper_bead_idx):
                    indic_funcs = self._protein.main_chain[k].turn_funcs()
                    if indic_funcs is None:
                        logger.debug(
                            f"Skipping distance calculation between MainBeads: {lower_bead.symbol} (index: {lower_bead.index}) and {upper_bead.symbol} (index: {upper_bead.index}) due to undefined turn functions"
                        )
                        continue

                    sub_lattice_sign: int = (-1) ** k

                    for axis_idx, indic_fun_x in enumerate(indic_funcs):
                        axes_vector[axis_idx] += sub_lattice_sign * indic_fun_x

                for axis_idx in range(len(axes_vector)):
                    axes_vector[axis_idx] = fix_qubits(axes_vector[axis_idx])
                    self._distance_map[lower_bead_idx][upper_bead_idx] += (
                        axes_vector[axis_idx] ** 2
                    )

                self._distance_map[lower_bead_idx][upper_bead_idx] = fix_qubits(
                    self._distance_map[lower_bead_idx][upper_bead_idx]
                )
                self._main_chain_distances_detected += 1

                logger.debug(
                    f"Calculated distance operator between MainBeads: {lower_bead.symbol} (index: {lower_bead.index}) and {upper_bead.symbol} (index: {upper_bead.index})"
                )

    def __getitem__(self, key: int) -> defaultdict[int, SparsePauliOp]:
        """
        Get the distance map entry for a given bead index.

        Args:
            key (int): Index of the bead.

        Returns:
            defaultdict[int, SparsePauliOp]: Distance map entry for the bead.

        """
        return self._distance_map[key]

    def __setitem__(self, key: int, value: defaultdict[int, SparsePauliOp]) -> None:
        """
        Set the distance map entry for a given bead index.

        Args:
            key (int): Index of the bead.
            value (defaultdict[int, SparsePauliOp]): Distance map entry to assign.

        """
        self._distance_map[key] = value

    def __len__(self) -> int:
        """
        Return the number of beads in the distance map.

        Returns:
            int: Number of beads stored in the distance map.

        """
        return len(self._distance_map)
