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

from constants import DIST_VECTOR_AXES
from logger import get_logger
from protein import Protein
from utils.qubit_utils import fix_qubits

logger = get_logger()


class DistanceMap:
    def __init__(self, protein: Protein):
        """
        Initializes the distance map for the given protein's main chain,
        setting up data structures to store distances along multiple axes
        and computing initial pairwise distances.
        """
        self.protein = protein
        self.distance_map = defaultdict(lambda: defaultdict(int))
        self._distances_vector: list[defaultdict] = [
            defaultdict(lambda: defaultdict(int)) for _ in range(DIST_VECTOR_AXES)
        ]
        self._calc_distances_main_chain()

    def _calc_distances_main_chain(self):
        """
        Computes the pairwise distances between all beads in the main chain.
        Distances are accumulated over multiple axes, squared, and adjusted
        using fixed qubit operators where applicable.
        """
        logger.debug("Creating distance map for main chain")

        main_chain_len = len(self.protein.main_chain)

        for lower_bead_idx in range(main_chain_len):
            for upper_bead_idx in range(lower_bead_idx + 1, main_chain_len):
                logger.debug(
                    "Processing pair main_chain_%s -> main_chain_%s",
                    lower_bead_idx,
                    upper_bead_idx,
                )
                for k in range(lower_bead_idx, upper_bead_idx):
                    indic_funcs = self.protein.main_chain[k].turn_funcs()
                    sub_lattice_sign = (-1) ** k

                    for indic_fun_x, dist_vector in zip(
                        indic_funcs, self._distances_vector
                    ):
                        dist_vector[lower_bead_idx][upper_bead_idx] += (
                            sub_lattice_sign * indic_fun_x
                        )

                for dist_vector in self._distances_vector:
                    dist_vector[lower_bead_idx][upper_bead_idx] = fix_qubits(
                        dist_vector[lower_bead_idx][upper_bead_idx]
                    )

                for dist_vector in self._distances_vector:
                    self.distance_map[lower_bead_idx][upper_bead_idx] += (
                        dist_vector[lower_bead_idx][upper_bead_idx] ** 2
                    )

                self.distance_map[lower_bead_idx][upper_bead_idx] = fix_qubits(
                    self.distance_map[lower_bead_idx][upper_bead_idx]
                )

                logger.debug(
                    f"main_chain_{lower_bead_idx} -> main_chain_{upper_bead_idx}: {self.distance_map[lower_bead_idx][upper_bead_idx]}"
                )
        logger.debug("Distance map initialized successfully")

    def __getitem__(self, key):
        """Returns the distance map entry for the given bead index."""
        return self.distance_map[key]

    def __setitem__(self, key, value):
        """Sets the distance map entry for the given bead index."""
        self.distance_map[key] = value

    def __len__(self):
        """Returns the number of beads in the distance map."""
        return len(self.distance_map)
