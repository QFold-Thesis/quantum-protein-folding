from collections import defaultdict
from typing import TYPE_CHECKING

from qiskit.quantum_info import SparsePauliOp

from constants import DIST_VECTOR_AXES
from logger import get_logger
from protein import Protein
from utils.qubit_utils import fix_qubits

if TYPE_CHECKING:
    from protein.bead import Bead

logger = get_logger()


class DistanceMap:
    def __init__(self, protein: Protein):
        self.protein = protein
        self.distance_map = defaultdict(lambda: defaultdict(int))
        self._distances_vector: list[defaultdict] = [
            defaultdict(lambda: defaultdict(int)) for _ in range(DIST_VECTOR_AXES)
        ]
        self._calc_distances_main_chain()

    def _calc_distances_main_chain(self):
        logger.debug("Creating distance map for main chain")

        main_chain_len = len(self.protein.main_chain)

        for lower_bead_idx in range(1, main_chain_len):
            for upper_bead_idx in range(lower_bead_idx + 1, main_chain_len + 1):
                lower_bead: Bead = self.protein.main_chain[lower_bead_idx - 1]
                upper_bead: Bead = self.protein.main_chain[upper_bead_idx - 1]

                logger.debug(
                    "Processing pair main_chain_%s -> main_chain_%s",
                    lower_bead_idx,
                    upper_bead_idx,
                )
                for k in range(lower_bead_idx, upper_bead_idx):
                    indic_funcs = self.protein.main_chain[k - 1].turn_funcs()
                    sub_lattice_sign = (-1) ** k

                    for indic_fun_x, dist_vector in zip(
                        indic_funcs, self._distances_vector
                    ):
                        dist_vector[lower_bead][upper_bead] += (
                            sub_lattice_sign * indic_fun_x
                        )

                for dist_vector in self._distances_vector:
                    dist_vector[lower_bead][upper_bead] = fix_qubits(
                        dist_vector[lower_bead][upper_bead]
                    )

                for dist_vector in self._distances_vector:
                    self.distance_map[lower_bead_idx][upper_bead_idx] += (
                        dist_vector[lower_bead][upper_bead] ** 2
                    )

                self.distance_map[lower_bead_idx][upper_bead_idx] = fix_qubits(
                    self.distance_map[lower_bead_idx][upper_bead_idx]
                )

                if not isinstance(
                    self.distance_map[lower_bead_idx][upper_bead_idx], SparsePauliOp
                ):
                    print(
                        40 * "-",
                        self.distance_map[lower_bead_idx][upper_bead_idx],
                        40 * "-",
                    )

                logger.debug(
                    f"main_chain_{lower_bead_idx} -> main_chain_{upper_bead_idx}: {self.distance_map[lower_bead_idx][upper_bead_idx]}"
                )
        logger.debug("Distance map initialized successfully")

    def __getitem__(self, key):
        return self.distance_map[key]

    def __setitem__(self, key, value):
        self.distance_map[key] = value

    def __len__(self):
        return len(self.distance_map)
