from collections import defaultdict

from constants import DIST_VECTOR_AXES
from logger import get_logger
from protein import Protein
from utils.qubit_utils import fix_qubits

logger = get_logger()


class DistanceMap:
    def __init__(self, protein: Protein):
        self._protein = protein
        self.distance_map = defaultdict(lambda: defaultdict(int))
        self._calc_distances_main_chain()

    def _calc_distances_main_chain(self):
        logger.debug("Creating distance map for main chain")

        main_chain_len = len(self._protein.main_chain)

        for lower_bead_idx in range(main_chain_len):
            for upper_bead_idx in range(lower_bead_idx + 1, main_chain_len):
                logger.debug(
                    "Processing pair main_chain_%s -> main_chain_%s",
                    lower_bead_idx,
                    upper_bead_idx,
                )

                axes_vector: list[int] = [0] * DIST_VECTOR_AXES

                for k in range(lower_bead_idx, upper_bead_idx):
                    indic_funcs = self._protein.main_chain[k].turn_funcs()
                    sub_lattice_sign = (-1) ** k

                    for axis_idx, indic_fun_x in enumerate(indic_funcs):
                        axes_vector[axis_idx] += sub_lattice_sign * indic_fun_x

                for axis_idx in range(len(axes_vector)):
                    axes_vector[axis_idx] = fix_qubits(axes_vector[axis_idx])
                    self.distance_map[lower_bead_idx][upper_bead_idx] += (
                        axes_vector[axis_idx] ** 2
                    )

                self.distance_map[lower_bead_idx][upper_bead_idx] = fix_qubits(
                    self.distance_map[lower_bead_idx][upper_bead_idx]
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
