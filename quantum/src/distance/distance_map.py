from collections import defaultdict

from qiskit.quantum_info import SparsePauliOp

from constants import DIST_VECTOR_AXES, QUBITS_PER_TURN
from logger import get_logger
from protein import Protein
from utils.qubit_utils import create_empty_sparse_pauli_op, fix_qubits

logger = get_logger()


class DistanceMap:
    def __init__(self, protein: Protein):
        self._protein: Protein = protein
        self._main_chain_len: int = len(self._protein.main_chain)

        self._pauli_op_len: int = (self._main_chain_len - 1) * QUBITS_PER_TURN
        self._distance_map: defaultdict[int, defaultdict[int, SparsePauliOp]] = (
            defaultdict(
                lambda: defaultdict(
                    lambda: create_empty_sparse_pauli_op(self._pauli_op_len)
                )
            )
        )

        self.distances_detected: int = 0

        self._calc_distances_main_chain()

    def _calc_distances_main_chain(self) -> None:
        logger.debug("Creating distance map for main chain")

        for lower_bead_idx in range(self._main_chain_len):
            for upper_bead_idx in range(lower_bead_idx + 1, self._main_chain_len):
                logger.debug(
                    f"Started calculating distance for main_chain_{lower_bead_idx} -> main_chain_{upper_bead_idx}"
                )

                axes_vector: list[SparsePauliOp] = [
                    create_empty_sparse_pauli_op(self._pauli_op_len)
                    for _ in range(DIST_VECTOR_AXES)
                ]

                for k in range(lower_bead_idx, upper_bead_idx):
                    indic_funcs = self._protein.main_chain[k].turn_funcs()
                    if indic_funcs is None:
                        logger.info("No turn functions for bead %s, skipping...", k)
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
                self.distances_detected += 1

                logger.info(
                    f"Calculated distance: main_chain_{lower_bead_idx} -> main_chain_{upper_bead_idx}: {self._distance_map[lower_bead_idx][upper_bead_idx]}"
                )
        logger.info("Distance map initialized successfully")

    def __getitem__(self, key: int) -> defaultdict[int, SparsePauliOp]:
        return self._distance_map[key]

    def __setitem__(self, key: int, value: defaultdict[int, SparsePauliOp]) -> None:
        self._distance_map[key] = value

    def __len__(self) -> int:
        return len(self._distance_map)
