from qiskit.quantum_info import SparsePauliOp

from constants import BOUNDING_CONSTANT, MJ_ENERGY_MULTIPLIER
from distance.distance_map import DistanceMap
from interaction.mj_interaction import MJInteraction
from logger import get_logger
from protein import Protein
from protein.bead import Bead
from utils.qubit_utils import build_full_identity, fix_qubits

logger = get_logger()


class HamiltonianBuilder:
    def __init__(self, protein: Protein):
        self.protein: Protein = protein
        self.mj: MJInteraction = MJInteraction(protein)
        self.distance_map: DistanceMap = DistanceMap(protein)

    def get_turn_operators(self, lower_bead: Bead, upper_bead: Bead) -> SparsePauliOp:
        turn_operators: SparsePauliOp = sum(
            lower_bead_idx @ upper_bead_idx
            for lower_bead_idx, upper_bead_idx in zip(
                lower_bead.turn_funcs(), upper_bead.turn_funcs()
            )
        )

        return fix_qubits(turn_operators)

    def get_first_neighbor_hamiltonian(
        self,
        lower_bead_idx: int,
        upper_bead_idx: int,
        lambda_1: float,
    ) -> SparsePauliOp:
        lambda_0: float = (
            BOUNDING_CONSTANT * (upper_bead_idx - lower_bead_idx + 1) * lambda_1
        )
        energy: float = self.mj.get_energy_by_indices(lower_bead_idx, upper_bead_idx)
        x: int = self.distance_map[lower_bead_idx][upper_bead_idx]
        expression: SparsePauliOp = lambda_0 * (
            x - build_full_identity(x.num_qubits)
        ) + (MJ_ENERGY_MULTIPLIER * energy * build_full_identity(x.num_qubits))
        return fix_qubits(expression)

    def get_second_neighbor_hamiltonian(
        self,
        lower_bead_idx: int,
        upper_bead_idx: int,
        lambda_1: float,
    ) -> SparsePauliOp:
        energy: float = self.mj.get_energy_by_indices(lower_bead_idx, upper_bead_idx)
        x: int = self.distance_map[lower_bead_idx][upper_bead_idx]
        expression: SparsePauliOp = lambda_1 * (
            2 * build_full_identity(x.num_qubits) - x
        ) + (MJ_ENERGY_MULTIPLIER * energy * build_full_identity(x.num_qubits))
        return fix_qubits(expression)
