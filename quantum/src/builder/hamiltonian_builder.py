from typing import TYPE_CHECKING

from qiskit.quantum_info import SparsePauliOp

from constants import BOUNDING_CONSTANT, MJ_ENERGY_MULTIPLIER
from contact.contact_map import ContactMap
from distance.distance_map import DistanceMap
from enums import Penalties
from interaction.interaction import Interaction
from logger import get_logger
from protein import Protein
from protein.bead import Bead
from utils.qubit_utils import build_full_identity, fix_qubits, pad_to_n_qubits

if TYPE_CHECKING:
    from protein.chain import MainChain

logger = get_logger()


class HamiltonianBuilder:
    def __init__(
        self,
        protein: Protein,
        interaction: Interaction,
        distance_map: DistanceMap,
        contact_map: ContactMap,
    ):
        self.protein: Protein = protein
        self.interaction: Interaction = interaction
        self.distance_map: DistanceMap = distance_map
        self.contact_map: ContactMap = contact_map

    def sum_hamiltonians(self) -> SparsePauliOp:
        """Build and sum all Hamiltonian components, padding to a common qubit size."""
        backbone: SparsePauliOp = self.build_backbone_contact_term()
        backtrack: SparsePauliOp = self.add_backtracking_penalty()

        part_hamiltonians: list[SparsePauliOp] = [backbone, backtrack]
        target_qubits: int = max(
            hamiltonian.num_qubits for hamiltonian in part_hamiltonians
        )
        padded_hamiltonians: list[SparsePauliOp] = [
            pad_to_n_qubits(hamiltonian, target_qubits)
            for hamiltonian in part_hamiltonians
        ]

        total_hamiltonian: SparsePauliOp = 0
        for hamiltonian in padded_hamiltonians:
            total_hamiltonian += hamiltonian

        return total_hamiltonian.simplify()

    def build_backbone_contact_term(self) -> SparsePauliOp:
        """
        Builds the Hamiltonian term corresponding to backbone_backbone (BB-BB) interactions.
        Includes both 1st neighbor and 2nd neighbor contributions (with shifts i±1, j±1).
        """
        logger.info("Creating h_bbbb term (BB-BB interactions)...")

        main_chain: MainChain = self.protein.main_chain
        hamiltonian: SparsePauliOp = 0
        chain_len: int = len(main_chain)

        for i in range(len(main_chain) - 4):
            for j in range(i + 4, len(main_chain)):
                if (j - i) % 2 == 0:
                    continue

                if 0 <= i < chain_len and 0 <= j < chain_len:
                    logger.debug(f"Adding BB-BB i={i}, j={j} (1st neighbor)")
                    hamiltonian += self.contact_map.main_main_contacts[i][
                        j
                    ] ^ self.get_first_neighbor_hamiltonian(
                        i, j, Penalties.OVERLAP_PENALTY
                    )

                for di, dj in [
                    (-1, 0),
                    (1, 0),
                    (0, -1),
                    (0, 1),
                ]:
                    ii, jj = i + di, j + dj
                    if 0 <= ii < chain_len and 0 <= jj < chain_len:
                        logger.debug(f"Adding BB-BB i={ii}, j={jj} (2nd neighbor)")
                        hamiltonian += self.contact_map.main_main_contacts[i][
                            j
                        ] ^ self.get_second_neighbor_hamiltonian(
                            ii, jj, Penalties.OVERLAP_PENALTY
                        )

                hamiltonian = fix_qubits(hamiltonian)

        logger.info(f"Finished creating h_bbbb term: {hamiltonian}")
        return hamiltonian

    def add_backtracking_penalty(self) -> SparsePauliOp:
        main_chain: MainChain = self.protein.main_chain
        h_back: SparsePauliOp = 0
        for i in range(1, len(main_chain) - 2):
            h_back += Penalties.BACK_PENALTY * self.get_turn_operators(
                main_chain[i], main_chain[i + 1]
            )

        return fix_qubits(h_back)

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
        symbol_lower: str = self.protein.main_chain.get_symbol_at(lower_bead_idx)
        symbol_upper: str = self.protein.main_chain.get_symbol_at(upper_bead_idx)
        energy: float = self.interaction.get_energy(symbol_lower, symbol_upper)
        x: SparsePauliOp = self.distance_map[lower_bead_idx][upper_bead_idx]
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
        symbol_lower: str = self.protein.main_chain.get_symbol_at(lower_bead_idx)
        symbol_upper: str = self.protein.main_chain.get_symbol_at(upper_bead_idx)
        energy: float = self.interaction.get_energy(symbol_lower, symbol_upper)
        x: SparsePauliOp = self.distance_map[lower_bead_idx][upper_bead_idx]
        expression: SparsePauliOp = lambda_1 * (
            2 * build_full_identity(x.num_qubits) - x
        ) + (MJ_ENERGY_MULTIPLIER * energy * build_full_identity(x.num_qubits))
        return fix_qubits(expression)
