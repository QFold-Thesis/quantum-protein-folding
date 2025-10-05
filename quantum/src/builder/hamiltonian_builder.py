from __future__ import annotations

from typing import TYPE_CHECKING

from constants import BOUNDING_CONSTANT, MJ_ENERGY_MULTIPLIER, QUBITS_PER_TURN
from enums import Penalties
from exceptions import InvalidOperatorError
from logger import get_logger
from utils.qubit_utils import (
    build_full_identity,
    create_empty_sparse_pauli_op,
    fix_qubits,
    pad_to_n_qubits,
)

if TYPE_CHECKING:
    from qiskit.quantum_info import SparsePauliOp

    from contact.contact_map import ContactMap
    from distance.distance_map import DistanceMap
    from interaction.interaction import Interaction
    from protein import Protein
    from protein.bead import Bead
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
        h_backbone: SparsePauliOp = self._build_backbone_contact_term()
        h_backtrack: SparsePauliOp = self._add_backtracking_penalty()

        part_hamiltonians: list[SparsePauliOp] = [h_backbone, h_backtrack]

        for hamiltonian in part_hamiltonians:
            if hamiltonian.num_qubits is None:
                msg = "One of the part Hamiltonians has num_qubits set to None."
                raise InvalidOperatorError(msg)

        target_qubits: int = max(
            int(hamiltonian.num_qubits)
            for hamiltonian in part_hamiltonians
            if hamiltonian.num_qubits is not None
        )

        padded_hamiltonians: list[SparsePauliOp] = [
            pad_to_n_qubits(hamiltonian, target_qubits)
            for hamiltonian in part_hamiltonians
        ]

        total_hamiltonian: SparsePauliOp = create_empty_sparse_pauli_op(target_qubits)
        for hamiltonian in padded_hamiltonians:
            total_hamiltonian += hamiltonian

        return total_hamiltonian.simplify()

    def _build_backbone_contact_term(self) -> SparsePauliOp:
        """
        Builds the Hamiltonian term corresponding to backbone_backbone (BB-BB) interactions.
        Includes both 1st neighbor and 2nd neighbor contributions (with shifts i±1, j±1).
        """
        logger.info("Creating h_backbone term (BB-BB interactions)")

        main_chain: MainChain = self.protein.main_chain
        chain_len: int = len(main_chain)

        h_backbone_num_qubits: int = (
            pow((chain_len - 1), 2) + (chain_len - 1) * QUBITS_PER_TURN
        )
        h_backbone: SparsePauliOp = create_empty_sparse_pauli_op(h_backbone_num_qubits)

        for i in range(len(main_chain) - 4):
            for j in range(i + 4, len(main_chain)):
                if (j - i) % 2 == 0:
                    continue

                if 0 <= i < chain_len and 0 <= j < chain_len:
                    logger.debug(f"Adding BB-BB i={i}, j={j} (1st neighbor)")
                    h_backbone += self.contact_map.main_main_contacts[i][
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
                        h_backbone += self.contact_map.main_main_contacts[i][
                            j
                        ] ^ self.get_second_neighbor_hamiltonian(
                            ii, jj, Penalties.OVERLAP_PENALTY
                        )

                h_backbone = fix_qubits(h_backbone)

        logger.info(
            f"Finished creating h_backbone term with {h_backbone.num_qubits} qubits."
        )
        return h_backbone

    def _add_backtracking_penalty(self) -> SparsePauliOp:
        logger.debug("Creating h_backtrack term")
        main_chain: MainChain = self.protein.main_chain

        h_backtrack_num_qubits: int = (len(main_chain) - 1) * QUBITS_PER_TURN
        h_backtrack: SparsePauliOp = create_empty_sparse_pauli_op(
            h_backtrack_num_qubits
        )

        for i in range(1, len(main_chain) - 2):
            logger.debug(f"Adding backtracking penalty between beads {i} and {i + 1}")
            h_backtrack += Penalties.BACK_PENALTY * self.get_turn_operators(
                main_chain[i], main_chain[i + 1]
            )

        logger.debug(
            f"Finished creating h_backtrack term with {h_backtrack.num_qubits} qubits."
        )
        return fix_qubits(h_backtrack)

    def get_turn_operators(self, lower_bead: Bead, upper_bead: Bead) -> SparsePauliOp:
        lower_turn_funcs: (
            None | tuple[SparsePauliOp, SparsePauliOp, SparsePauliOp, SparsePauliOp]
        ) = lower_bead.turn_funcs()
        upper_turn_funcs: (
            None | tuple[SparsePauliOp, SparsePauliOp, SparsePauliOp, SparsePauliOp]
        ) = upper_bead.turn_funcs()

        if lower_turn_funcs is None or upper_turn_funcs is None:
            logger.debug(
                f"One of the beads {lower_bead.symbol}|{lower_bead.index} or {upper_bead.symbol}|{upper_bead.index} has no turn functions. Skipping turn operator calculation."
            )
            return create_empty_sparse_pauli_op(
                (len(self.protein.main_chain) - 1) * QUBITS_PER_TURN
            )

        turn_operators: SparsePauliOp = create_empty_sparse_pauli_op(
            (len(self.protein.main_chain) - 1) * QUBITS_PER_TURN
        )

        for lower_bead_idx, upper_bead_idx in zip(lower_turn_funcs, upper_turn_funcs):
            turn_operators += lower_bead_idx @ upper_bead_idx

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

        if x.num_qubits is None:
            msg = "x.num_qubits is None, cannot build first neighbor Hamiltonian."
            raise InvalidOperatorError(msg)

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

        if x.num_qubits is None:
            msg = "x.num_qubits is None, cannot build second neighbor Hamiltonian."
            raise InvalidOperatorError(msg)

        expression: SparsePauliOp = lambda_1 * (
            2 * build_full_identity(x.num_qubits) - x
        ) + (MJ_ENERGY_MULTIPLIER * energy * build_full_identity(x.num_qubits))

        return fix_qubits(expression)
