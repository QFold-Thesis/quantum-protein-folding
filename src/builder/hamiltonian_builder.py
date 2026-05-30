"""Utilities for building the hamiltonian of a protein for quantum simulations.

This module provides the HamiltonianBuilder class, which constructs hamiltonian
operators for a given protein, including backbone interactions, backtracking
penalties, neighbor-based contact terms, and an optional external-field bias term.

The external-field term (Variant A) models a position-dependent energy landscape
by mapping each bead's sequential index ``i`` to a coordinate key ``(i,)`` and
looking up the corresponding energy from an :class:`~particle.ExternalField`
instance.  The resulting contribution is a scalar multiple of the identity:

    H_field = Σ_i  E_field((i,)) · I

This shifts the global energy proportionally to the sum of field energies over
the chain, differentiating between conformations only when combined with other
terms that break the chain's symmetry.  Full spatial coupling δ(r_i, r_field)
that depends on the actual lattice position of each bead belongs to Variant B.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from constants import (
    BOUNDING_CONSTANT,
    EMPTY_OP_COEFF,
    MJ_ENERGY_MULTIPLIER,
    QUBITS_PER_TURN,
)
from enums import Penalties
from exceptions import InvalidOperatorError
from logger import get_logger
from utils.qubit_utils import (
    build_identity_op,
    fix_qubits,
    pad_to_n_qubits,
)

if TYPE_CHECKING:
    from qiskit.quantum_info import SparsePauliOp

    from contact.contact_map import ContactMap
    from distance.distance_map import DistanceMap
    from interaction.interaction import Interaction
    from particle.external_field import ExternalField
    from protein import Protein
    from protein.bead import Bead
    from protein.chain import _MainChain

logger = get_logger()


class HamiltonianBuilder:
    """Constructs hamiltonian operators for a given protein, including backbone interactions and backtracking penalties.

    Optionally accepts an :class:`~particle.ExternalField` to include an
    external-field bias term ``H_field`` in the total Hamiltonian.  When
    *external_field* is ``None`` the builder behaves exactly as before and
    produces identical results to the field-free case.

    Attributes:
        protein (Protein): The Protein object that includes all information about protein.
        interaction (Interaction): Interaction model between beads of the protein.
        distance_map (DistanceMap): Matrix of pairwise distances between residues.
        contact_map (ContactMap): Matrix indicating residue-residue contacts.
        external_field (ExternalField | None): Optional external interaction field.

    """

    def __init__(
        self,
        protein: Protein,
        interaction: Interaction,
        distance_map: DistanceMap,
        contact_map: ContactMap,
        external_field: ExternalField | None = None,
    ) -> None:
        """Initializes the HamiltonianBuilder with required protein data and interaction maps.

        Args:
            protein (Protein): The Protein object that includes all information about protein.
            interaction (Interaction): Interaction model between beads of the protein.
            distance_map (DistanceMap): Matrix of pairwise distances between residues.
            contact_map (ContactMap): Matrix indicating residue-residue contacts.
            external_field (ExternalField | None, optional): External interaction field that
                contributes a per-bead energy bias to the Hamiltonian.  Pass ``None`` (the
                default) to omit the field term entirely, preserving backward compatibility.

        """
        self.protein: Protein = protein
        self.interaction: Interaction = interaction
        self.distance_map: DistanceMap = distance_map
        self.contact_map: ContactMap = contact_map
        self.external_field: ExternalField | None = external_field

    def sum_hamiltonians(self) -> SparsePauliOp:
        """Build and sum all hamiltonian components, padding to a common qubit size.

        Constructs the backbone, backtracking, and (optionally) external-field
        terms, checks qubit consistency, pads them to the same qubit count, and
        sums them into a single hamiltonian.

        Note:
            The padding step ensures that all SparsePauliOp operators have the same
            number of qubits, which is required for valid operator addition.
            The total hamiltonian is initialized with an identity operator and then
            each padded component is added sequentially.

            When :attr:`external_field` is ``None`` the result is identical to the
            field-free case; ``H_field`` contributes a zero operator in that case.

        Returns:
            SparsePauliOp: The total hamiltonian operator, simplified and ready for use.

        Raises:
            InvalidOperatorError: If any part hamiltonian has `num_qubits` set to None.

        """
        logger.debug("Started process of building total hamiltonian...")
        h_backbone: SparsePauliOp = self._build_backbone_contact_term()
        h_backtrack: SparsePauliOp = self._add_backtracking_penalty()
        h_field: SparsePauliOp = self._build_external_field_term()

        part_hamiltonians: list[SparsePauliOp] = [h_backbone, h_backtrack, h_field]

        for idx, hamiltonian in enumerate(part_hamiltonians):
            if hamiltonian.num_qubits is None:
                msg: str = f"Hamiltonian of part {idx} has num_qubits set to None"
                raise InvalidOperatorError(msg)

        target_qubits: int = max(
            int(hamiltonian.num_qubits) for hamiltonian in part_hamiltonians
        )
        logger.debug(
            "Target qubits count for the final hamiltonian to be padded to: %s",
            target_qubits,
        )

        padded_hamiltonians: list[SparsePauliOp] = [
            pad_to_n_qubits(hamiltonian, target_qubits)
            for hamiltonian in part_hamiltonians
        ]

        total_hamiltonian: SparsePauliOp = build_identity_op(
            target_qubits, EMPTY_OP_COEFF
        )
        for hamiltonian in padded_hamiltonians:
            total_hamiltonian += hamiltonian

        logger.info("Finished building total hamiltonian.")
        return total_hamiltonian.simplify()

    def _build_backbone_contact_term(self) -> SparsePauliOp:
        """Builds the hamiltonian term corresponding to backbone_backbone (BB-BB) interactions. Includes both 1st neighbor and 2nd neighbor contributions (with shifts i±1, j±1).

        Note:
            Only pairs that belong to different sublattices are considered for first-neighbor interactions.
            For each valid pair, two contributions are added: one for first-neighbor interactions and one for
            second-neighbor interactions with nearby beads. The second-neighbor contribution applies a penalty
            to avoid overlaps.

        Returns:
            SparsePauliOp: hamiltonian term representing BB-BB interactions.

        """
        logger.debug(
            "Creating hamiltonian term of backbone-backbone (BB-BB) contacts..."
        )

        main_chain: _MainChain = self.protein.main_chain
        chain_len: int = len(main_chain)

        h_backbone_num_qubits: int = (
            pow((chain_len - 1), 2) + (chain_len - 1) * QUBITS_PER_TURN
        )
        h_backbone: SparsePauliOp = build_identity_op(
            h_backbone_num_qubits, EMPTY_OP_COEFF
        )

        for i in range(len(main_chain) - 4):
            for j in range(i + 4, len(main_chain)):
                if (j - i) % 2 == 0:
                    continue

                if 0 <= i < chain_len and 0 <= j < chain_len:
                    logger.debug(
                        "Adding backbone-backbone contact between Bead (index %s) and Bead (index %s) [1st neighbor contact]",
                        i,
                        j,
                    )
                    h_backbone += self.contact_map.main_main_contacts[i][
                        j
                    ] ^ self.get_first_neighbor_hamiltonian(
                        i, j, float(Penalties.OVERLAP_PENALTY)
                    )

                for di, dj in [
                    (-1, 0),
                    (1, 0),
                    (0, -1),
                    (0, 1),
                ]:
                    ii, jj = i + di, j + dj
                    if 0 <= ii < chain_len and 0 <= jj < chain_len:
                        logger.debug(
                            "Adding backbone-backbone contact between Bead (index %s) and Bead (index %s) [2nd neighbor contact]",
                            ii,
                            jj,
                        )
                        h_backbone += self.contact_map.main_main_contacts[i][
                            j
                        ] ^ self.get_second_neighbor_hamiltonian(
                            ii, jj, float(Penalties.OVERLAP_PENALTY)
                        )

                h_backbone = fix_qubits(h_backbone)

        logger.info(
            "Finished creating hamiltonian term of backbone-backbone (BB-BB) contacts with %s qubits.",
            h_backbone.num_qubits,
        )
        return h_backbone

    def _build_external_field_term(self) -> SparsePauliOp:
        """Builds the Hamiltonian bias term from the external field (Variant A).

        When :attr:`external_field` is ``None`` this method returns a zero identity
        operator matching the turn-qubit register size and the builder behaves
        identically to the field-free case.

        Otherwise, the method iterates over all beads in the main chain and sums
        per-bead energy contributions::

            H_field = Σ_i  E_field((i,)) · I

        Each bead's sequential index ``i`` is used as its lattice coordinate key
        ``(i,)``.  This means a non-uniform :class:`~particle.ExternalField` can
        assign different energies to different positions along the sequence (e.g.
        a stronger field at the active site residues), while a uniform field
        simply shifts the total energy by a constant.

        Note:
            This implementation corresponds to Variant A of the external-field
            extension.  Variant B will encode the actual spatial position of each
            bead as a quantum degree of freedom, enabling a true
            δ(r_i, r_field) contact term.

        Returns:
            SparsePauliOp: Hamiltonian bias term H_field, or a zero identity
            operator when no field is present.

        """
        n_turn_qubits: int = (len(self.protein.main_chain) - 1) * QUBITS_PER_TURN

        if self.external_field is None:
            logger.debug(
                "No external field provided - H_field set to zero identity (%d qubits).",
                n_turn_qubits,
            )
            return build_identity_op(n_turn_qubits, EMPTY_OP_COEFF)

        logger.debug("Building external field Hamiltonian term (Variant A)...")

        h_field: SparsePauliOp = build_identity_op(n_turn_qubits, EMPTY_OP_COEFF)
        total_field_energy: float = 0.0

        for bead in self.protein.main_chain:
            bead_coord: tuple[int] = (bead.index,)
            energy: float = self.external_field.get_energy(bead_coord)
            h_field = h_field + energy * build_identity_op(n_turn_qubits)
            total_field_energy += energy
            logger.debug(
                "Bead %s (index=%d): field energy = %s",
                bead.symbol,
                bead.index,
                energy,
            )

        logger.info(
            "Finished building H_field: total field energy = %s over %d beads.",
            total_field_energy,
            len(self.protein.main_chain),
        )
        return h_field.simplify()

    def _add_backtracking_penalty(self) -> SparsePauliOp:
        """Adds a penalty term to the hamiltonian to discourage backtracking in the main chain configuration.

        Returns:
            SparsePauliOp: hamiltonian term representing backtracking penalties.

        """
        logger.debug("Creating hamiltonian term of backtracking penalty...")

        main_chain: _MainChain = self.protein.main_chain
        h_backtrack_num_qubits: int = (len(main_chain) - 1) * QUBITS_PER_TURN
        h_backtrack: SparsePauliOp = build_identity_op(
            h_backtrack_num_qubits, EMPTY_OP_COEFF
        )

        for i in range(1, len(main_chain) - 2):
            logger.debug(
                "Adding backtracking penalty between Bead (index %s) and Bead (index %s)",
                i,
                i + 1,
            )
            h_backtrack += Penalties.BACK_PENALTY * self.get_turn_operators(
                main_chain[i], main_chain[i + 1]
            )

        logger.info(
            "Finished creating hamiltonian term of backtracking penalty with %s qubits.",
            h_backtrack.num_qubits,
        )
        return fix_qubits(h_backtrack)

    def get_turn_operators(self, lower_bead: Bead, upper_bead: Bead) -> SparsePauliOp:
        """Builds the combined turn operators for two consecutive beads in the main chain.

        Generates a quantum operator representing allowed directional turns
        between two beads based on their turn functions. If either bead lacks
        defined turn functions, an identity operator is returned.

        Args:
            lower_bead (Bead): The bead from the main chain at the lower index.
            upper_bead (Bead): The bead from the main chain at the upper index.

        Returns:
            SparsePauliOp: Combined turn operator describing the interaction between the two beads.

        """
        lower_turn_funcs: (
            None | tuple[SparsePauliOp, SparsePauliOp, SparsePauliOp, SparsePauliOp]
        ) = lower_bead.turn_funcs()
        upper_turn_funcs: (
            None | tuple[SparsePauliOp, SparsePauliOp, SparsePauliOp, SparsePauliOp]
        ) = upper_bead.turn_funcs()

        if lower_turn_funcs is None or upper_turn_funcs is None:
            logger.info(
                "One of the beads has no turn functions defined. Returning identity operator instead"
            )
            return build_identity_op(
                (len(self.protein.main_chain) - 1) * QUBITS_PER_TURN,
                EMPTY_OP_COEFF,
            )

        turn_operators: SparsePauliOp = build_identity_op(
            (len(self.protein.main_chain) - 1) * QUBITS_PER_TURN, EMPTY_OP_COEFF
        )

        for lower_bead_idx, upper_bead_idx in zip(
            lower_turn_funcs, upper_turn_funcs, strict=True
        ):
            turn_operators += lower_bead_idx @ upper_bead_idx

        return fix_qubits(turn_operators)

    def get_first_neighbor_hamiltonian(
        self,
        lower_bead_idx: int,
        upper_bead_idx: int,
        lambda_1: float,
    ) -> SparsePauliOp:
        """Computes the hamiltonian contribution for first-neighbor bead pairs, combining distance-based and interaction contact energies.

        Note:
             lambda_0 combines the bounding constant, bead separation, and lambda_1
            to scale the distance-based penalty. MJ_ENERGY_MULTIPLIER scales the
            contribution from the Miyazawa-Jernigan interaction energy.

        Args:
            lower_bead_idx (int): Index of the lower bead in the main chain.
            upper_bead_idx (int): Index of the upper bead in the main chain.
            lambda_1 (float): Penalty coefficient for first neighbor interaction.

        Returns:
            SparsePauliOp: Quantum operator representing the first neighbor hamiltonian term.

        Raises:
            InvalidOperatorError: If the number of qubits in the operator is None.

        """
        lambda_0: float = (
            BOUNDING_CONSTANT * (upper_bead_idx - lower_bead_idx + 1) * lambda_1
        )
        symbol_lower: str = self.protein.main_chain.get_symbol_at(lower_bead_idx)
        symbol_upper: str = self.protein.main_chain.get_symbol_at(upper_bead_idx)

        energy: float = self.interaction.get_energy(symbol_lower, symbol_upper)
        x: SparsePauliOp = self.distance_map[lower_bead_idx][upper_bead_idx]

        if x.num_qubits is None:
            msg: str = "x.num_qubits is None, cannot build first neighbor hamiltonian."
            raise InvalidOperatorError(msg)

        expression: SparsePauliOp = lambda_0 * (x - build_identity_op(x.num_qubits)) + (
            MJ_ENERGY_MULTIPLIER * energy * build_identity_op(x.num_qubits)
        )

        return fix_qubits(expression)

    def get_second_neighbor_hamiltonian(
        self,
        lower_bead_idx: int,
        upper_bead_idx: int,
        lambda_1: float,
    ) -> SparsePauliOp:
        """Computes the hamiltonian contribution for second-neighbor bead pairs, including distance-based and interaction terms.

        Args:
            lower_bead_idx (int): Index of the lower bead in the main chain.
            upper_bead_idx (int): Index of the upper bead in the main chain.
            lambda_1 (float): Penalty coefficient for second neighbor interaction.

        Returns:
            SparsePauliOp: Quantum operator representing the second neighbor hamiltonian term.

        Raises:
            InvalidOperatorError: If the number of qubits in the operator is None.

        """
        symbol_lower: str = self.protein.main_chain.get_symbol_at(lower_bead_idx)
        symbol_upper: str = self.protein.main_chain.get_symbol_at(upper_bead_idx)

        energy: float = self.interaction.get_energy(symbol_lower, symbol_upper)
        x: SparsePauliOp = self.distance_map[lower_bead_idx][upper_bead_idx]

        if x.num_qubits is None:
            msg: str = "x.num_qubits is None, cannot build second neighbor hamiltonian."
            raise InvalidOperatorError(msg)

        expression: SparsePauliOp = lambda_1 * (
            2 * build_identity_op(x.num_qubits) - x
        ) + (MJ_ENERGY_MULTIPLIER * energy * build_identity_op(x.num_qubits))

        return fix_qubits(expression)
