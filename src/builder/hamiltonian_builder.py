"""Utilities for building the hamiltonian of a protein for quantum simulations.

This module provides the HamiltonianBuilder class, which constructs hamiltonian
operators for a given protein, including backbone interactions, backtracking
penalties, and neighbor-based contact terms, using distance and interaction maps.

Two optional terms extend the folding Hamiltonian:

    H_total = H_backbone + H_backtrack + H_field + H_ligand

``H_field`` couples an :class:`~particle.external_field.ExternalField` to the
lattice positions of the residues, and ``H_ligand`` couples a free
:class:`~particle.ligand.Ligand` to the residues it sits next to. Both are
*spatial*: they read the conformation rather than adding a constant, so they
reorder the spectrum and genuinely move the optimal fold. Omitting both
reproduces the plain folding Hamiltonian bit for bit.

Register layout
---------------
Optional registers are appended above the existing ones, which leaves the
backbone operators untouched and keeps results backward compatible::

    [ protein turns | backbone contacts | ligand walk | ligand contacts ]
      (N-1)*Q         (N-1)^2             steps*Q       eligible residues

Ligand coupling
---------------
For each residue the ligand may bind, a boolean contact qubit ``c_i`` states
whether the ligand claims to touch residue ``i``. The claim is scored against
the geometry actually encoded in the turn qubits::

    H_ligand = sum_i c_i * [ E(L, aa_i) + lambda_c * (d2(i, L) - 1)^2 ]
             + lambda_u * (sum_i c_i - 1)^2

The squared deviation is zero exactly when the ligand is a lattice nearest
neighbour of residue ``i`` and positive otherwise, so a false claim is penalised
rather than rewarded - unlike a linear ``(d2 - 1)`` penalty, which pays out when
the two particles overlap. The uniqueness term pins the ligand to exactly one
binding partner, removing the degenerate "ligand drifts away at zero energy"
states.

Excluded volume needs a separate term. Contact qubits only exist for residues on
the *opposite* sublattice, because those are the only ones the ligand can be a
nearest neighbour of; residues on its *own* sublattice are exactly the ones it
can sit on top of, and nothing above forbids that. Their squared distance to the
ligand is an even number, so the Lagrange polynomial through the reachable
shells isolates the overlapping one::

    H_exclusion = lambda_x * sum_j prod_{s in {2, 4}} (d2(j, L) - s) / (0 - s)

This evaluates to 1 when the ligand occupies residue ``j``'s site and to 0 on the
two nearest allowed shells. Beyond those shells it grows again, which acts as a
weak confinement keeping the ligand in the neighbourhood of the chain - harmless
here, since a ligand encoded with a handful of steps cannot travel far anyway.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from constants import (
    BOUNDING_CONSTANT,
    EMPTY_OP_COEFF,
    LATTICE_CONTACT_DISTANCE,
    LIGAND_CONTACT_PENALTY,
    LIGAND_ENERGY_MULTIPLIER,
    LIGAND_EXCLUSION_PENALTY,
    LIGAND_UNIQUENESS_PENALTY,
    MJ_ENERGY_MULTIPLIER,
    QUBITS_PER_TURN,
    SAME_SUBLATTICE_SHELLS,
)
from enums import Penalties
from exceptions import InvalidOperatorError
from logger import get_logger
from utils.lattice_utils import build_chain_position, build_squared_distance
from utils.qubit_utils import (
    build_identity_op,
    build_turn_qubit,
    fix_qubits,
    pad_to_n_qubits,
)

if TYPE_CHECKING:
    from qiskit.quantum_info import SparsePauliOp

    from contact.contact_map import ContactMap
    from distance.distance_map import DistanceMap
    from interaction.interaction import Interaction
    from interaction.ligand_interaction import LigandInteraction
    from particle.external_field import ExternalField
    from particle.ligand import Ligand
    from protein import Protein
    from protein.bead import Bead
    from protein.chain import _MainChain

logger = get_logger()


class HamiltonianBuilder:
    """Constructs hamiltonian operators for a given protein, including backbone interactions and backtracking penalties.

    Optionally couples an external field and a free ligand to the chain. Both are
    off by default, in which case the builder reproduces the plain folding
    Hamiltonian exactly.

    Attributes:
        protein (Protein): The Protein object that includes all information about protein.
        interaction (Interaction): Interaction model between beads of the protein.
        distance_map (DistanceMap): Matrix of pairwise distances between residues.
        contact_map (ContactMap): Matrix indicating residue-residue contacts.
        external_field (ExternalField | None): Field acting on residue positions.
        ligand (Ligand | None): Free particle sharing the lattice with the chain.
        ligand_interaction (LigandInteraction | None): Ligand-residue energies.

    """

    def __init__(
        self,
        protein: Protein,
        interaction: Interaction,
        distance_map: DistanceMap,
        contact_map: ContactMap,
        external_field: ExternalField | None = None,
        ligand: Ligand | None = None,
        ligand_interaction: LigandInteraction | None = None,
    ) -> None:
        """Initializes the HamiltonianBuilder with required protein data and interaction maps.

        Args:
            protein (Protein): The Protein object that includes all information about protein.
            interaction (Interaction): Interaction model between beads of the protein.
            distance_map (DistanceMap): Matrix of pairwise distances between residues.
            contact_map (ContactMap): Matrix indicating residue-residue contacts.
            external_field (ExternalField | None, optional): Field acting on the
                residues' lattice positions. Defaults to None.
            ligand (Ligand | None, optional): Free particle to place on the
                lattice alongside the chain. Defaults to None.
            ligand_interaction (LigandInteraction | None, optional): Energies
                between the ligand and each residue. Required when *ligand* is
                given. Defaults to None.

        Raises:
            ValueError: If a ligand is supplied without its interaction model.

        """
        if ligand is not None and ligand_interaction is None:
            msg: str = "ligand_interaction is required when a ligand is supplied"
            raise ValueError(msg)

        self.protein: Protein = protein
        self.interaction: Interaction = interaction
        self.distance_map: DistanceMap = distance_map
        self.contact_map: ContactMap = contact_map
        self.external_field: ExternalField | None = external_field
        self.ligand: Ligand | None = ligand
        self.ligand_interaction: LigandInteraction | None = ligand_interaction

    @property
    def num_protein_qubits(self) -> int:
        """int: Width of the register used by the protein-only terms."""
        chain_len: int = len(self.protein.main_chain)
        return pow(chain_len - 1, 2) + (chain_len - 1) * QUBITS_PER_TURN

    @property
    def ligand_walk_offset(self) -> int:
        """int: Index of the first qubit encoding the ligand's lattice position."""
        return self.num_protein_qubits

    @property
    def ligand_contact_offset(self) -> int:
        """int: Index of the first qubit flagging a claimed ligand contact."""
        if self.ligand is None:
            return self.num_protein_qubits
        return self.num_protein_qubits + self.ligand.num_walk_qubits

    @property
    def num_qubits(self) -> int:
        """int: Total register width, including any ligand registers."""
        if self.ligand is None:
            return self.num_protein_qubits

        eligible: list[int] = self.ligand.eligible_bead_indices(
            len(self.protein.main_chain)
        )
        return self.ligand_contact_offset + len(eligible)

    def sum_hamiltonians(self) -> SparsePauliOp:
        """Build and sum all hamiltonian components, padding to a common qubit size.

        Constructs the backbone and backtracking terms plus any optional field and
        ligand terms, checks qubit consistency, pads them to the same qubit count,
        and sums them into a single hamiltonian.

        Note:
            The padding step ensures that all SparsePauliOp operators have the same
            number of qubits, which is required for valid operator addition.
            The total hamiltonian is initialized with an identity operator and then
            each padded component is added sequentially.

            Optional registers sit above the protein's own qubits, so padding the
            protein-only terms with identities keeps them exactly as they were.

        Returns:
            SparsePauliOp: The total hamiltonian operator, simplified and ready for use.

        Raises:
            InvalidOperatorError: If any part hamiltonian has `num_qubits` set to None.

        """
        logger.debug("Started process of building total hamiltonian...")
        h_backbone: SparsePauliOp = self._build_backbone_contact_term()
        h_backtrack: SparsePauliOp = self._add_backtracking_penalty()

        part_hamiltonians: list[SparsePauliOp] = [h_backbone, h_backtrack]

        for idx, hamiltonian in enumerate(part_hamiltonians):
            if hamiltonian.num_qubits is None:
                msg: str = f"Hamiltonian of part {idx} has num_qubits set to None"
                raise InvalidOperatorError(msg)

        target_qubits: int = max(
            *(int(hamiltonian.num_qubits) for hamiltonian in part_hamiltonians),
            self.num_qubits,
        )
        logger.debug(
            "Target qubits count for the final hamiltonian to be padded to: %s",
            target_qubits,
        )

        padded_hamiltonians: list[SparsePauliOp] = [
            pad_to_n_qubits(hamiltonian, target_qubits)
            for hamiltonian in part_hamiltonians
        ]

        padded_hamiltonians.append(self._build_external_field_term(target_qubits))
        padded_hamiltonians.append(self._build_ligand_term(target_qubits))

        total_hamiltonian: SparsePauliOp = build_identity_op(
            target_qubits, EMPTY_OP_COEFF
        )
        for hamiltonian in padded_hamiltonians:
            total_hamiltonian += hamiltonian

        logger.info("Finished building total hamiltonian.")
        return total_hamiltonian.simplify()

    def _build_external_field_term(self, num_qubits: int) -> SparsePauliOp:
        """Builds the external field contribution to the hamiltonian.

        Args:
            num_qubits (int): Width of the register the operator must span.

        Returns:
            SparsePauliOp: The field term, or a zero operator when no field is
            configured.

        """
        if self.external_field is None:
            logger.debug("No external field configured - H_field omitted.")
            return build_identity_op(num_qubits, EMPTY_OP_COEFF)

        field_term: SparsePauliOp = self.external_field.build_hamiltonian(
            chain_length=len(self.protein.main_chain), num_qubits=num_qubits
        )

        return fix_qubits(field_term)

    def _build_ligand_term(self, num_qubits: int) -> SparsePauliOp:
        """Builds the ligand contribution to the hamiltonian.

        Couples the ligand's lattice position to the residues it claims to touch,
        scoring each claim against the geometry encoded in the turn qubits. See
        the module docstring for the exact form.

        Args:
            num_qubits (int): Width of the register the operator must span.

        Returns:
            SparsePauliOp: The ligand term, or a zero operator when no ligand is
            configured.

        """
        if self.ligand is None or self.ligand_interaction is None:
            logger.debug("No ligand configured - H_ligand omitted.")
            return build_identity_op(num_qubits, EMPTY_OP_COEFF)

        logger.debug("Creating hamiltonian term of ligand-residue contacts...")

        chain_len: int = len(self.protein.main_chain)
        eligible: list[int] = self.ligand.eligible_bead_indices(chain_len)

        ligand_position: list[SparsePauliOp] = self.ligand.position_operators(
            num_qubits=num_qubits, walk_qubit_offset=self.ligand_walk_offset
        )

        identity: SparsePauliOp = build_identity_op(num_qubits)
        contact_distance: SparsePauliOp = LATTICE_CONTACT_DISTANCE * identity

        ligand_term: SparsePauliOp = build_identity_op(num_qubits, EMPTY_OP_COEFF)
        contact_count: SparsePauliOp = build_identity_op(num_qubits, EMPTY_OP_COEFF)

        for slot, bead_index in enumerate(eligible):
            contact_qubit: SparsePauliOp = build_turn_qubit(
                z_index=self.ligand_contact_offset + slot, num_qubits=num_qubits
            )
            contact_count = contact_count + contact_qubit

            residue_position: list[SparsePauliOp] = build_chain_position(
                bead_index=bead_index, num_qubits=num_qubits
            )
            squared_distance: SparsePauliOp = build_squared_distance(
                residue_position, ligand_position
            )

            symbol: str = self.protein.main_chain.get_symbol_at(bead_index)
            energy: float = self.ligand_interaction.get_energy(symbol)

            deviation: SparsePauliOp = squared_distance - contact_distance
            claim_penalty: SparsePauliOp = LIGAND_CONTACT_PENALTY * (
                deviation @ deviation
            )

            ligand_term = ligand_term + (
                contact_qubit
                @ ((LIGAND_ENERGY_MULTIPLIER * energy * identity) + claim_penalty)
            )

            logger.debug(
                "Ligand %s may bind residue %s (index %s) with energy %s",
                self.ligand.symbol,
                symbol,
                bead_index,
                energy,
            )

        uniqueness_deviation: SparsePauliOp = contact_count - identity
        ligand_term = ligand_term + LIGAND_UNIQUENESS_PENALTY * (
            uniqueness_deviation @ uniqueness_deviation
        )

        ligand_term = ligand_term + self._build_exclusion_term(
            num_qubits=num_qubits, ligand_position=ligand_position
        )

        logger.info(
            "Finished creating H_ligand for %s over %d candidate residues on %s qubits.",
            self.ligand.symbol,
            len(eligible),
            num_qubits,
        )
        return fix_qubits(ligand_term.simplify())

    def _build_exclusion_term(
        self, num_qubits: int, ligand_position: list[SparsePauliOp]
    ) -> SparsePauliOp:
        """Builds the penalty keeping the ligand off the residues' lattice sites.

        Only residues sharing the ligand's sublattice can coincide with it, and
        those are precisely the residues that have no contact qubit to constrain
        them. See the module docstring for the polynomial used.

        Args:
            num_qubits (int): Width of the register the operator must span.
            ligand_position (list[SparsePauliOp]): Axis-coefficient operators of
                the ligand's position.

        Returns:
            SparsePauliOp: The excluded-volume penalty term.

        """
        if self.ligand is None:
            return build_identity_op(num_qubits, EMPTY_OP_COEFF)

        chain_len: int = len(self.protein.main_chain)
        identity: SparsePauliOp = build_identity_op(num_qubits)
        exclusion_term: SparsePauliOp = build_identity_op(num_qubits, EMPTY_OP_COEFF)

        overlapping_residues: list[int] = [
            index
            for index in range(chain_len)
            if index % 2 == self.ligand.sublattice_parity
        ]

        for bead_index in overlapping_residues:
            residue_position: list[SparsePauliOp] = build_chain_position(
                bead_index=bead_index, num_qubits=num_qubits
            )
            squared_distance: SparsePauliOp = build_squared_distance(
                residue_position, ligand_position
            )

            overlap_indicator: SparsePauliOp = identity
            for shell in SAME_SUBLATTICE_SHELLS[1:]:
                overlap_indicator = overlap_indicator @ (
                    (squared_distance - shell * identity) * (1.0 / -shell)
                )

            exclusion_term = exclusion_term + (
                LIGAND_EXCLUSION_PENALTY * overlap_indicator
            )

        logger.debug(
            "Built ligand excluded-volume penalty over residues %s.",
            overlapping_residues,
        )
        return exclusion_term.simplify()

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
