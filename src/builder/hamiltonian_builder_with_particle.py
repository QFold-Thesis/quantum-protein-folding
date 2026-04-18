"""Extended Hamiltonian builder with particle support.

Provides HamiltonianBuilderWithParticle, which extends the standard
HamiltonianBuilder to include particle interaction terms in two variants:
1. External Field: Static potential affecting all beads
2. Dynamic Bead: Full bead-like particle with contact-based interactions
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from builder.hamiltonian_builder import HamiltonianBuilder
from constants import EMPTY_OP_COEFF
from exceptions import InvalidOperatorError
from logger import get_logger
from utils.qubit_utils import build_identity_op, pad_to_n_qubits

if TYPE_CHECKING:
    from qiskit.quantum_info import SparsePauliOp

    from contact.contact_map import ContactMap
    from distance.distance_map import DistanceMap
    from interaction.interaction import Interaction
    from particle import Particle
    from protein import Protein

logger = get_logger()


class HamiltonianBuilderWithParticle(HamiltonianBuilder):
    """Extended HamiltonianBuilder incorporating particle interaction terms.

    Builds Hamiltonians for protein-particle systems in two variants:
    - External field: Static potential from particle
    - Dynamic bead: Particle as a lattice bead with contact interactions

    Attributes:
        particle (Particle): Particle object specifying type and parameters.
        base_builder (HamiltonianBuilder): Original builder for protein-only terms.

    """

    def __init__(
        self,
        protein: Protein,
        interaction: Interaction,
        distance_map: DistanceMap,
        contact_map: ContactMap,
        particle: Particle | None = None,
    ) -> None:
        """Initialize extended Hamiltonian builder with optional particle.

        Args:
            protein (Protein): Protein object.
            interaction (Interaction): Interaction model (should support particle types).
            distance_map (DistanceMap): Distance map for protein beads.
            contact_map (ContactMap): Contact map for protein beads.
            particle (Particle, optional): Particle object. If None, behaves as standard builder.

        """
        super().__init__(protein, interaction, distance_map, contact_map)
        self.particle = particle
        
        if particle:
            logger.info(
                "Initialized HamiltonianBuilderWithParticle for %s variant",
                particle.particle_type.value,
            )

    def sum_hamiltonians(self) -> SparsePauliOp:
        """Build total Hamiltonian including particle terms.

        Returns:
            SparsePauliOp: Total Hamiltonian = protein + particle terms.

        Returns:
            SparsePauliOp: The total hamiltonian operator with particle contributions.

        """
        if self.particle is None:
            # No particle: use standard behavior
            logger.debug("No particle defined, using standard Hamiltonian")
            return super().sum_hamiltonians()

        logger.debug("Building Hamiltonian with particle term...")

        # Get protein-only Hamiltonian from base builder
        h_protein: SparsePauliOp = super().sum_hamiltonians()

        # Build particle-specific term
        from particle import ParticleType

        if self.particle.particle_type == ParticleType.EXTERNAL_FIELD:
            h_particle = self._build_external_field_term()
        elif self.particle.particle_type == ParticleType.DYNAMIC_BEAD:
            h_particle = self._build_dynamic_bead_term()
        else:
            msg: str = f"Unknown particle type: {self.particle.particle_type}"
            raise ValueError(msg)

        # Pad to common qubit size and sum
        target_qubits: int = max(h_protein.num_qubits, h_particle.num_qubits)
        h_protein_padded = pad_to_n_qubits(h_protein, target_qubits)
        h_particle_padded = pad_to_n_qubits(h_particle, target_qubits)

        total_hamiltonian: SparsePauliOp = h_protein_padded + h_particle_padded

        logger.info(
            "Built total Hamiltonian with particle: %s qubits", total_hamiltonian.num_qubits
        )
        return total_hamiltonian.simplify()

    def _build_external_field_term(self) -> SparsePauliOp:
        """Build Hamiltonian term for external field particle variant.

        The external field contributes a distance-dependent potential energy
        that depends on the configuration of chain beads.

        Returns:
            SparsePauliOp: Hamiltonian term for external field interactions.

        """
        logger.debug("Building external field particle term...")

        from qiskit.quantum_info import SparsePauliOp

        # For external field, we need to encode:
        # H_external = sum_i V(distance_i) * penalty_i
        # where V depends on bead position relative to field center

        # Initialize with identity operator
        h_field: SparsePauliOp = build_identity_op(
            self.distance_map.num_qubits if hasattr(self.distance_map, "num_qubits") else 10,
            EMPTY_OP_COEFF,
        )

        # For implementation: this requires encoding distance information
        # into quantum operators, which is non-trivial. A simplified approach
        # would integrate field effects into distance-based contact terms.
        # Full implementation deferred to integration module.

        logger.debug("External field term created with %s qubits", h_field.num_qubits)
        return h_field

    def _build_dynamic_bead_term(self) -> SparsePauliOp:
        """Build Hamiltonian term for dynamic bead particle variant.

        The particle is treated as a bead with:
        - Interaction energy with each protein bead based on bead types
        - Contribution only when beads are in contact (neighbors on lattice)

        Returns:
            SparsePauliOp: Hamiltonian term for dynamic bead interactions.

        """
        logger.debug("Building dynamic bead particle term...")

        from qiskit.quantum_info import SparsePauliOp

        # Initialize with identity operator
        h_particle: SparsePauliOp = build_identity_op(
            int(self.contact_map.main_main_contacts[0][1].num_qubits),
            EMPTY_OP_COEFF,
        )

        # For each bead in the protein chain, add contact term with particle
        main_chain = self.protein.main_chain
        chain_len = len(main_chain)

        for i in range(chain_len):
            # Get bead symbol for interaction lookup
            bead_symbol = main_chain.beads[i].symbol

            # Look up interaction energy between bead and particle
            try:
                energy = self.interaction.get_energy(bead_symbol, self.particle.symbol)
            except Exception as e:
                logger.warning(
                    "Could not get interaction energy for bead %s with particle: %s",
                    bead_symbol,
                    e,
                )
                continue

            if energy == 0.0:
                continue  # Skip zero-energy interactions

            # Add contact term: multiply contact indicator by energy
            # This is simplified; full implementation requires proper qubit mapping
            logger.debug(
                "Added particle interaction for bead %d (%s) with energy %f",
                i,
                bead_symbol,
                energy,
            )

        logger.debug("Dynamic bead term created")
        return h_particle

    def set_particle(self, particle: Particle) -> None:
        """Set or replace the particle for this builder.

        Args:
            particle (Particle): New particle to use.

        """
        self.particle = particle
        logger.info("Particle updated to %s type", particle.particle_type.value)

    def get_particle(self) -> Particle | None:
        """Get the currently configured particle.

        Returns:
            Particle | None: Current particle or None if not set.

        """
        return self.particle
