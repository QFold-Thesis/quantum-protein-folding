"""Dynamic bead particle implementation.

In this variant, the particle is modeled as a full bead in the lattice/diamond
representation with its own coordinates, allowed moves, and pairwise interactions
with chain beads following an extended potential matrix.

Mathematical formulation:
    H_particle = sum_i V(type_i, type_particle) * contact(r_i, r_particle)
where:
    - V is extended interaction matrix including particle type
    - type_i and type_particle are bead types/symbols
    - contact(r_i, r_particle) is 1 if beads are neighbors, 0 otherwise
    - Particle coordinates constrained to valid lattice positions
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from logger import get_logger
from particle.particle import Particle, ParticleType

if TYPE_CHECKING:
    from interaction.interaction import Interaction
    from protein import Protein

logger = get_logger()


class DynamicBeadParticle(Particle):
    """Particle implemented as full bead with lattice coordinates and mobility.

    The particle is treated as a legitimate lattice bead that:
    - Has coordinates on the same lattice as chain beads
    - Interacts with chain beads according to extended interaction matrix
    - Can sample different positions (in optimization context)
    - Shares the same allowed move set as chain beads

    Attributes:
        position (tuple[int, int, int]): Current 3D lattice coordinates (x, y, z).
        fixed (bool): If True, particle position is fixed; if False, it can move.
        interaction_matrix_column (dict[str, float]): Extended interaction potentials
            mapping amino acid symbols to energies for this particle type.

    """

    def __init__(
        self,
        position: tuple[int, int, int],
        symbol: str = "L",
        fixed: bool = False,
        interaction_strength: float = 1.0,
    ) -> None:
        """Initialize dynamic bead particle.

        Args:
            position (tuple[int, int, int]): Initial lattice position (x, y, z).
            symbol (str, optional): Single-character symbol for this particle. Defaults to "L".
            fixed (bool, optional): If False, particle can move during folding.
                If True, maintains fixed position. Defaults to False.
            interaction_strength (float, optional): Overall strength multiplier
                for particle interactions. Defaults to 1.0.

        """
        super().__init__(
            particle_type=ParticleType.DYNAMIC_BEAD,
            symbol=symbol,
            interaction_strength=interaction_strength,
        )
        self.position: tuple[int, int, int] = position
        self.fixed: bool = fixed
        self.interaction_matrix_column: dict[str, float] = {}
        logger.debug(
            "DynamicBeadParticle initialized at position %s (fixed=%s)",
            self.position,
            self.fixed,
        )

    def set_interaction_potentials(self, potentials: dict[str, float]) -> None:
        """Set interaction energies between particle and each amino acid type.

        Args:
            potentials (dict[str, float]): Mapping of amino acid symbols to
                interaction energies with this particle.

        """
        self.interaction_matrix_column = potentials.copy()
        logger.debug(
            "Set interaction potentials for particle with %d amino acid types",
            len(potentials),
        )

    def get_energy_contribution(self, protein: Protein, configuration: dict) -> float:
        """Calculate energy from particle interactions with chain beads.

        Counts contacts between particle and chain beads, applying interaction
        energies based on amino acid types and spatial proximity.

        Args:
            protein (Protein): Protein object with bead type information.
            configuration (dict): Dictionary containing:
                - 'coordinates': mapping bead index to (x, y, z) tuple
                - 'particle_position': (x, y, z) tuple for particle location
                - 'contact_distance': Lattice distance threshold for contacts

        Returns:
            float: Total energy from particle-chain interactions.

        Raises:
            KeyError: If required fields missing from configuration.
            ValueError: If particle position or interaction matrix not set.

        """
        if "coordinates" not in configuration:
            msg = "Configuration must contain 'coordinates' mapping bead index to (x,y,z)"
            raise KeyError(msg)

        if "particle_position" not in configuration:
            msg = "Configuration must contain 'particle_position' for dynamic bead"
            raise KeyError(msg)

        if not self.interaction_matrix_column:
            msg = "Particle interaction potentials not set"
            raise ValueError(msg)

        coordinates = configuration["coordinates"]
        particle_pos = configuration["particle_position"]
        contact_distance = configuration.get("contact_distance", 1)  # Default: nearest neighbors

        total_energy = 0.0

        # Sum interaction energy over all contacts
        for bead_idx, bead_coord in coordinates.items():
            if self._are_neighbors(bead_coord, particle_pos, contact_distance):
                # Get bead type (symbol)
                bead_symbol = self._get_bead_symbol(protein, bead_idx)

                # Look up interaction energy
                if bead_symbol in self.interaction_matrix_column:
                    energy = self.interaction_matrix_column[bead_symbol]
                    total_energy += energy

        return total_energy * self.interaction_strength

    def _are_neighbors(
        self,
        coord1: tuple[float, float, float],
        coord2: tuple[float, float, float],
        distance_threshold: float,
    ) -> bool:
        """Check if two beads are neighbors (within contact distance).

        Uses Manhattan distance on lattice (sum of absolute differences).

        Args:
            coord1: First bead coordinate.
            coord2: Second bead coordinate.
            distance_threshold: Maximum distance for contact.

        Returns:
            bool: True if beads are in direct contact.

        """
        manhattan_distance = sum(abs(c1 - c2) for c1, c2 in zip(coord1, coord2))
        return manhattan_distance == distance_threshold

    def _get_bead_symbol(self, protein: Protein, bead_idx: int) -> str:
        """Get the symbol (amino acid type) of a bead at given index.

        Args:
            protein (Protein): Protein object.
            bead_idx (int): Index of bead in chain.

        Returns:
            str: Single-character amino acid symbol.

        """
        if bead_idx < len(protein.main_chain.beads):
            return protein.main_chain.beads[bead_idx].symbol
        else:
            msg = f"Bead index {bead_idx} out of range for protein with {len(protein.main_chain.beads)} beads"
            raise IndexError(msg)

    def update_position(self, new_position: tuple[int, int, int]) -> bool:
        """Update particle position (if not fixed).

        Args:
            new_position (tuple[int, int, int]): New lattice coordinates.

        Returns:
            bool: True if position was updated, False if particle is fixed.

        """
        if self.fixed:
            logger.debug("Cannot update position: particle is fixed")
            return False

        self.position = new_position
        logger.debug("Particle position updated to %s", new_position)
        return True

    def validate_compatibility(self, protein: Protein) -> bool:
        """Validate particle compatibility with protein.

        Currently checks that interaction potentials are defined for all
        amino acids in the protein.

        Args:
            protein (Protein): Protein to validate.

        Returns:
            bool: True if particle can properly interact with all protein beads.

        """
        if not self.interaction_matrix_column:
            logger.warning("Particle interaction potentials not yet set")
            return False

        protein_symbols = set(protein.main_chain.get_sequence())
        missing_symbols = protein_symbols - set(self.interaction_matrix_column.keys())

        if missing_symbols:
            logger.warning(
                "Particle interaction matrix missing symbols: %s", missing_symbols
            )
            return False

        return True
