"""External field particle implementation.

In this variant, the particle is modeled as a static external potential field
that affects the energy of the protein based on distance from field center.
The particle does not have mobility - it represents a fixed external potential.

Mathematical formulation:
    H_external = sum_i V(r_i - r_particle) * c_i
where:
    - V(r) is the potential energy function
    - r_i is position of bead i
    - r_particle is fixed particle location
    - c_i characterizes bead i (hydrophobic/polar in HP, amino acid type in MJ)
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from logger import get_logger
from particle.particle import Particle, ParticleType

if TYPE_CHECKING:
    from protein import Protein

logger = get_logger()


class ExternalFieldParticle(Particle):
    """Particle implemented as static external field affecting chain conformation.

    In this model, the particle is placed at a fixed position and exerts a
    distance-dependent potential on chain beads. The potential can be:
    - Distance-dependent (e.g., Lennard-Jones like)
    - Based on bead type (hydrophobic attraction, etc.)
    - Uniform across all bead types

    Attributes:
        position (tuple[int, int, int]): Fixed 3D lattice coordinates (x, y, z).
        potential_type (str): Type of potential function ('distance', 'depth_dependent', 'well').
        well_depth (float): Depth of potential well (energy scale).
        well_radius (float): Characteristic range of potential (lattice units).

    """

    def __init__(
        self,
        position: tuple[int, int, int],
        potential_type: str = "distance",
        well_depth: float = -1.0,
        well_radius: float = 2.0,
        symbol: str = "L",
        interaction_strength: float = 1.0,
    ) -> None:
        """Initialize external field particle.

        Args:
            position (tuple[int, int, int]): Fixed 3D position in lattice (x, y, z).
            potential_type (str, optional): Shape of potential.
                Options: 'distance', 'depth_dependent', 'well'. Defaults to "distance".
            well_depth (float, optional): Depth of potential well. Defaults to -1.0.
            well_radius (float, optional): Characteristic range in lattice units.
                Defaults to 2.0.
            symbol (str, optional): Particle symbol for interactions. Defaults to "L".
            interaction_strength (float, optional): Overall strength multiplier.
                Defaults to 1.0.

        """
        super().__init__(
            particle_type=ParticleType.EXTERNAL_FIELD,
            symbol=symbol,
            interaction_strength=interaction_strength,
        )
        self.position: tuple[int, int, int] = position
        self.potential_type: str = potential_type
        self.well_depth: float = well_depth
        self.well_radius: float = well_radius
        logger.debug(
            "ExternalFieldParticle initialized at position %s with %s potential",
            self.position,
            self.potential_type,
        )

    def get_energy_contribution(self, protein: Protein, configuration: dict) -> float:
        """Calculate energy of protein in external field.

        Args:
            protein (Protein): Protein object.
            configuration (dict): Dictionary containing 'coordinates' mapping bead indices
                to their (x, y, z) positions.

        Returns:
            float: Total energy contribution from external field interactions.

        """
        if "coordinates" not in configuration:
            msg = "Configuration must contain 'coordinates' mapping bead index to (x,y,z) tuple"
            raise KeyError(msg)

        coordinates = configuration["coordinates"]
        total_energy = 0.0

        # Sum potential energy over all beads
        for bead_idx, bead_coord in coordinates.items():
            distance = self._euclidean_distance(bead_coord, self.position)
            bead_energy = self._evaluate_potential(distance)
            total_energy += bead_energy

        return total_energy * self.interaction_strength

    def _euclidean_distance(
        self, coord1: tuple[float, float, float], coord2: tuple[float, float, float]
    ) -> float:
        """Calculate Euclidean distance between two 3D coordinates."""
        return sum((c1 - c2) ** 2 for c1, c2 in zip(coord1, coord2)) ** 0.5

    def _evaluate_potential(self, distance: float) -> float:
        """Evaluate potential energy as function of distance.

        Args:
            distance (float): Distance from particle.

        Returns:
            float: Potential energy at this distance.

        """
        if self.potential_type == "distance":
            # Simple distance-dependent: V(r) = depth / (1 + r)
            return self.well_depth / (1.0 + distance)

        elif self.potential_type == "depth_dependent":
            # Gaussian-like: V(r) = depth * exp(-(r/radius)^2)
            exponent = -(distance / self.well_radius) ** 2
            return self.well_depth * (2.71828 ** exponent)  # e^x approximation

        elif self.potential_type == "well":
            # Square well: V(r) = depth if r < radius else 0
            if distance < self.well_radius:
                return self.well_depth
            return 0.0

        else:
            msg = f"Unknown potential type: {self.potential_type}"
            raise ValueError(msg)

    def validate_compatibility(self, protein: Protein) -> bool:
        """External field particles are compatible with all proteins.

        Args:
            protein (Protein): Protein to check.

        Returns:
            bool: Always True for external field particles.

        """
        return True
