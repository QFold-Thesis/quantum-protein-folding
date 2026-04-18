"""Abstract base class for particle models in protein folding simulations.

Defines the interface for particles that can interact with protein chains.
Supports two main variants: external fields and dynamic beads.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from enum import Enum
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from protein import Protein


class ParticleType(Enum):
    """Enumeration of supported particle types."""

    EXTERNAL_FIELD = "external_field"
    DYNAMIC_BEAD = "dynamic_bead"


class Particle(ABC):
    """Abstract base class for particles that interact with protein chains.

    Attributes:
        particle_type (ParticleType): Type of particle (field or dynamic bead).
        interaction_strength (float): Overall strength coefficient for particle interactions.
        symbol (str): Single-character symbol representing the particle for interaction matrices.

    """

    def __init__(
        self,
        particle_type: ParticleType,
        symbol: str = "L",
        interaction_strength: float = 1.0,
    ) -> None:
        """Initialize a particle.

        Args:
            particle_type (ParticleType): Type of this particle.
            symbol (str, optional): Symbol for particle in interaction calculations. Defaults to "L" (ligand).
            interaction_strength (float, optional): Overall interaction strength multiplier.
                Defaults to 1.0.

        """
        self.particle_type: ParticleType = particle_type
        self.symbol: str = symbol
        self.interaction_strength: float = interaction_strength

    @abstractmethod
    def get_energy_contribution(self, protein: Protein, configuration: dict) -> float:
        """Calculate the energy contribution of this particle to the system.

        Args:
            protein (Protein): Protein object with current folding configuration.
            configuration (dict): Dictionary containing coordinate information and folding state.
                Structure depends on particle type.

        Returns:
            float: Energy contribution from particle interactions, scaled by interaction_strength.

        """
        pass

    @abstractmethod
    def validate_compatibility(self, protein: Protein) -> bool:
        """Validate that this particle configuration is compatible with the protein.

        Args:
            protein (Protein): Protein to check compatibility with.

        Returns:
            bool: True if particle can interact with this protein, False otherwise.

        """
        pass
