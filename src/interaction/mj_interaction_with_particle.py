"""Extended MJ interaction model with particle support.

Extends the Miyazawa-Jernigan model to include a particle as an additional
interaction partner. The particle is treated as a pseudo-amino-acid with
configurable energies against all standard residues.
"""

from __future__ import annotations

from pathlib import Path

from constants import MJ_INTERACTION_MATRIX_FILEPATH
from exceptions import UnsupportedAminoAcidSymbolError
from interaction.mj_interaction import MJInteraction
from logger import get_logger

logger = get_logger()


class MJInteractionWithParticle(MJInteraction):
    """MJ interaction model extended with particle support.

    Adds a particle type as a pseudo-amino-acid with its own row/column
    in the interaction matrix. The particle can have distinct energies
    with each protein residue type.

    Attributes:
        particle_symbol (str): Single-character symbol for the particle (typically "L").
        particle_interactions (dict[str, float]): Mapping of residue symbols to
            interaction energies with the particle.

    """

    def __init__(
        self,
        interaction_matrix_path: Path = MJ_INTERACTION_MATRIX_FILEPATH,
        particle_symbol: str = "L",
    ) -> None:
        """Initialize MJ model with particle support.

        Args:
            interaction_matrix_path (Path): Path to MJ interaction matrix file.
            particle_symbol (str, optional): Symbol for particle. Defaults to "L".

        """
        super().__init__(interaction_matrix_path)

        self.particle_symbol: str = particle_symbol
        self.particle_interactions: dict[str, float] = {}

        # Add particle to valid symbols
        self.valid_symbols.add(particle_symbol)

        logger.info(
            "MJInteractionWithParticle initialized with particle symbol '%s'",
            particle_symbol,
        )

    def get_energy(self, symbol_i: str, symbol_j: str) -> float:
        """Get interaction energy (extended to include particle).

        Args:
            symbol_i (str): First symbol (amino acid or particle).
            symbol_j (str): Second symbol (amino acid or particle).

        Returns:
            float: Interaction energy.

        Raises:
            UnsupportedAminoAcidSymbolError: If symbols not recognized.

        """
        # Check if either symbol is the particle
        if symbol_i == self.particle_symbol or symbol_j == self.particle_symbol:
            return self._get_particle_energy(symbol_i, symbol_j)

        # Fall back to standard MJ energy for residue-residue interactions
        return super().get_energy(symbol_i, symbol_j)

    def _get_particle_energy(self, symbol_i: str, symbol_j: str) -> float:
        """Calculate energy when particle is involved.

        Args:
            symbol_i (str): First symbol.
            symbol_j (str): Second symbol.

        Returns:
            float: Interaction energy.

        Raises:
            UnsupportedAminoAcidSymbolError: If residue symbol not recognized.

        """
        # Determine which symbol is the particle
        if symbol_i == self.particle_symbol:
            amino_acid = symbol_j
        else:
            amino_acid = symbol_i

        if amino_acid not in self.valid_symbols:
            msg = f"Unsupported amino acid symbol: {amino_acid}"
            raise UnsupportedAminoAcidSymbolError(msg)

        # Return particle interaction energy if set, otherwise 0
        if amino_acid in self.particle_interactions:
            return self.particle_interactions[amino_acid]

        return 0.0  # Default: no interaction if not explicitly set

    def set_particle_interactions(self, interaction_dict: dict[str, float]) -> None:
        """Set particle-residue interaction energies.

        Args:
            interaction_dict (dict[str, float]): Mapping of residue symbols to
                interaction energies with particle.

        """
        self.particle_interactions = interaction_dict.copy()
        logger.debug(
            "Set particle interactions for %d residue types", len(interaction_dict)
        )

    def set_particle_interaction_for_residue(
        self, residue_symbol: str, energy: float
    ) -> None:
        """Set interaction energy for a single residue type.

        Args:
            residue_symbol (str): Single-character residue symbol.
            energy (float): Interaction energy with particle.

        Raises:
            UnsupportedAminoAcidSymbolError: If residue symbol not in MJ matrix.

        """
        if residue_symbol not in self.valid_symbols or residue_symbol == self.particle_symbol:
            msg = f"Cannot set particle interaction for symbol: {residue_symbol}"
            raise UnsupportedAminoAcidSymbolError(msg)

        self.particle_interactions[residue_symbol] = energy
        logger.debug(
            "Set particle-%s interaction energy to %f", residue_symbol, energy
        )

    def get_particle_interaction_column(self) -> dict[str, float]:
        """Get the complete particle interaction column.

        Returns:
            dict[str, float]: Current particle-residue interactions.

        """
        return self.particle_interactions.copy()
