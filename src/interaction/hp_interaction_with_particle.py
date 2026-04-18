"""Extended HP interaction model with particle support.

Extends the classical HP model to include interactions between chain beads and
an external particle (ligand). The particle can interact attractively or repulsively
with hydrophobic, polar, or both types of residues.
"""

from __future__ import annotations

from pathlib import Path

from constants import (
    HP_HH_CONTACT_ENERGY,
    HP_INTERACTION_MATRIX_FILEPATH,
    HP_NON_HH_CONTACT_ENERGY,
)
from exceptions import UnsupportedAminoAcidSymbolError
from interaction.hp_interaction import HPInteraction
from logger import get_logger

logger = get_logger()


class HPInteractionWithParticle(HPInteraction):
    """HP interaction model extended to include particle-chain interactions.

    Adds a particle type symbol with configurable interactions to hydrophobic
    and/or polar residues. Maintains backward compatibility with standard HP calculations.

    Attributes:
        particle_symbol (str): Single-character symbol for the particle (e.g., "L" for ligand).
        particle_hh_energy (float): Energy for particle-hydrophobic contacts.
        particle_pp_energy (float): Energy for particle-polar contacts.

    """

    def __init__(
        self,
        interaction_matrix_path: Path = HP_INTERACTION_MATRIX_FILEPATH,
        particle_symbol: str = "L",
        particle_hh_energy: float = -1.0,
        particle_pp_energy: float = 0.0,
    ) -> None:
        """Initialize HP model with particle support.

        Args:
            interaction_matrix_path (Path): Path to HP matrix file.
            particle_symbol (str, optional): Particle symbol. Defaults to "L".
            particle_hh_energy (float, optional): Energy for particle-H contacts. Defaults to -1.0.
            particle_pp_energy (float, optional): Energy for particle-P contacts. Defaults to 0.0.

        """
        super().__init__(interaction_matrix_path)

        self.particle_symbol: str = particle_symbol
        self.particle_hh_energy: float = particle_hh_energy
        self.particle_pp_energy: float = particle_pp_energy

        # Add particle to valid symbols
        self.valid_symbols.add(particle_symbol)

        logger.info(
            "HPInteractionWithParticle initialized with particle symbol '%s'",
            particle_symbol,
        )

    def get_energy(self, symbol_i: str, symbol_j: str) -> float:
        """Get interaction energy between two symbols (extended to include particle).

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

        # Fall back to standard HP energy for residue-residue interactions
        return super().get_energy(symbol_i, symbol_j)

    def _get_particle_energy(self, symbol_i: str, symbol_j: str) -> float:
        """Calculate energy when particle is involved.

        Args:
            symbol_i (str): First symbol.
            symbol_j (str): Second symbol.

        Returns:
            float: Interaction energy.

        Raises:
            UnsupportedAminoAcidSymbolError: If non-particle symbol not recognized.

        """
        # Determine which symbol is the particle
        if symbol_i == self.particle_symbol:
            amino_acid = symbol_j
        else:
            amino_acid = symbol_i

        if amino_acid not in self.valid_symbols:
            msg = f"Unsupported amino acid symbol: {amino_acid}"
            raise UnsupportedAminoAcidSymbolError(msg)

        # Return energy based on whether amino acid is hydrophobic or polar
        if amino_acid in self._hydrophobic_symbols:
            return self.particle_hh_energy
        else:
            return self.particle_pp_energy

    def set_particle_interaction(
        self,
        particle_hh_energy: float | None = None,
        particle_pp_energy: float | None = None,
    ) -> None:
        """Update particle interaction energies.

        Args:
            particle_hh_energy (float, optional): New energy for particle-H. If None, keep current.
            particle_pp_energy (float, optional): New energy for particle-P. If None, keep current.

        """
        if particle_hh_energy is not None:
            self.particle_hh_energy = particle_hh_energy

        if particle_pp_energy is not None:
            self.particle_pp_energy = particle_pp_energy

        logger.debug(
            "Updated particle interactions: H=%f, P=%f",
            self.particle_hh_energy,
            self.particle_pp_energy,
        )
