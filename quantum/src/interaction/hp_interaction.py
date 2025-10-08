"""
Hydrophobic / Polar (HP) interaction model.

This module provides ``HPInteraction`` which:

- loads a simple text file mapping residue one-letter symbols to a binary
    hydrophobic(1)/polar(0) flag,
- exposes the list of hydrophobic symbols, and
- provides a pairwise energy function where only hydrophobic-hydrophobic
    contacts contribute (default energy -1.0) while any pair involving a polar
    residue contributes 0.0.

HP matrix file format
---------------------
Expected file (e.g. ``hp_matrix.txt``) structure:

        A 1
        R 0
        N 0
        ...

Columns are separated by arbitrary whitespace. Only the first two tokens on a line are considered:
the residue symbol (single character) and a flag (``0`` or ``1``).

Symmetry / Energies
-------------------
The classical lattice HP model assigns the same energy to any hydrophobic-
hydrophobic contact independent of residue identities (here set to -1.0).
All other pairs have energy 0.0.
"""

from pathlib import Path

import numpy as np

from constants import (
    HP_HH_CONTACT_ENERGY,
    HP_INTERACTION_MATRIX_FILEPATH,
    HP_NON_HH_CONTACT_ENERGY,
)
from interaction import Interaction
from logger import get_logger

logger = get_logger()


class HPInteraction(Interaction):
    def __init__(
        self,
        interaction_matrix_path: Path = HP_INTERACTION_MATRIX_FILEPATH,
    ) -> None:
        super().__init__(interaction_matrix_path)
        self._hydrophobic_symbols: list[str] = self._load_hydrophobic_symbols(
            self._interaction_matrix_path
        )

    def _load_hydrophobic_symbols(
        self, hp_filepath: Path = HP_INTERACTION_MATRIX_FILEPATH
    ) -> list[str]:
        """
        Load hydrophobic amino acid symbols from a HP interaction matrix file.

        Parses a small, strict matrix file where each line has the format '<SYMBOL> <0|1>'.
        Ignores blank lines and lines starting with '#'. Only symbols with value '1' are returned.

        Args:
            hp_filepath (Path, optional): Path to the HP interaction matrix file.
                Defaults to HP_INTERACTION_MATRIX_FILEPATH.

        Returns:
            list[str]: List of hydrophobic amino acid symbols.

        Raises:
            Exception: If the HP matrix file cannot be read or parsed.

        """
        hydros: list[str] = []
        try:
            hp_matrix = np.loadtxt(hp_filepath, dtype=str)
            for row in range(1, np.shape(hp_matrix)[0]):
                for col in range(1, np.shape(hp_matrix)[1]):
                    if hp_matrix[row, col] == "1":
                        hydros.extend(hp_matrix[0, col])
        except Exception:
            logger.exception("Error loading HP matrix")
            raise
        else:
            logger.debug(
                f"Successfully loaded {len(hydros)} hydrophobic symbols from HP matrix at: {hp_filepath}"
            )
            return hydros

    def _is_hydrophobic(self, symbol: str) -> bool:
        """
        Check if an amino acid symbol is hydrophobic.

        Args:
            symbol (str): Single-letter amino acid symbol.

        Returns:
            bool: True if the symbol is hydrophobic, False otherwise.

        """
        return symbol in self._hydrophobic_symbols

    def get_energy(self, symbol_i: str, symbol_j: str) -> float:
        """
        Return the HP model pair energy.

        Returns the hydrophobic-hydrophobic contact energy (HP_HH_CONTACT_ENERGY) if both residues are hydrophobic,
        otherwise returns the non-hydrophobic contact energy (HP_NON_HH_CONTACT_ENERGY).

        Args:
            symbol_i (str): Single-letter amino acid symbol for the first residue.
            symbol_j (str): Single-letter amino acid symbol for the second residue.

        Returns:
            float: Interaction energy according to the HP model.

        Raises:
            RuntimeError: If an error occurs while computing the pair energy.

        """
        try:
            return (
                HP_HH_CONTACT_ENERGY
                if (self._is_hydrophobic(symbol_i) and self._is_hydrophobic(symbol_j))
                else HP_NON_HH_CONTACT_ENERGY
            )
        except Exception as e:
            msg: str = f"Error computing HP energy for pair: {symbol_i}, {symbol_j}"
            logger.exception(msg)
            raise RuntimeError(msg) from e
