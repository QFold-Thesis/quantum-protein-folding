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
from exceptions import UnsupportedAminoAcidSymbolError
from interaction import Interaction
from logger import get_logger

logger = get_logger()


class HPInteraction(Interaction):
    def __init__(
        self,
        interaction_matrix_path: Path = HP_INTERACTION_MATRIX_FILEPATH,
    ) -> None:
        super().__init__(interaction_matrix_path)
        self._hydrophobic_symbols, _polar_symbols = self._load_hp_symbols(
            self._interaction_matrix_path
        )
        self.valid_symbols = set(self._hydrophobic_symbols) | set(_polar_symbols)


        logger.debug(
            f"HPInteraction initialized with {len(self.valid_symbols)} valid amino acid symbols."
        )

    def _load_hp_symbols(
        self, hp_filepath: Path = HP_INTERACTION_MATRIX_FILEPATH
    ) -> tuple[list[str], list[str]]:
        """
        Load hydrophobic amino acid symbols from a HP interaction matrix file.

        Parses a small, strict matrix file where each line has the format '<SYMBOL> <0|1>'.
        Ignores blank lines and lines starting with '#'. Lines with '1' indicate hydrophobic residues.

        Args:
            hp_filepath (Path, optional): Path to the HP interaction matrix file.
                Defaults to HP_INTERACTION_MATRIX_FILEPATH.

        Returns:
            tuple[list[str], list[str]]: Tuple containing lists of hydrophobic and polar amino acid symbols.

        Raises:
            Exception: If the HP matrix file cannot be read or parsed.

        """
        hydrophobic: list[str] = []
        polar: list[str] = []
        try:
            hp_matrix = np.loadtxt(hp_filepath, dtype=str)

            for line in hp_matrix:
                if line[1] == "1":
                    hydrophobic.extend(line[0])
                else:
                    polar.extend(line[0])
        except Exception:
            logger.exception("Error loading HP matrix")
            raise
        else:
            logger.debug(
                f"Successfully loaded {len(hydrophobic)} hydrophobic and {len(polar)} polar symbols from HP matrix at: {hp_filepath}"
            )
            return hydrophobic, polar

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
            UnsupportedAminoAcidSymbolError: If either residue symbol is not in the HP matrix.

        """
        if symbol_i not in self.valid_symbols or symbol_j not in self.valid_symbols:
            msg: str = f"Amino acid symbols of {symbol_i}, {symbol_j} not supported in loaded HP interaction model."
            logger.error(msg)
            raise UnsupportedAminoAcidSymbolError(msg)

        return (
            HP_HH_CONTACT_ENERGY
            if (self._is_hydrophobic(symbol_i) and self._is_hydrophobic(symbol_j))
            else HP_NON_HH_CONTACT_ENERGY
        )
