"""
Miyazawa-Jernigan (MJ) pairwise contact potentials.

This module provides a small utility class, ``MJInteraction``, that:

- loads an MJ interaction matrix from a plain-text file,
- records the set of valid residue symbols (one-letter codes), and
- builds a dictionary that maps residue-residue pairs (e.g., "AL", "LA")
  to their contact energies.

Matrix format assumptions
-------------------------
The file is expected to contain a header row with one-letter residue symbols.
Starting from the second row, each row begins with a residue symbol and the
upper triangle (including the diagonal) of energies is provided. The loader
mirrors values to create a symmetric mapping, so both "AB" and "BA" keys are
present with the same energy.
"""

from pathlib import Path

import numpy as np

from constants import MJ_INTERACTION_MATRIX_FILEPATH
from exceptions import UnsupportedAminoAcidSymbolError
from interaction import Interaction
from logger import get_logger

logger = get_logger()


class MJInteraction(Interaction):
    def __init__(
        self,
        interaction_matrix_path: Path = MJ_INTERACTION_MATRIX_FILEPATH,
    ) -> None:
        """
        Initialize MJInteraction instance.

        Loads the MJ interaction matrix from a file and prepares a mapping of
        residue-residue pairs to their contact energies.

        Args:
            interaction_matrix_path (Path): Path to the MJ interaction matrix file.

        """
        super().__init__(interaction_matrix_path)
        logger.debug("Initializing MJInteraction...")

        self._energy_pairs: dict[str, float] = self._prepare_mj_interaction_matrix(
            self._interaction_matrix_path
        )

        self.valid_symbols = {symbol for pair in self._energy_pairs for symbol in pair}

        logger.info(
            f"MJInteraction initialized with {len(self.valid_symbols)} valid amino acid symbols"
        )

    def _prepare_mj_interaction_matrix(
        self, mj_filepath: Path = MJ_INTERACTION_MATRIX_FILEPATH
    ) -> dict[str, float]:
        """
        Prepare the MJ interaction matrix.

        Reads the MJ matrix file, records valid residue symbols, and builds a
        dictionary of symmetric residue-residue contact energies.

        Args:
            mj_filepath (Path): Path to the MJ matrix file.

        Returns:
            dict[str, float]: Dictionary mapping residue pair codes to energies.

        Raises:
            Exception: If an error occurs while loading or parsing the matrix.

        """
        try:
            mj_matrix = np.loadtxt(mj_filepath, dtype=str)

            header = mj_matrix[0, :]
            energy_pairs: dict[str, float] = {}

            for row in range(1, np.shape(mj_matrix)[0]):
                row_symbol = header[row - 1]
                for col in range(row - 1, np.shape(mj_matrix)[1]):
                    col_symbol = header[col]
                    energy = float(mj_matrix[row, col])
                    key = f"{col_symbol}{row_symbol}"
                    energy_pairs[key] = energy
                    energy_pairs[key[::-1]] = energy
        except Exception:
            logger.exception("Error loading MJ matrix")
            raise
        else:
            logger.info(
                f"Successfully loaded {len(energy_pairs)} energy pairs from MJ matrix at: {mj_filepath}"
            )
            return energy_pairs

    def get_energy(self, symbol_i: str, symbol_j: str) -> float:
        """
        Return MJ interaction energy for a pair of residue symbols.

        Args:
            symbol_i (str): One-letter code of the first residue.
            symbol_j (str): One-letter code of the second residue.

        Returns:
            float: MJ interaction energy between the two residues.

        Raises:
            UnsupportedAminoAcidSymbolError: If either residue symbol is not in the MJ matrix.

        """
        key = f"{symbol_i}{symbol_j}"
        try:
            return self._energy_pairs[key]
        except KeyError as e:
            msg: str = f"Energy pair of '{key}' not supported in loaded MJ interaction model"
            raise UnsupportedAminoAcidSymbolError(msg) from e
