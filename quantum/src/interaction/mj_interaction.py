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

from logger import get_logger
from exceptions import InvalidAminoAcidError
from protein import Protein

logger = get_logger()


class MJInteraction:

    def __init__(
        self,
        protein: Protein,
        interaction_matrix_path: Path = MJ_INTERACTION_MATRIX_FILEPATH,
    ) -> None:

        self._interaction_matrix_path: Path = interaction_matrix_path

        self.valid_symbols: list[str] = []

        self.energy_pairs: dict[str, float] = self._prepare_mj_interaction_matrix(
            self._interaction_matrix_path
        )

        protein_sequence: str = str(protein.main_chain) + str(protein.side_chain)
        
        self._check_if_valid_sequence(protein_sequence)

    def _prepare_mj_interaction_matrix(
        self, mj_filepath: Path = MJ_INTERACTION_MATRIX_FILEPATH
    ) -> dict[str, float]:
        try:
            mj_matrix = np.loadtxt(mj_filepath, dtype=str)


            self.valid_symbols: list[str] = [str(symbol) for symbol in mj_matrix[0, :]]

            energy_pairs: dict[str, float] = {}

            for row in range(1, np.shape(mj_matrix)[0]):
                for col in range(row - 1, np.shape(mj_matrix)[1]):
                    bead_1: str = self.valid_symbols[col]
                    bead_2: str = self.valid_symbols[row - 1]
                    energy: float = float(mj_matrix[row, col])

                    key: str = f"{bead_1}{bead_2}"

                    energy_pairs[key] = energy
                    energy_pairs[key[::-1]] = energy
        except Exception as e:
            logger.error(f"Error loading MJ matrix: {e}")
            raise e
        else:
            logger.debug(f"Successfully loaded {len(energy_pairs)} energy pairs from MJ matrix at: {mj_filepath}")
            return energy_pairs


    def _check_if_valid_sequence(self, protein_sequence: str) -> bool:
        for symbol in protein_sequence:
            if symbol not in self.valid_symbols:
                msg = f"Invalid amino acid symbol '{symbol}' found in the protein sequence."
                logger.error(msg)
                raise InvalidAminoAcidError(msg)
        return True
