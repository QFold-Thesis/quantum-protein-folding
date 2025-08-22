from pathlib import Path
from constants import MJ_INTERACTION_MATRIX_FILEPATH
import numpy as np


class MJInteraction:
    """Class representing Miyazawa-Jernigan interactions between amino acids."""

    def __init__(
        self,
        protein_sequence: str,
        interaction_matrix_path: Path = MJ_INTERACTION_MATRIX_FILEPATH,
    ) -> None:
        self._protein_sequence: str = protein_sequence
        self._interaction_matrix_path: Path = interaction_matrix_path

        self.valid_symbols: list[str] = []

        self.energy_pairs: dict[str, float] = self._prepare_mj_interaction_matrix(
            self._interaction_matrix_path
        )
        self._check_if_valid_sequence()

    def _prepare_mj_interaction_matrix(
        self, mj_filepath: Path = MJ_INTERACTION_MATRIX_FILEPATH
    ) -> dict[str, float]:
        """Prepare the MJ interaction matrix from the specified file."""

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

        return energy_pairs

    def _check_if_valid_sequence(self) -> bool:
        """Check if the protein sequence contains only valid amino acid symbols."""

        for symbol in self._protein_sequence:
            if symbol not in self.valid_symbols:
                raise ValueError(
                    f"Invalid amino acid symbol '{symbol}' found in the protein sequence."
                )

        return True
