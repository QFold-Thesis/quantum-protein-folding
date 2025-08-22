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


class MJInteraction:
    """
    Represent Miyazawa-Jernigan interactions between amino acids.

    The class reads an MJ interaction matrix and exposes:

    - ``valid_symbols``: the list of supported residue symbols
    - ``energy_pairs``: a dictionary mapping two-letter residue pairs to
      their contact energy (e.g., ``{"AL": -0.6, "LA": -0.6, ...}``).

    Attributes
    ----------
    _protein_sequence:
        The input protein sequence as a string of one-letter residue codes.
    _interaction_matrix_path:
        Path to the MJ interaction matrix file.
    valid_symbols:
        List of one-letter residue symbols parsed from the matrix header.
    energy_pairs:
        Mapping from two-letter residue pairs to contact energies. Keys are
        symmetric (e.g., "AB" and "BA" both exist with the same value).

    """

    def __init__(
        self,
        protein_sequence: str,
        interaction_matrix_path: Path = MJ_INTERACTION_MATRIX_FILEPATH,
    ) -> None:
        """
        Initialize the interaction table and validate the sequence.

        Parameters
        ----------
        protein_sequence:
            Protein sequence consisting of one-letter residue symbols that must
            be present in the interaction matrix.
        interaction_matrix_path:
            Path to the MJ matrix file. Defaults to
            :data:`constants.MJ_INTERACTION_MATRIX_FILEPATH`.

        Raises
        ------
        ValueError
            If ``protein_sequence`` contains a symbol not present in the loaded
            matrix.

        """
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
        """
        Load and symmetrize the MJ interaction matrix from a file.

        The loader expects a header row of residue symbols and the upper
        triangular (including diagonal) numeric values beginning from the
        second row. For each parsed value, two keys are inserted into the
        resulting dictionary so that pair energies are symmetric.

        Parameters
        ----------
        mj_filepath:
            Path to the text file containing the MJ matrix.

        Returns
        -------
        dict[str, float]
            A mapping from two-letter residue pairs (e.g., "AL") to contact
            energies, with symmetric entries present for reversed pairs.

        Notes
        -----
        - ``valid_symbols`` is populated from the matrix header during the
          load.
        - Any I/O or parsing errors raised by NumPy are propagated.

        """
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
        """
        Validate that the sequence contains only supported residue symbols.

        Returns
        -------
        bool
            ``True`` if the sequence is valid.

        Raises
        ------
        ValueError
            If an unknown residue symbol is encountered in the sequence.

        """
        for symbol in self._protein_sequence:
            if symbol not in self.valid_symbols:
                msg = f"Invalid amino acid symbol '{symbol}' found in the protein sequence."
                raise ValueError(msg)

        return True
