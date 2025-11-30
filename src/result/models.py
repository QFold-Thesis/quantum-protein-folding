"""Data models for interpreting and visualizing protein folding results."""

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class SparseVQEOutput:
    """Data class for storing sparse VQE output results.

    Attributes:
        bitstring (str): The bitstring representation of the quantum state.
        probability (float): The probability associated with the bitstring.
        state (str): The quantum state in string format.
        energy_value (np.complex128): The energy value corresponding to the quantum state.

    """

    bitstring: str
    probability: float
    state: str
    energy_value: np.complex128

    def __repr__(self) -> str:
        """Returns a string representation of the SparseVQEOutput instance.

        Returns:
            str: Concatenated sequence of bead symbols.

        """
        return (
            f"SparseVQEOutput(bitstring={self.bitstring},\n"
            f"  probability={self.probability},\n"
            f"  state={self.state},\n"
            f"  energy_value={self.energy_value})"
        )


@dataclass(frozen=True)
class BeadPosition:
    """Data class for storing the position of a bead in 3D space.

    Attributes:
        index (int): The index of the bead.
        symbol (str): The chemical symbol of the bead.
        x (float): The x-coordinate of the bead.
        y (float): The y-coordinate of the bead.
        z (float): The z-coordinate of the bead.

    """

    index: int
    symbol: str
    x: float
    y: float
    z: float
