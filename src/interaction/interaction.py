"""
Interaction models for protein folding.

Defines the abstract base class `Interaction`, which loads an interaction matrix
and computes interaction energies used in folding models such as HP or MJ.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from pathlib import Path


class Interaction(ABC):
    """
    Abstract base class for interaction models.

    Subclasses must implement `get_energy` to return a numeric energy value and define their own initialization logic.
    Subclasses should also manage to properly set the `valid_symbols` set to determine which amino acid symbols they support.

    Attributes:
        valid_symbols (set[str]): Set of valid amino acid symbols for the interaction model.

    """

    @abstractmethod
    def __init__(self, interaction_matrix_path: Path) -> None:
        """
        Initialize the interaction model.

        Args:
            interaction_matrix_path (Path): Path to the file containing the interaction matrix.

        """
        self._interaction_matrix_path: Path = interaction_matrix_path
        self.valid_symbols: set[str] = set()

    @abstractmethod
    def get_energy(self, *args: Any, **kwargs: Any) -> float:
        """
        Compute and return the interaction energy.

        This method must be implemented by all subclasses.

        Returns:
            float: Calculated interaction energy.

        """
        raise NotImplementedError
