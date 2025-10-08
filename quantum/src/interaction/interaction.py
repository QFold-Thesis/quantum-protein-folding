from __future__ import annotations

from abc import ABC, abstractmethod
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from pathlib import Path


class Interaction(ABC):
    """
    Abstract base class for interaction models.

    Subclasses must implement :meth:`get_energy` to return a numeric energy value
    and define their own initialization logic.
    """

    @abstractmethod
    def __init__(self, interaction_matrix_path: Path) -> None:
        """
        Initialize the interaction model.

        Args:
            interaction_matrix_path (Path): Path to the file containing the interaction matrix.

        """
        self._interaction_matrix_path: Path = interaction_matrix_path

    @abstractmethod
    def get_energy(self, *args: Any, **kwargs: Any) -> float:
        """
        Compute and return the interaction energy.

        This method must be implemented by all subclasses.

        Returns:
            float: Calculated interaction energy.

        """
        raise NotImplementedError
