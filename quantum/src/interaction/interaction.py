from __future__ import annotations

from abc import ABC, abstractmethod
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from pathlib import Path


class Interaction(ABC):
    """
    Abstract base class for interaction models.
    Subclasses must implement get_energy returning a numeric energy value and __init__.
    """

    @abstractmethod
    def __init__(self, interaction_matrix_path: Path) -> None:
        self.interaction_matrix_path: Path = interaction_matrix_path

    @abstractmethod
    def get_energy(self, *args: Any, **kwargs: Any) -> float:
        """
        Compute and return the interaction energy.
        Must be implemented by subclasses.
        """
        raise NotImplementedError
