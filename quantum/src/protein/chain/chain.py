from abc import ABC, abstractmethod
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from protein.bead import Bead


class Chain(ABC):
    """
    Abstract base class for protein chains, defining shared behavior for all chain types.
    """

    @abstractmethod
    def __init__(self, protein_sequence: str) -> None:
        """
        Initializes the chain with an empty list of beads.
        """
        self.beads: list[Bead] = []

    def get_symbol_at(self, index: int) -> str:
        """
        Returns the bead symbol at the given chain index.

        """
        return self.beads[index].symbol
