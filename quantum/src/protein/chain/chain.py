from abc import ABC, abstractmethod
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from protein.bead import Bead


class Chain(ABC):
    """Abstract base class for protein chains, defining shared behavior for all chain types."""

    @abstractmethod
    def __init__(self, protein_sequence: str) -> None:
        """
        Initialize the chain with an empty list of beads.

        Args:
            protein_sequence (str): The amino acid sequence representing the protein chain.

        """
        self.beads: list[Bead] = []

    def get_symbol_at(self, index: int) -> str:
        """
        Return the symbol of the bead at the given chain index.

        Args:
            index (int): Index of the bead in the chain.

        Returns:
            str: The symbol of the bead located at the specified index.

        """
        return self.beads[index].symbol
