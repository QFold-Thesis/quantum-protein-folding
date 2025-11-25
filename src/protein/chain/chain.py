"""
Base abstraction for protein chains.

Provides the abstract `Chain` class, representing a sequence of `Bead` objects
and defining common operations such as indexing, iteration, and string
conversion.
"""

from abc import ABC, abstractmethod
from collections.abc import Iterator

from protein.bead import Bead


class Chain(ABC):
    """
    Abstract base class for protein chains, defining shared behavior for all chain types.
    
    Attributes:
        beads (list[Bead]): List of beads comprising the protein chain.
    """

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

    def __iter__(self) -> Iterator[Bead]:
        """
        Return an iterator over the beads in the chain.

        Returns:
            Iterator[Bead]: Iterator over the chain's beads.

        """
        return iter(self.beads)

    def __getitem__(self, index: int) -> Bead:
        """
        Return the bead at the specified index in the chain.

        Args:
            index (int): Position of the bead in the chain.

        Returns:
            Bead: Bead instance at the given index.

        """
        return self.beads[index]

    def __len__(self) -> int:
        """
        Return the number of beads in the chain.

        Returns:
            int: Total number of beads in the chain.

        """
        return len(self.beads)

    def __str__(self) -> str:
        """
        Returns a string representation of the chain as a sequence of bead symbols.

        Returns:
            str: Concatenated sequence of bead symbols.

        """
        return "".join(bead.symbol for bead in self.beads)
