from collections.abc import Iterator

from logger import get_logger
from protein.bead import Bead
from protein.bead.main_bead import MainBead
from protein.chain import Chain

logger = get_logger()


class MainChain(Chain):
    """Represents the main chain of a protein, a linear sequence of amino acids forming its backbone."""

    def __init__(self, protein_sequence: str) -> None:
        """
        Initialize the main chain with beads corresponding to the protein sequence.

        Args:
            protein_sequence (str): Amino acid sequence of the protein backbone.

        """
        super().__init__(protein_sequence=protein_sequence)
        logger.debug(
            "Initializing MainChain with protein sequence: %s", protein_sequence
        )

        self.beads = [
            MainBead(
                symbol=bead,
                index=index,
                parent_chain_len=len(protein_sequence),
            )
            for index, bead in enumerate(protein_sequence)
        ]

    def __iter__(self) -> Iterator[Bead]:
        """
        Return an iterator over the beads in the main chain.

        Returns:
            Iterator[Bead]: Iterator over the chain's beads.

        """
        return iter(self.beads)

    def __getitem__(self, index: int) -> Bead:
        """
        Return the bead at the specified index in the main chain.

        Args:
            index (int): Position of the bead in the chain.

        Returns:
            Bead: Bead instance at the given index.

        """
        return self.beads[index]

    def __len__(self) -> int:
        """
        Return the number of beads in the main chain.

        Returns:
            int: Total number of beads in the chain.

        """
        return len(self.beads)

    def __str__(self) -> str:
        """
        Returns a string representation of the main chain as a sequence of bead symbols.

        Returns:
            str: Concatenated sequence of bead symbols.

        """
        return "".join(bead.symbol for bead in self.beads)
