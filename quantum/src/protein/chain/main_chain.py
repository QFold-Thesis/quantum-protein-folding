from collections.abc import Iterator

from logger import get_logger
from protein.bead import Bead
from protein.bead.main_bead import MainBead
from protein.chain import Chain

logger = get_logger()


class MainChain(Chain):
    """
    Represents the main chain of a protein, a linear sequence of amino acids forming its backbone.
    """

    def __init__(self, protein_sequence: str) -> None:
        """
        Initializes the main chain with beads corresponding to the protein sequence.
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
        """Returns an iterator over the beads in the main chain."""
        return iter(self.beads)

    def __getitem__(self, index: int) -> Bead:
        """Returns the bead at the specified index in the main chain."""
        return self.beads[index]

    def __len__(self) -> int:
        """Returns the number of beads in the main chain."""
        return len(self.beads)

    def __str__(self) -> str:
        """Returns a string representation of the main chain as a sequence of bead symbols."""
        return "".join(bead.symbol for bead in self.beads)
