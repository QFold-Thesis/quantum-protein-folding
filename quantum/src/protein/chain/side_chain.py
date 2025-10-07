from constants import EMPTY_SIDECHAIN_PLACEHOLDER
from logger import get_logger
from protein.bead.side_bead import SideBead
from protein.chain import Chain

logger = get_logger()


class SideChain(Chain):
    """Represents the side chain of a protein, consisting of amino acids attached to residues of the main chain."""

    def __init__(self, protein_sequence: str) -> None:
        """
        Initialize the side chain with beads corresponding to the protein sequence.

        Only residues that are not EMPTY_SIDECHAIN_PLACEHOLDER are included.

        Args:
            protein_sequence (str): Amino acid sequence of the protein side chain.

        """
        super().__init__(protein_sequence=protein_sequence)
        logger.debug(
            "Initializing SideChain with protein sequence: %s", protein_sequence
        )

        self.beads = [
            SideBead(
                symbol=bead,
                index=index,
                parent_chain_len=len(protein_sequence),
            )
            for index, bead in enumerate(protein_sequence)
            if bead != EMPTY_SIDECHAIN_PLACEHOLDER
        ]

    def __str__(self) -> str:
        """
        Return a string representation of the side chain as a sequence of bead symbols.

        Returns:
            str: Concatenated sequence of bead symbols in the side chain.

        """
        return "".join(bead.symbol for bead in self.beads)
