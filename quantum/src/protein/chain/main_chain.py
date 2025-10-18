from logger import get_logger
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
