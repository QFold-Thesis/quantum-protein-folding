from logger import get_logger
from protein.bead.main_bead import _MainBead
from protein.chain import Chain

logger = get_logger()


class _MainChain(Chain):
    """
    Represents the main chain of a protein, a linear sequence of amino acids forming its backbone.

    Attributes:
        beads (list[_MainBead]): List of main beads in the protein's backbone.

    """

    def __init__(self, protein_sequence: str) -> None:
        """
        Initialize the main chain with beads corresponding to the protein sequence.

        Args:
            protein_sequence (str): Amino acid sequence of the protein backbone.

        """
        logger.debug(
            f"Initializing MainChain based from protein sequence: {protein_sequence}..."
        )
        super().__init__(protein_sequence=protein_sequence)

    def _initialize_beads(self, protein_sequence: str) -> None:
        """
        Initialize main beads (_MainBead) based on the protein sequence.

        Args:
            protein_sequence (str): The amino acid sequence representing the protein chain.

        """
        self.beads = [
            _MainBead(
                symbol=bead,
                index=index,
                parent_chain_len=len(protein_sequence),
            )
            for index, bead in enumerate(protein_sequence)
        ]

        logger.info(
            f"MainChain for {protein_sequence} initialized with {len(self.beads)} MainBeads."
        )
