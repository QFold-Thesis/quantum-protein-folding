"""Defines the `SideChain` class for representing the side chain of a protein."""

from constants import EMPTY_SIDECHAIN_PLACEHOLDER
from logger import get_logger
from protein.bead.placeholder_side_bead import _PlaceholderSideBead
from protein.bead.side_bead import _SideBead
from protein.chain import Chain

logger = get_logger()


class _SideChain(Chain):
    """Represents the side chain of a protein, consisting of amino acids attached to residues of the main chain.

    Attributes:
        beads (list[_SideBead | _PlaceholderSideBead]): List of side beads in the protein's side chain.

    """

    def __init__(self, protein_sequence: str) -> None:
        """Initialize the side chain with beads corresponding to the protein sequence.

        If a bead symbol matches the EMPTY_SIDECHAIN_PLACEHOLDER, a PlaceholderSideBead is created instead.

        Args:
            protein_sequence (str): Amino acid sequence of the protein side chain.

        """
        logger.debug(
            "Initializing SideChain from protein sequence: %s...", protein_sequence
        )
        super().__init__(protein_sequence=protein_sequence)

    def _initialize_beads(self, protein_sequence: str) -> None:
        """
        Initialize side beads (_SideBead) and placeholder side beads (_PlaceholderSideBead) based on the protein sequence.

        Args:
            protein_sequence (str): The amino acid sequence representing the protein chain.

        """
        self.beads = [
            _SideBead(
                _symbol=bead,
                _index=index,
                _parent_chain_len=len(protein_sequence),
            )
            if bead != EMPTY_SIDECHAIN_PLACEHOLDER
            else _PlaceholderSideBead(
                symbol=bead,
                index=index,
                parent_chain_len=len(protein_sequence),
            )
            for index, bead in enumerate(protein_sequence)
        ]

        side_bead_count: int = len(
            [bead for bead in self.beads if not isinstance(bead, _PlaceholderSideBead)]
        )
        placeholder_bead_count: int = len(self.beads) - side_bead_count

        logger.info(
            "SideChain for %s initialized with %d SideBeads and %d PlaceholderSideBeads.",
            protein_sequence,
            side_bead_count,
            placeholder_bead_count,
        )
