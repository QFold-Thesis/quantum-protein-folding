from protein.chain import Chain
from protein.bead.side_bead import SideBead
from constants import EMPTY_SIDECHAIN_PLACEHOLDER


class SideChain(Chain):
    def __init__(self, protein_sequence: str) -> None:
        self.beads = [
            SideBead(symbol=bead, index=i)
            for i, bead in enumerate(protein_sequence)
            if bead != EMPTY_SIDECHAIN_PLACEHOLDER
        ]
