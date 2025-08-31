from constants import EMPTY_SIDECHAIN_PLACEHOLDER
from protein.bead.side_bead import SideBead
from protein.chain import Chain


class SideChain(Chain):
    def __init__(self, protein_sequence: str) -> None:
        self.beads = [
            SideBead(symbol=bead, index=i, chain_length=len(protein_sequence))
            for i, bead in enumerate(protein_sequence)
            if bead != EMPTY_SIDECHAIN_PLACEHOLDER
        ]

    @staticmethod
    def build_turn_qubit() -> None:
        pass
