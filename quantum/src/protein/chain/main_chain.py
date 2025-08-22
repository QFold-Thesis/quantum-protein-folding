from protein.chain import Chain
from protein.bead.main_bead import MainBead


class MainChain(Chain):
    def __init__(self, protein_sequence: str) -> None:
        self.beads = [
            MainBead(symbol=bead, index=i) for i, bead in enumerate(protein_sequence)
        ]
