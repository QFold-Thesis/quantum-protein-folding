from protein.bead import Bead


class Chain:
    def __init__(self, protein_sequence: str) -> None:
        self.beads: list[Bead] = [
            Bead(symbol=bead, index=i) for i, bead in enumerate(protein_sequence)
        ]
