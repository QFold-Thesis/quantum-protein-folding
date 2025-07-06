import itertools
from encoding import SubLattice


class Bead:
    def __init__(self, symbol: str, index: int):
        self.symbol: str = symbol
        self.index: int = index
        self.sublattice: SubLattice = (
            SubLattice.B if index % 2 == 1 else SubLattice.A
        )
    
    def __repr__(self):
        return f"Bead (index={self.index}, symbol={self.symbol}, sublattice={self.sublattice.name})"
    

def all_turn_combinations(seq_len: int):
    n_turns = seq_len - 3
    fixed_turns = [1, 0]
    for turns in itertools.product(range(4), repeat=n_turns):
        yield fixed_turns + list(turns)

def create_protein_sequence(sequence: list[str]) -> list[Bead]:
    beads: list[Bead] = []
    for i, aa in enumerate(sequence):
        bead = Bead(index=i, symbol=aa)
        beads.append(bead)
    return beads