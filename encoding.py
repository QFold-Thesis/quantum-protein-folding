from enum import IntEnum, Enum
import numpy as np
import itertools


class SubLattice(IntEnum):
    A = 0
    B = 1


class Turn(Enum):
    DIR_0 = (1, 1, 1)
    DIR_1 = (1, -1, -1)
    DIR_2 = (-1, 1, -1)
    DIR_3 = (-1, -1, 1)


class Bead:
    def __init__(self, symbol: str, index: int):
        self.symbol: str = symbol
        self.index: int = index

        self.sublattice: SubLattice = (
            SubLattice.B.value if index % 2 == 1 else SubLattice.A.value
        )


class Chain:
    def __init__(self, beads: list[Bead]):
        self.beads: list[Bead] = beads
        self.length: int = len(beads)

    def __getitem__(self, index):
        return self.beads[index]

    def __len__(self):
        return self.length


class TetrahedralLattice:
    def __init__(self, size_x, size_y, size_z):
        self.size_x = size_x
        self.size_y = size_y
        self.size_z = size_z
        self.lattice = self._generate_lattice()

        self.directions = [
            np.array(Turn.DIR_0.value),
            np.array(Turn.DIR_1.value),
            np.array(Turn.DIR_2.value),
            np.array(Turn.DIR_3.value),
        ]

    def _generate_lattice(self):
        return [
            [[{} for _ in range(self.size_z)] for _ in range(self.size_y)]
            for _ in range(self.size_x)
        ]

    def is_valid(self, pos):
        x, y, z = pos
        return 0 <= x < self.size_x and 0 <= y < self.size_y and 0 <= z < self.size_z

    def get_neighbors(self, pos):
        x, y, z = pos
        neighbors = []
        for d in self.directions:
            nx, ny, nz = (np.array([x, y, z]) + d).tolist()
            if self.is_valid((nx, ny, nz)):
                neighbors.append((nx, ny, nz))
        return neighbors

    def encode_turn_sequence(self, turns):
        qubit_string = ""
        for turn in turns:
            binary = format(turn, "02b")
            qubit_string += binary
        return qubit_string

    def generate_positions(self, start_pos, turns):
        pos_list = [np.array(start_pos)]
        for t in turns:
            new_pos = pos_list[-1] + self.directions[t]
            if not self.is_valid(new_pos):
                raise ValueError(f"Position out of bounds: {new_pos}")
            pos_list.append(new_pos)
        return pos_list


def all_turn_combinations(seq_len):
    n = seq_len - 1
    return list(itertools.product(range(4), repeat=n))
