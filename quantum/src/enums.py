from enum import IntEnum


class ConformationEncoding(IntEnum):
    """Enum representing a map of encoding types and qubit counts."""

    SPARSE = 4
    DENSE = 2


class SubLattice(IntEnum):
    """Enum representing the sublattices in the protein chain."""

    A = 0
    B = 1


class Penalties(IntEnum):
    """Enum representing penalty types for protein folding constraints."""

    OVERLAP_PENALTY = 10
    CHIRALITY_PENALTY = 10
    BACK_PENALTY = 10


class InteractionType(IntEnum):
    """Enum representing interaction types."""

    MJ = 0
    HP = 1


class TurnDirection(IntEnum):
    """Enum representing turn directions on a tetrahedral lattice."""

    DIR_1 = 0
    DIR_2 = 1
    DIR_3 = 2
    DIR_4 = 3
