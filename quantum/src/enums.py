from enum import IntEnum


class ConformationEncoding(IntEnum):
    """Enum representing a map of encoding types and qubit counts."""

    SPARSE = 4
    DENSE = 2


class SubLattice(IntEnum):
    """Enum representing the sublattices in the protein chain."""

    A = 0
    B = 1
