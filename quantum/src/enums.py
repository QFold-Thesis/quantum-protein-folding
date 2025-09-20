from enum import IntEnum, FloatEnum


class ConformationEncoding(IntEnum):
    """Enum representing a map of encoding types and qubit counts."""

    SPARSE = 4
    DENSE = 2


class SubLattice(IntEnum):
    """Enum representing the sublattices in the protein chain."""

    A = 0
    B = 1

class Penalties(FloatEnum):
    """Enum representing penalty types for protein folding constraints."""

    OVERLAP_PENALTY = 10.0
    CHIRALITY_PENALTY = 10.0
    BACK_PENALTY = 10.0