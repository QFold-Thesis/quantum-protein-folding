"""Defines key enumerations for modeling protein folding constraints and encodings."""

from enum import Enum, IntEnum


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

    DIR_0 = 0
    DIR_1 = 1
    DIR_2 = 2
    DIR_3 = 3

    def __str__(self) -> str:
        """Return a string representation of the turn direction."""
        return f"Direction {self.value}"


class BackendType(Enum):
    """Enum representing quantum backend types for VQE execution."""

    LOCAL_STATEVECTOR = "local_statevector"
    IBM_QUANTUM = "ibm_quantum"


class FieldMode(Enum):
    """Enum representing how an external field couples to the peptide chain.

    UNIFORM assigns the same energy to every bead regardless of where it sits,
    which makes the term a multiple of the identity: it shifts every
    conformation equally and therefore cannot change the ground state. It is
    kept as an explicit baseline to contrast against GRADIENT.

    GRADIENT couples to the actual lattice position of each bead, so different
    conformations feel different energies and the optimum does move with the
    field strength.
    """

    UNIFORM = "uniform"
    GRADIENT = "gradient"


class LigandInteractionMode(Enum):
    """Enum representing how ligand-residue contact energies are defined.

    HP_LIKE reuses the hydrophobic/polar matrix by treating the ligand itself as
    either an H or a P residue. CUSTOM takes an explicit per-residue energy map.
    """

    HP_LIKE = "hp_like"
    CUSTOM = "custom"
