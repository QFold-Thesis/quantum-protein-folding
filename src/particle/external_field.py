"""External field acting on the peptide chain.

Two coupling modes are offered, and the difference between them is the whole
point of this module.

UNIFORM assigns every residue the same energy irrespective of where it sits.
Summed over the chain this produces a multiple of the identity operator, so it
shifts every conformation by the same constant and the ground state cannot
move. It is retained deliberately as a *null baseline*: any analysis that shows
a response to a uniform field is measuring an artefact.

GRADIENT couples to the actual lattice position of each residue,

    H_field = -strength * sum_i (direction . r_i)

where ``r_i`` is bead ``i``'s position operator built from the turn qubits.
Because ``r_i`` genuinely depends on the conformation, this term reorders the
energy spectrum and the optimal fold does change with *strength* - a chain in a
strong gradient stretches along the field, while a weak field leaves the
compact hydrophobic core intact. It costs no additional qubits: the field reads
the turn register that already exists.
"""

from __future__ import annotations

import math
from typing import TYPE_CHECKING

import numpy as np

from constants import (
    DEFAULT_FIELD_DIRECTION,
    DIST_VECTOR_AXES,
    EMPTY_OP_COEFF,
    FCC_BASIS,
)
from enums import FieldMode
from logger import get_logger
from utils.lattice_utils import build_chain_position
from utils.qubit_utils import build_identity_op

if TYPE_CHECKING:
    from numpy.typing import NDArray
    from qiskit.quantum_info import SparsePauliOp

logger = get_logger()


class ExternalField:
    """A position-dependent or uniform external field acting on the chain.

    Attributes:
        mode (FieldMode): Whether the field couples to bead positions (GRADIENT)
            or applies a flat per-bead shift (UNIFORM).
        strength (float): Magnitude of the field.
        direction (NDArray[np.float64]): Unit Cartesian direction of the
            gradient. Unused in UNIFORM mode.

    """

    def __init__(
        self,
        mode: FieldMode,
        strength: float,
        direction: NDArray[np.float64] | None = None,
    ) -> None:
        """Initialise an external field.

        Prefer the :meth:`uniform` and :meth:`gradient` factory methods.

        Args:
            mode (FieldMode): Coupling mode.
            strength (float): Field magnitude. Negative values invert the
                preferred direction.
            direction (NDArray[np.float64] | None, optional): Cartesian
                direction of the gradient, normalised on assignment. Defaults to
                DEFAULT_FIELD_DIRECTION.

        Raises:
            ValueError: If *strength* is not finite, or if *direction* has zero
                length or the wrong shape.

        """
        if not math.isfinite(strength):
            msg: str = f"strength must be a finite number, got {strength!r}"
            raise ValueError(msg)

        resolved_direction: NDArray[np.float64] = np.asarray(
            DEFAULT_FIELD_DIRECTION if direction is None else direction,
            dtype=np.float64,
        )

        if resolved_direction.shape != (FCC_BASIS.shape[1],):
            msg: str = (
                f"direction must have shape {(FCC_BASIS.shape[1],)}, "
                f"got {resolved_direction.shape}"
            )
            raise ValueError(msg)

        norm: float = float(np.linalg.norm(resolved_direction))
        if norm == 0.0:
            msg: str = "direction must be a non-zero vector"
            raise ValueError(msg)

        self.mode: FieldMode = mode
        self.strength: float = strength
        self.direction: NDArray[np.float64] = resolved_direction / norm

        logger.info(
            "ExternalField created [mode=%s, strength=%s, direction=%s]",
            mode.value,
            strength,
            self.direction.tolist(),
        )

    @classmethod
    def uniform(cls, strength: float) -> ExternalField:
        """Create a flat per-bead field that cannot change the ground state.

        Args:
            strength (float): Energy applied to every residue.

        Returns:
            ExternalField: A field in UNIFORM mode.

        """
        return cls(mode=FieldMode.UNIFORM, strength=strength)

    @classmethod
    def gradient(
        cls,
        strength: float,
        direction: NDArray[np.float64] | None = None,
    ) -> ExternalField:
        """Create a spatially varying field that couples to bead positions.

        Args:
            strength (float): Field magnitude. Positive values pull the chain
                along *direction*.
            direction (NDArray[np.float64] | None, optional): Cartesian
                direction of the gradient. Defaults to DEFAULT_FIELD_DIRECTION.

        Returns:
            ExternalField: A field in GRADIENT mode.

        """
        return cls(mode=FieldMode.GRADIENT, strength=strength, direction=direction)

    def axis_couplings(self) -> NDArray[np.float64]:
        """Project the field direction onto the four tetrahedral lattice axes.

        A bead position is stored as coefficients along the tetrahedral basis
        vectors, so the scalar product with the field direction reduces to a
        weighted sum of those coefficients.

        Returns:
            NDArray[np.float64]: Coupling weight for each tetrahedral axis.

        """
        basis: NDArray[np.float64] = FCC_BASIS / np.linalg.norm(FCC_BASIS[0])
        return basis @ self.direction

    def build_hamiltonian(self, chain_length: int, num_qubits: int) -> SparsePauliOp:
        """Build the field contribution to the Hamiltonian.

        Args:
            chain_length (int): Number of residues in the main chain.
            num_qubits (int): Width of the register the operator must span.

        Returns:
            SparsePauliOp: The field term. In UNIFORM mode this is a multiple of
            the identity; in GRADIENT mode it depends on the turn qubits.

        """
        if self.mode == FieldMode.UNIFORM:
            return self._build_uniform_term(chain_length, num_qubits)

        return self._build_gradient_term(chain_length, num_qubits)

    def _build_uniform_term(self, chain_length: int, num_qubits: int) -> SparsePauliOp:
        """Build the conformation-independent identity shift.

        Args:
            chain_length (int): Number of residues in the main chain.
            num_qubits (int): Width of the register the operator must span.

        Returns:
            SparsePauliOp: ``chain_length * strength`` times the identity.

        """
        total_energy: float = self.strength * chain_length

        logger.info(
            "Built uniform H_field: %s * I on %d qubits (ground state unaffected by construction)",
            total_energy,
            num_qubits,
        )
        return build_identity_op(num_qubits, total_energy)

    def _build_gradient_term(self, chain_length: int, num_qubits: int) -> SparsePauliOp:
        """Build the position-coupled field term.

        Args:
            chain_length (int): Number of residues in the main chain.
            num_qubits (int): Width of the register the operator must span.

        Returns:
            SparsePauliOp: Operator equal to ``-strength * sum_i direction . r_i``.

        """
        couplings: NDArray[np.float64] = self.axis_couplings()
        field_term: SparsePauliOp = build_identity_op(num_qubits, EMPTY_OP_COEFF)

        for bead_index in range(chain_length):
            position: list[SparsePauliOp] = build_chain_position(
                bead_index=bead_index, num_qubits=num_qubits
            )

            for axis in range(DIST_VECTOR_AXES):
                field_term = field_term - (
                    self.strength * float(couplings[axis]) * position[axis]
                )

        logger.info(
            "Built gradient H_field over %d beads [strength=%s, axis couplings=%s]",
            chain_length,
            self.strength,
            np.round(couplings, 4).tolist(),
        )
        return field_term.simplify()

    def __repr__(self) -> str:
        """Return a developer-readable representation of the field."""
        return (
            f"ExternalField(mode={self.mode.value!r}, strength={self.strength}, "
            f"direction={self.direction.tolist()})"
        )
