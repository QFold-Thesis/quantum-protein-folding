"""Free ligand particle on the tetrahedral lattice.

The ligand is not bonded to the peptide, so it carries no backbone constraints.
Its position is nevertheless expressed in exactly the same language as the
chain: a walk of ``num_steps`` unit steps away from the chain origin, each step
choosing one of four tetrahedral directions. That shared language is what makes
a genuine spatial coupling possible - the squared distance between a residue and
the ligand is an operator built from both walks, so the Hamiltonian can reward
the ligand for actually sitting next to a residue rather than merely shifting
every energy by a constant.

Sublattice parity
-----------------
On the diamond lattice a site reached in ``m`` steps lies on the same sublattice
as backbone bead ``m``, and only sites of *opposite* sublattice can be nearest
neighbours. A ligand encoded with ``num_steps = m`` can therefore only contact
residues whose index has parity opposite to ``m``. This is a property of the
lattice, not a limitation of the encoding: running both parities explores the
two complementary sets of binding sites.

Register layout
---------------
The ligand occupies two contiguous slices above the protein's own qubits::

    [ protein turns | backbone contacts | ligand walk | ligand contacts ]

The walk slice holds ``num_steps * QUBITS_PER_TURN`` qubits encoding the
position. The contact slice holds one qubit per eligible residue, flagging which
residue the ligand claims to bind - the same boolean-contact device the backbone
term already uses.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from constants import (
    DEFAULT_LIGAND_STEPS,
    DEFAULT_LIGAND_SYMBOL,
    QUBITS_PER_TURN,
)
from logger import get_logger
from utils.lattice_utils import build_walk_position

if TYPE_CHECKING:
    from qiskit.quantum_info import SparsePauliOp

logger = get_logger()


class Ligand:
    """A free particle placed on the same tetrahedral lattice as the peptide.

    Attributes:
        symbol (str): Symbol identifying the ligand in results and visualisations.
        num_steps (int): Number of lattice steps encoding the ligand's position.

    """

    def __init__(
        self,
        symbol: str = DEFAULT_LIGAND_SYMBOL,
        num_steps: int = DEFAULT_LIGAND_STEPS,
    ) -> None:
        """Initialise the ligand.

        Args:
            symbol (str, optional): Symbol identifying the ligand. Defaults to
                DEFAULT_LIGAND_SYMBOL.
            num_steps (int, optional): Lattice steps encoding the ligand's
                position. More steps let the ligand roam further at the cost of
                two qubits each. Defaults to DEFAULT_LIGAND_STEPS.

        Raises:
            ValueError: If *num_steps* is not a positive integer.

        """
        if num_steps < 1:
            msg: str = f"num_steps must be a positive integer, got {num_steps!r}"
            raise ValueError(msg)

        self.symbol: str = symbol
        self.num_steps: int = num_steps

        logger.info(
            "Ligand %s initialised with %d lattice steps (%d walk qubits, sublattice %d)",
            self.symbol,
            self.num_steps,
            self.num_walk_qubits,
            self.sublattice_parity,
        )

    @property
    def num_walk_qubits(self) -> int:
        """int: Number of qubits encoding the ligand's lattice position."""
        return self.num_steps * QUBITS_PER_TURN

    @property
    def sublattice_parity(self) -> int:
        """int: Sublattice the ligand occupies, matching backbone bead ``num_steps``."""
        return self.num_steps % 2

    def eligible_bead_indices(self, chain_length: int) -> list[int]:
        """Return the residues the ligand can physically contact.

        Only residues on the opposite sublattice can be lattice nearest
        neighbours of the ligand.

        Args:
            chain_length (int): Number of residues in the main chain.

        Returns:
            list[int]: Indices of residues the ligand may bind, in chain order.

        """
        eligible: list[int] = [
            index
            for index in range(chain_length)
            if index % 2 != self.sublattice_parity
        ]

        logger.debug(
            "Ligand %s (%d steps) may contact residues %s of %d",
            self.symbol,
            self.num_steps,
            eligible,
            chain_length,
        )
        return eligible

    def position_operators(
        self, num_qubits: int, walk_qubit_offset: int
    ) -> list[SparsePauliOp]:
        """Build the ligand's axis-coefficient position operators.

        Args:
            num_qubits (int): Total width of the register the operators live on.
            walk_qubit_offset (int): Index of the first qubit of the ligand's
                walk slice.

        Returns:
            list[SparsePauliOp]: Four operators holding the ligand's coefficient
            along each tetrahedral direction.

        """
        return build_walk_position(
            qubit_bases=[
                walk_qubit_offset + QUBITS_PER_TURN * step
                for step in range(self.num_steps)
            ],
            num_qubits=num_qubits,
        )

    def __repr__(self) -> str:
        """Return a developer-readable representation of the ligand."""
        return f"Ligand(symbol={self.symbol!r}, num_steps={self.num_steps})"
