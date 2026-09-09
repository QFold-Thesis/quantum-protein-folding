"""Ligand bead representation for quantum protein folding.

A *ligand* (small molecule) is modelled as a single, free-floating bead that
occupies a node of the same lattice as the protein chain but is **not** covalently
bound to it. Unlike protein beads, a ligand does not have turn qubits; instead its
quantum state is described exclusively by *position qubits* that encode which
lattice node it currently occupies.

Two encoding strategies are supported and selectable at construction time via the
:class:`PositionEncoding` enum:

* **BINARY** (default) - uses ceil(log2 N) qubits to address N lattice nodes in
  standard binary (computational-basis) encoding.  This is the recommended choice:
  it is maximally qubit-efficient and maps naturally onto the Pauli-Z basis used by
  the rest of the Hamiltonian.  For a 4x4 2-D square lattice (16 nodes) only 4
  qubits are required, compared to 16 for unary.

* **UNARY** (one-hot) - uses exactly N qubits, one per node, where the
  computational state |0...010...0> with a single '1' at position k means "ligand is
  at node k".  Unary encoding simplifies the construction of interaction terms
  (the projector onto node k is simply (I - Z_k)/2) at the cost of many more
  qubits for large lattices.  It may be preferable for small lattices or when
  interaction Hamiltonians are the primary concern.

Design decision
---------------
``LigandBead`` is implemented as an **independent class** rather than a subclass of
:class:`~protein.bead.bead.Bead`.  The reasons are:

1. ``Bead`` is tightly coupled to the concept of a *turn qubit* and to the chain's
   ``parent_chain_len`` parameter; forcing a ligand to fit that interface would
   require either raising ``NotImplementedError`` for every abstract method or
   carrying dummy state.
2. The ligand's position qubits are structurally different from turn qubits
   (different physical meaning, different Hilbert-space dimension, no index
   dependence within a chain).
3. Keeping ``LigandBead`` independent avoids accidental coupling and keeps the door
   open for richer ligand models (multiple position qubits for flexible ligands,
   orientation qubits, …) without touching the protein-chain hierarchy.

Example usage::

    from particle.ligand_bead import LigandBead, PositionEncoding

    # Binary encoding for a 3-D cubic lattice with 8x8x8 = 512 nodes
    lig = LigandBead(symbol="L", index=0, num_lattice_nodes=512)
    print(lig.num_position_qubits)   # -> 9  (ceil(log2 512))

    # Unary encoding for a tiny 2x2 lattice (4 nodes)
    lig2 = LigandBead(
        symbol="L",
        index=1,
        num_lattice_nodes=4,
        encoding=PositionEncoding.UNARY,
    )
    print(lig2.num_position_qubits)  # → 4
    op = lig2.position_projector(node_index=2)  # projector onto node 2
"""

from __future__ import annotations

import math
from enum import Enum, auto
from typing import TYPE_CHECKING

from qiskit.quantum_info import SparsePauliOp

from logger import get_logger
from utils.qubit_utils import build_identity_op, build_turn_qubit

if TYPE_CHECKING:
    pass

logger = get_logger()

_MIN_LATTICE_NODES: int = 2


class PositionEncoding(Enum):
    """Encoding strategy for the ligand's position qubits.

    Attributes:
        BINARY: Standard binary encoding - ceil(log2 N) qubits for N nodes.
            Recommended for qubit efficiency; maps onto the same Pauli-Z basis
            used by protein turn qubits.
        UNARY: One-hot (unary) encoding - N qubits for N nodes.
            Simplifies projector construction; suitable for small lattices or
            when interaction Hamiltonian terms are the primary concern.

    """

    BINARY = auto()
    UNARY = auto()


class LigandBead:
    """A free-floating ligand particle on the lattice.

    The ligand occupies a single lattice node at any given time. Its quantum
    state is described by *position qubits* rather than turn qubits: there is
    no covalent bond to the protein chain.

    Attributes:
        symbol (str): Identifier symbol for the ligand (e.g. ``"L"``).
        index (int): Unique integer index distinguishing this ligand from
            others (if multiple ligands are present in future extensions).
        encoding (PositionEncoding): The chosen position-encoding strategy.
        num_lattice_nodes (int): Total number of addressable lattice nodes.
        num_position_qubits (int): Number of qubits used to represent the
            ligand's position under the selected encoding.
        position_qubits (tuple[SparsePauliOp, ...]): Pauli operators
            representing the individual position qubits.  For BINARY encoding
            these are single-qubit Z operators embedded in the full register;
            for UNARY encoding each operator is the projector onto the
            corresponding basis state.

    """

    # ------------------------------------------------------------------
    # Construction
    # ------------------------------------------------------------------

    def __init__(
        self,
        symbol: str,
        index: int,
        num_lattice_nodes: int,
        encoding: PositionEncoding = PositionEncoding.BINARY,
    ) -> None:
        """Initialise a :class:`LigandBead`.

        Args:
            symbol (str): Ligand identifier symbol.  Any non-empty string is
                accepted (e.g. ``"L"``, ``"ATP"``, ``"LIG"``).
            index (int): Non-negative integer index for this ligand.
            num_lattice_nodes (int): Number of distinct lattice nodes the
                ligand can occupy.  Must be ≥ 2.
            encoding (PositionEncoding, optional): Position-encoding strategy.
                Defaults to :attr:`PositionEncoding.BINARY`.

        Raises:
            ValueError: If *symbol* is empty, *index* is negative, or
                *num_lattice_nodes* is less than 2.
            TypeError: If *encoding* is not a :class:`PositionEncoding` member.

        """
        if not symbol:
            msg = "symbol must be a non-empty string."
            raise ValueError(msg)
        if index < 0:
            msg = f"index must be a non-negative integer, got {index!r}."
            raise ValueError(msg)
        if num_lattice_nodes < _MIN_LATTICE_NODES:
            msg = (
                f"num_lattice_nodes must be at least 2, got {num_lattice_nodes!r}. "
                "A ligand with fewer than 2 reachable nodes cannot move."
            )
            raise ValueError(msg)
        if not isinstance(encoding, PositionEncoding):
            msg = (
                f"encoding must be a PositionEncoding member, "
                f"got {type(encoding).__name__!r}."
            )
            raise TypeError(msg)

        self.symbol: str = symbol
        self.index: int = index
        self.encoding: PositionEncoding = encoding
        self.num_lattice_nodes: int = num_lattice_nodes
        self.num_position_qubits: int = self._calc_num_position_qubits()

        self.position_qubits: tuple[SparsePauliOp, ...] = self._build_position_qubits()

        logger.info(
            "LigandBead '%s' (index=%d) initialised: encoding=%s, "
            "lattice_nodes=%d, position_qubits=%d",
            self.symbol,
            self.index,
            self.encoding.name,
            self.num_lattice_nodes,
            self.num_position_qubits,
        )

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _calc_num_position_qubits(self) -> int:
        """Return the number of position qubits required for the chosen encoding.

        Returns:
            int: Number of position qubits.

        """
        if self.encoding is PositionEncoding.BINARY:
            # ceil(log2(num_lattice_nodes))  -  at least 1 qubit even for N=2
            return max(1, math.ceil(math.log2(self.num_lattice_nodes)))
        # UNARY: one qubit per node
        return self.num_lattice_nodes

    def _build_position_qubits(self) -> tuple[SparsePauliOp, ...]:
        """Build Pauli operators representing position qubits.

        For **BINARY** encoding each qubit k corresponds to the k-th bit of
        the binary representation of the node index.  The operator is the
        standard turn-qubit form ``half*(I - Z_k)`` already used by the protein
        chain, so the ligand naturally integrates into the same operator
        algebra.

        For **UNARY** encoding each qubit k corresponds to a projector onto
        the k-th one-hot basis state.  The number of qubits equals
        ``num_lattice_nodes``.

        Returns:
            tuple[SparsePauliOp, ...]: Tuple of position Pauli operators.

        """
        n = self.num_position_qubits

        if self.encoding is PositionEncoding.BINARY:
            # Re-use build_turn_qubit: half*(I - Z_k) at position k
            return tuple(build_turn_qubit(z_index=k, num_qubits=n) for k in range(n))

        # UNARY: projector half*(I - Z_k) is already a natural one-hot qubit
        # operator when each physical qubit is constrained to the one-hot
        # subspace.  We build one operator per lattice node.
        return tuple(build_turn_qubit(z_index=k, num_qubits=n) for k in range(n))

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def position_projector(self, node_index: int) -> SparsePauliOp:
        """Return the Pauli projector onto the lattice node *node_index*.

        For **UNARY** encoding this is simply ``position_qubits[node_index]``,
        i.e. the one-qubit projector ``half*(I - Z_{node_index})``.

        For **BINARY** encoding the projector is the tensor product of
        single-qubit projectors corresponding to the bits of *node_index* in
        the binary representation::

            P(k) = ∏_{b=0}^{n-1}  P_b(bit_b(k))

        where ``P_b(0) = half*(I - Z_b)`` and ``P_b(1) = half*(I + Z_b)``.

        Args:
            node_index (int): Index of the target lattice node.
                Must satisfy ``0 ≤ node_index < num_lattice_nodes``.

        Returns:
            SparsePauliOp: Projector onto the quantum state representing
                *node_index*.

        Raises:
            ValueError: If *node_index* is out of range.

        """
        if not (0 <= node_index < self.num_lattice_nodes):
            msg = (
                f"node_index must be in [0, {self.num_lattice_nodes}), "
                f"got {node_index!r}."
            )
            raise ValueError(msg)

        n = self.num_position_qubits

        if self.encoding is PositionEncoding.UNARY:
            return self.position_qubits[node_index]

        # BINARY: build product of single-qubit projectors
        full_identity: SparsePauliOp = build_identity_op(num_qubits=n)
        projector: SparsePauliOp = full_identity

        for bit_pos in range(n):
            bit_val = (node_index >> bit_pos) & 1
            z_op = SparsePauliOp.from_sparse_list([("Z", [bit_pos], 1.0)], num_qubits=n)
            if bit_val == 0:
                # projector onto |0>: half*(I - Z)
                single_proj = (0.5 * (full_identity - z_op)).simplify()
            else:
                # projector onto |1⟩: ½(I + Z)
                single_proj = (0.5 * (full_identity + z_op)).simplify()

            projector = (projector @ single_proj).simplify()

        return projector

    def __repr__(self) -> str:
        """Return a developer-readable string representation."""
        return (
            f"LigandBead(symbol={self.symbol!r}, index={self.index}, "
            f"encoding={self.encoding.name}, "
            f"num_lattice_nodes={self.num_lattice_nodes}, "
            f"num_position_qubits={self.num_position_qubits})"
        )
