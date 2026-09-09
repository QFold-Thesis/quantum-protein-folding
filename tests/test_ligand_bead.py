"""Tests for :class:`~particle.ligand_bead.LigandBead`.

Covers:
* Correct initialisation for both encoding strategies.
* Correct ``num_position_qubits`` for a range of lattice sizes.
* Type and shape of ``position_qubits`` tuple.
* Validity of Pauli operators (correct qubit count, SparsePauliOp type).
* ``position_projector`` for UNARY and BINARY encodings.
* Validation errors (bad symbol, negative index, too-few nodes, wrong encoding type).
* ``__repr__`` sanity check.
"""

from __future__ import annotations

import math

import pytest
from qiskit.quantum_info import SparsePauliOp

from particle.ligand_bead import LigandBead, PositionEncoding

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_binary(num_nodes: int, symbol: str = "L", index: int = 0) -> LigandBead:
    return LigandBead(
        symbol=symbol,
        index=index,
        num_lattice_nodes=num_nodes,
        encoding=PositionEncoding.BINARY,
    )


def _make_unary(num_nodes: int, symbol: str = "L", index: int = 0) -> LigandBead:
    return LigandBead(
        symbol=symbol,
        index=index,
        num_lattice_nodes=num_nodes,
        encoding=PositionEncoding.UNARY,
    )


# ---------------------------------------------------------------------------
# Initialisation - basic attributes
# ---------------------------------------------------------------------------


class TestLigandBeadInit:
    """Verify that LigandBead stores all constructor arguments correctly."""

    def test_symbol_stored(self):
        lig = _make_binary(8, symbol="ATP")
        assert lig.symbol == "ATP"

    def test_index_stored(self):
        lig = _make_binary(8, index=3)
        assert lig.index == 3

    def test_encoding_stored_binary(self):
        lig = _make_binary(8)
        assert lig.encoding is PositionEncoding.BINARY

    def test_encoding_stored_unary(self):
        lig = _make_unary(8)
        assert lig.encoding is PositionEncoding.UNARY

    def test_num_lattice_nodes_stored(self):
        lig = _make_binary(16)
        assert lig.num_lattice_nodes == 16

    def test_default_encoding_is_binary(self):
        lig = LigandBead(symbol="L", index=0, num_lattice_nodes=4)
        assert lig.encoding is PositionEncoding.BINARY


# ---------------------------------------------------------------------------
# num_position_qubits - BINARY encoding
# ---------------------------------------------------------------------------


class TestNumPositionQubitsBinary:
    """Verify ⌈log₂ N⌉ qubit calculation for BINARY encoding."""

    @pytest.mark.parametrize(
        "num_nodes, expected_qubits",
        [
            (2, 1),  # ⌈log₂ 2⌉ = 1
            (3, 2),  # ⌈log₂ 3⌉ = 2
            (4, 2),  # ⌈log₂ 4⌉ = 2
            (5, 3),  # ⌈log₂ 5⌉ = 3
            (8, 3),  # ⌈log₂ 8⌉ = 3
            (9, 4),  # ⌈log₂ 9⌉ = 4
            (16, 4),  # ⌈log₂ 16⌉ = 4
            (32, 5),
            (64, 6),
            (512, 9),
        ],
    )
    def test_num_qubits(self, num_nodes: int, expected_qubits: int):
        lig = _make_binary(num_nodes)
        assert lig.num_position_qubits == expected_qubits

    def test_matches_ceil_log2(self):
        """Cross-check against math.ceil(math.log2(N)) for various N."""
        for n in range(2, 50):
            lig = _make_binary(n)
            expected = max(1, math.ceil(math.log2(n)))
            assert lig.num_position_qubits == expected, (
                f"Failed for num_nodes={n}: got {lig.num_position_qubits}, expected {expected}"
            )


# ---------------------------------------------------------------------------
# num_position_qubits - UNARY encoding
# ---------------------------------------------------------------------------


class TestNumPositionQubitsUnary:
    """Verify that UNARY encoding uses exactly N qubits for N nodes."""

    @pytest.mark.parametrize("num_nodes", [2, 3, 4, 8, 16, 32])
    def test_num_qubits_equals_num_nodes(self, num_nodes: int):
        lig = _make_unary(num_nodes)
        assert lig.num_position_qubits == num_nodes


# ---------------------------------------------------------------------------
# position_qubits tuple
# ---------------------------------------------------------------------------


class TestPositionQubitsTuple:
    """Verify the structure and content of the position_qubits tuple."""

    def test_tuple_length_binary(self):
        lig = _make_binary(8)  # 3 qubits
        assert len(lig.position_qubits) == 3

    def test_tuple_length_unary(self):
        lig = _make_unary(6)
        assert len(lig.position_qubits) == 6

    def test_all_elements_are_sparse_pauli_op_binary(self):
        lig = _make_binary(16)
        for op in lig.position_qubits:
            assert isinstance(op, SparsePauliOp), (
                f"Expected SparsePauliOp, got {type(op)}"
            )

    def test_all_elements_are_sparse_pauli_op_unary(self):
        lig = _make_unary(4)
        for op in lig.position_qubits:
            assert isinstance(op, SparsePauliOp)

    def test_qubit_register_size_binary(self):
        """Each Pauli op must act on num_position_qubits qubits."""
        lig = _make_binary(8)  # 3 qubits
        for op in lig.position_qubits:
            assert op.num_qubits == lig.num_position_qubits

    def test_qubit_register_size_unary(self):
        lig = _make_unary(5)
        for op in lig.position_qubits:
            assert op.num_qubits == lig.num_position_qubits


# ---------------------------------------------------------------------------
# position_projector
# ---------------------------------------------------------------------------


class TestPositionProjector:
    """Verify position_projector for both encoding strategies."""

    # -- UNARY --

    def test_unary_projector_is_position_qubit(self):
        lig = _make_unary(4)
        for k in range(4):
            proj = lig.position_projector(node_index=k)
            # In UNARY mode the projector IS the k-th position qubit
            assert proj == lig.position_qubits[k]

    def test_unary_projector_type(self):
        lig = _make_unary(4)
        assert isinstance(lig.position_projector(0), SparsePauliOp)

    def test_unary_projector_num_qubits(self):
        lig = _make_unary(6)
        for k in range(6):
            assert lig.position_projector(k).num_qubits == 6

    # -- BINARY --

    def test_binary_projector_type(self):
        lig = _make_binary(4)
        assert isinstance(lig.position_projector(0), SparsePauliOp)

    def test_binary_projector_num_qubits(self):
        lig = _make_binary(8)  # 3 qubits
        for k in range(8):
            proj = lig.position_projector(k)
            assert proj.num_qubits == lig.num_position_qubits

    def test_binary_projectors_distinct(self):
        """Projectors onto different nodes must be different operators."""
        lig = _make_binary(4)
        proj_0 = lig.position_projector(0)
        proj_1 = lig.position_projector(1)
        proj_2 = lig.position_projector(2)
        proj_3 = lig.position_projector(3)
        # The operators should be mutually distinguishable
        assert proj_0 != proj_1
        assert proj_0 != proj_2
        assert proj_1 != proj_3

    def test_binary_projectors_sum_to_identity(self):
        """Sum of all projectors over a power-of-2 node space = identity."""
        num_nodes = 4
        lig = _make_binary(num_nodes)
        total = sum(lig.position_projector(k) for k in range(num_nodes)).simplify()
        identity = SparsePauliOp.from_list([("II", 1.0)])
        # Coefficients and Pauli strings must match
        assert total.equiv(identity)

    # -- Boundary / error --

    def test_projector_raises_for_negative_index(self):
        lig = _make_binary(4)
        with pytest.raises(ValueError, match="node_index"):
            lig.position_projector(-1)

    def test_projector_raises_for_out_of_range_index_unary(self):
        lig = _make_unary(4)
        with pytest.raises(ValueError, match="node_index"):
            lig.position_projector(4)

    def test_projector_raises_for_out_of_range_index_binary(self):
        lig = _make_binary(4)
        with pytest.raises(ValueError, match="node_index"):
            lig.position_projector(5)


# ---------------------------------------------------------------------------
# Validation / error handling
# ---------------------------------------------------------------------------


class TestLigandBeadValidation:
    """Verify that invalid arguments raise appropriate exceptions."""

    def test_empty_symbol_raises_value_error(self):
        with pytest.raises(ValueError, match="symbol"):
            LigandBead(symbol="", index=0, num_lattice_nodes=4)

    def test_negative_index_raises_value_error(self):
        with pytest.raises(ValueError, match="index"):
            LigandBead(symbol="L", index=-1, num_lattice_nodes=4)

    def test_num_nodes_less_than_2_raises_value_error(self):
        with pytest.raises(ValueError, match="num_lattice_nodes"):
            LigandBead(symbol="L", index=0, num_lattice_nodes=1)

    def test_num_nodes_zero_raises_value_error(self):
        with pytest.raises(ValueError, match="num_lattice_nodes"):
            LigandBead(symbol="L", index=0, num_lattice_nodes=0)

    def test_wrong_encoding_type_raises_type_error(self):
        with pytest.raises(TypeError, match="encoding"):
            LigandBead(symbol="L", index=0, num_lattice_nodes=4, encoding="binary")  # type: ignore[arg-type]


# ---------------------------------------------------------------------------
# Multiple ligands - independent qubit registers
# ---------------------------------------------------------------------------


class TestMultipleLigands:
    """Verify that separate LigandBead instances have independent state."""

    def test_two_ligands_independent(self):
        lig_a = LigandBead(symbol="A", index=0, num_lattice_nodes=4)
        lig_b = LigandBead(symbol="B", index=1, num_lattice_nodes=8)
        assert lig_a.num_position_qubits != lig_b.num_position_qubits
        assert lig_a.position_qubits is not lig_b.position_qubits

    def test_different_indices(self):
        lig0 = LigandBead(symbol="L", index=0, num_lattice_nodes=4)
        lig1 = LigandBead(symbol="L", index=1, num_lattice_nodes=4)
        assert lig0.index == 0
        assert lig1.index == 1


# ---------------------------------------------------------------------------
# __repr__
# ---------------------------------------------------------------------------


class TestRepr:
    def test_repr_contains_symbol(self):
        lig = LigandBead(symbol="XYZ", index=2, num_lattice_nodes=8)
        assert "XYZ" in repr(lig)

    def test_repr_contains_encoding_name(self):
        lig_b = _make_binary(8)
        assert "BINARY" in repr(lig_b)

        lig_u = _make_unary(8)
        assert "UNARY" in repr(lig_u)

    def test_repr_contains_index(self):
        lig = LigandBead(symbol="L", index=7, num_lattice_nodes=8)
        assert "7" in repr(lig)
