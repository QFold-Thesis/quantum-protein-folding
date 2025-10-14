from qiskit import QiskitError
import pytest
from qiskit.quantum_info import SparsePauliOp

from src.utils.qubit_utils import (
    build_identity_op,
    build_turn_qubit,
    build_pauli_z_operator,
    convert_to_qubits,
    fix_qubits,
    pad_to_n_qubits,
)

# ---- Fixtures --------------------------------------------------------------


@pytest.fixture
def sample_pauli():
    return SparsePauliOp.from_list([("ZI", 1.0)])


@pytest.fixture
def single_term_op():
    return SparsePauliOp.from_list([("ZZZZZZ", 1.0)])


@pytest.fixture
def multi_term_op():
    return SparsePauliOp.from_list([("ZZZZZZ", 1.0), ("IIIIII", 0.5)])


# ---- build_identity_op -----------------------------------------------------


def test_build_identity_op_returns_correct_identity():
    op = build_identity_op(3, coeff=2.0)
    assert isinstance(op, SparsePauliOp)
    assert op.num_qubits == 3
    label, coeff = op.to_list()[0]
    assert label == "III"
    assert coeff == 2.0


# ---- build_turn_qubit ------------------------------------------------------


def test_build_turn_qubit_creates_valid_operator():
    z_index = 1
    num_qubits = 4
    op = build_turn_qubit(z_index, num_qubits)

    labels = [label for label, _ in op.to_list()]
    print("Pauli label:", labels)

    assert labels[0] == "IIII"
    assert labels[1] == "IIZI"


def test_build_turn_qubit_creates_invalid_operator():
    z_index = 3
    num_qubits = 2
    with pytest.raises(QiskitError):
        build_turn_qubit(z_index, num_qubits)


# ---- build_pauli_z_operator ------------------------------------------------


def test_build_pauli_z_operator_nonempty():
    op = build_pauli_z_operator(4, {0, 2})
    label, _ = op.to_list()[0]

    print(op.to_list())
    assert label == "IZIZ"
    assert op.num_qubits == 4


def test_build_pauli_z_operator_empty():
    op = build_pauli_z_operator(3, set())
    label, coeff = op.to_list()[0]
    assert label == "III"
    assert coeff == 1.0


# ---- convert_to_qubits -----------------------------------------------------


def test_convert_to_qubits_valid(sample_pauli):
    converted = convert_to_qubits(sample_pauli)
    assert isinstance(converted, SparsePauliOp)
    assert converted.num_qubits == 2


# ---------------------------------------------------------------------------
# pad_to_n_qubits
# ---------------------------------------------------------------------------


def test_pad_to_n_qubits_adds_padding(sample_pauli):
    padded = pad_to_n_qubits(sample_pauli, 4)
    assert isinstance(padded, SparsePauliOp)
    assert padded.num_qubits == 4


def test_pad_to_n_qubits_same_size(sample_pauli):
    same = pad_to_n_qubits(sample_pauli, 2)
    assert same.num_qubits == 2


# ---------------------------------------------------------------------------
# fix_qubits
# ---------------------------------------------------------------------------


def test_fix_qubits_single_term(single_term_op):
    pass


def test_fix_qubits_multi_term(multi_term_op):
    pass


def test_fix_qubits_with_side_chain(single_term_op):
    pass


def test_fix_qubits_raises_on_invalid_operator():
    with pytest.raises(AttributeError):
        fix_qubits(None)


# ---------------------------------------------------------------------------
# unused_qubits
# ---------------------------------------------------------------------------

# todo
