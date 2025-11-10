import pytest
from qiskit import QiskitError
from qiskit.quantum_info import SparsePauliOp

from src.constants import (
    MAIN_CHAIN_FIXED_POSITIONS,
    NORM_FACTOR,
    SIGN_FLIP_SECOND_QUBIT_INDEX,
)
from src.utils.qubit_utils import (
    build_identity_op,
    build_pauli_z_operator,
    build_turn_qubit,
    convert_to_qubits,
    find_unused_qubits,
    fix_qubits,
    pad_to_n_qubits,
    remove_unused_qubits,
)

# ---- Fixtures --------------------------------------------------------------


@pytest.fixture
def sample_pauli() -> SparsePauliOp:
    return SparsePauliOp.from_list([("ZI", 1.0)])


@pytest.fixture
def single_term_op() -> SparsePauliOp:
    return SparsePauliOp.from_list([("ZZZZZZ", 1.0)])


@pytest.fixture
def multi_term_op() -> SparsePauliOp:
    return SparsePauliOp.from_list([("ZZZZZZ", 1.0), ("IIIIII", 0.5)])


# ---- build_identity_op -----------------------------------------------------


def test_build_identity_op_returns_correct_identity():
    op = build_identity_op(3, coeff=2.0)
    assert isinstance(op, SparsePauliOp)
    assert op.num_qubits == 3  # noqa: PLR2004
    expected = SparsePauliOp.from_list([("III", 2.0)])
    assert op == expected


# ---- build_turn_qubit ------------------------------------------------------


def test_build_turn_qubit_creates_valid_operator():
    z_index = 1
    num_qubits = 4
    op = build_turn_qubit(z_index, num_qubits)

    expected = NORM_FACTOR * SparsePauliOp.from_list([("IIII", 1.0), ("IIZI", -1.0)])

    assert isinstance(op, SparsePauliOp)
    assert op.num_qubits == 4  # noqa: PLR2004
    assert op == expected


def test_build_turn_qubit_creates_invalid_operator():
    z_index = 3
    num_qubits = 2
    with pytest.raises(QiskitError):
        build_turn_qubit(z_index, num_qubits)


# ---- build_pauli_z_operator ------------------------------------------------


def test_build_pauli_z_operator_nonempty():
    op = build_pauli_z_operator(4, {0, 2})
    expected = SparsePauliOp.from_list([("IZIZ", 1.0)])
    assert op == expected


def test_build_pauli_z_operator_empty():
    op = build_pauli_z_operator(3, set())
    expected = SparsePauliOp.from_list([("III", 1.0)])
    assert op == expected


# ---- convert_to_qubits -----------------------------------------------------


def test_convert_to_qubits_valid(sample_pauli):
    converted = convert_to_qubits(sample_pauli)
    assert isinstance(converted, SparsePauliOp)
    assert converted.num_qubits == 2  # noqa: PLR2004


# ---------------------------------------------------------------------------
# pad_to_n_qubits
# ---------------------------------------------------------------------------


def test_pad_to_n_qubits_adds_padding(sample_pauli):
    padded = pad_to_n_qubits(sample_pauli, 4)
    expected = SparsePauliOp.from_list([("IIZI", 1.0)])
    assert isinstance(padded, SparsePauliOp)
    assert padded.num_qubits == 4  # noqa: PLR2004
    assert padded == expected


def test_pad_to_n_qubits_same_size(sample_pauli):
    same = pad_to_n_qubits(sample_pauli, 2)
    expected = SparsePauliOp.from_list([("ZI", 1.0)])
    assert same == expected


# ---------------------------------------------------------------------------
# fix_qubits
# ---------------------------------------------------------------------------


def test_fix_qubits_single_term(single_term_op):
    fixed = fix_qubits(single_term_op)
    assert isinstance(fixed, SparsePauliOp)
    assert fixed.num_qubits == 6  # noqa: PLR2004

    label, _ = fixed.to_list()[0]
    label_chars = list(label[::-1])

    assert all(label_chars[i] == "I" for i in MAIN_CHAIN_FIXED_POSITIONS)


def test_fix_qubits_multi_term(multi_term_op):
    fixed = fix_qubits(multi_term_op)
    assert isinstance(fixed, SparsePauliOp)
    assert fixed.num_qubits == 6  # noqa: PLR2004

    expected_terms = []
    for label, coeff in multi_term_op.to_list():
        if (
            len(label) > SIGN_FLIP_SECOND_QUBIT_INDEX
            and label[SIGN_FLIP_SECOND_QUBIT_INDEX] == "Z"
        ):
            coeff = -coeff  # noqa: PLW2901

        label_chars = list(label[::-1])
        for idx in MAIN_CHAIN_FIXED_POSITIONS:
            if idx < len(label_chars):
                label_chars[idx] = "I"
        new_label = "".join(label_chars[::-1])

        expected_terms.append((new_label, coeff))

    expected = SparsePauliOp.from_list(expected_terms)
    assert fixed == expected


def test_fix_qubits_raises_on_invalid_operator():
    with pytest.raises(AttributeError):
        fix_qubits(None)


# ---------------------------------------------------------------------------
# unused_qubits
# ---------------------------------------------------------------------------
def test_find_unused_qubits_sample(sample_pauli):
    unused = find_unused_qubits(sample_pauli)
    assert unused == [0]


def test_find_unused_qubits_single_term(single_term_op):
    unused = find_unused_qubits(single_term_op)
    assert unused == []


def test_remove_unused_qubits_removes_correctly():
    op = SparsePauliOp.from_list([("IZII", 1.0), ("IIII", 0.5)])
    reduced = remove_unused_qubits(op)
    expected = SparsePauliOp.from_list([("Z", 1.0), ("I", 0.5)])
    assert reduced == expected
