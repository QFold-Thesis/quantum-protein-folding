from dataclasses import dataclass

import numpy as np
from qiskit_algorithms import SamplingMinimumEigensolverResult

from constants import (
    CONFORMATION_ENCODING,
    DENSE_TURN_INDICATORS,
    QUBITS_PER_TURN,
    SPARSE_TURN_INDICATORS,
)
from enums import ConformationEncoding, TurnDirection
from exceptions import ConformationEncodingError
from logger import get_logger

logger = get_logger()


@dataclass
class VQEOutput:
    bitstring: str
    probability: float
    state: str
    energy_value: np.complex128

    def __repr__(self) -> str:
        return (
            f"VQEOutput(bitstring={self.bitstring},\n"
            f"  probability={self.probability},\n"
            f"  state={self.state},\n"
            f"  energy_value={self.energy_value})"
        )


def interpret_raw_vqe_output(raw_output: SamplingMinimumEigensolverResult) -> VQEOutput:
    best_measurement = raw_output.best_measurement

    if not best_measurement:
        msg = "No best measurement found in VQE output."
        raise ValueError(msg)

    bitstring: str | None = best_measurement.get("bitstring")
    probability: float | None = best_measurement.get("probability")
    state: str | None = best_measurement.get("state")
    energy_value: np.complex128 | None = best_measurement.get("value")

    if None in (bitstring, probability, state, energy_value):
        msg = "Incomplete best measurement data in VQE output."
        raise ValueError(msg)

    return VQEOutput(
        bitstring=bitstring,
        probability=probability,
        state=state,
        energy_value=energy_value,
    )


def _preprocess_bitstring(
    bitstring: str, turn_encoding: dict[TurnDirection, str]
) -> str:
    """Preprocesses the bitstring by appending initial turns and reversing it."""
    return "".join(
        reversed(
            bitstring
            + turn_encoding[TurnDirection.DIR_1]
            + turn_encoding[TurnDirection.DIR_2]
        )
    )


def generate_coords_from_bitstring(bitstring: str):
    if CONFORMATION_ENCODING == ConformationEncoding.DENSE:
        turn_encoding = DENSE_TURN_INDICATORS
    elif CONFORMATION_ENCODING == ConformationEncoding.SPARSE:
        turn_encoding = SPARSE_TURN_INDICATORS
    else:
        raise ConformationEncodingError

    bitstring = _preprocess_bitstring(bitstring, turn_encoding)

    tetra_dirs = np.array(
        [[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]], dtype=float
    )
    tetra_dirs /= np.linalg.norm(tetra_dirs[0])

    bitstring_to_direction = {
        bitstring: direction.value for direction, bitstring in turn_encoding.items()
    }

    turns_length = len(bitstring) // QUBITS_PER_TURN
    chunks = [
        bitstring[i * QUBITS_PER_TURN : (i + 1) * QUBITS_PER_TURN]
        for i in range(turns_length)
    ]

    pos = np.array([0.0, 0.0, 0.0])
    coords = [tuple(float(x) for x in pos)]

    for ch in chunks:
        if ch not in bitstring_to_direction:
            logger.warning(f"Unknown turn encoding: {ch}")
            continue
        direction_idx = bitstring_to_direction[ch]
        direction = tetra_dirs[direction_idx]
        pos = pos + direction
        coords.append(tuple(float(x) for x in pos))

    return coords
