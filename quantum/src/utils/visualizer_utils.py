import numpy as np

from constants import CONFORMATION_ENCODING, QUBITS_PER_TURN
from enums import ConformationEncoding
from exceptions import ConformationEncodingError
from logger import get_logger

logger = get_logger()


def generate_coords_from_bitstring(bitstring: str):
    tetra_dirs = np.array([
        [ 1,  1,  1],
        [ 1, -1, -1],
        [-1,  1, -1],
        [-1, -1,  1]
    ], dtype=float)

    tetra_dirs /= np.linalg.norm(tetra_dirs[0])

    if CONFORMATION_ENCODING == ConformationEncoding.DENSE:
        bit_to_dir = {
            "00": 0,
            "01": 1,
            "10": 2,
            "11": 3
        }
    elif CONFORMATION_ENCODING == ConformationEncoding.SPARSE:
        bit_to_dir = {
            "0001": 0,
            "0010": 1,
            "0100": 2,
            "1000": 3,
        }
    else:
        raise ConformationEncodingError

    chunks = [bitstring[i:i+QUBITS_PER_TURN] for i in range(0, len(bitstring), QUBITS_PER_TURN)]
    if len(chunks[-1]) != QUBITS_PER_TURN:
        chunks = chunks[:-1]

    pos = np.array([0.0, 0.0, 0.0])
    coords = [tuple(float(x) for x in pos)]

    for ch in chunks:
        direction = tetra_dirs[bit_to_dir[ch]]

        pos = pos + direction
        coords.append(tuple(float(x) for x in pos))

    return coords
