import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

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


def _preprocess_bitstring(
    bitstring: str, turn_encoding: dict[TurnDirection, str]
) -> str:
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


def visualize_3d(coords, color="blue", figsize=(8, 8)):
    coords = np.array(coords)

    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111, projection="3d")

    base_dirs = np.array(
        [[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]], dtype=float
    )
    base_dirs /= np.linalg.norm(base_dirs[0])

    lattice_points = []
    n = 4
    for i in range(-n, n + 1):
        for j in range(-n, n + 1):
            for k in range(-n, n + 1):
                point = i * base_dirs[0] + j * base_dirs[1] + k * base_dirs[2]
                lattice_points.append(point)
    lattice_points = np.unique(np.round(lattice_points, 3), axis=0)

    center = np.mean(coords, axis=0)
    lattice_points += center

    ax.scatter(
        lattice_points[:, 0],
        lattice_points[:, 1],
        lattice_points[:, 2],
        s=25,
        color="lightgray",
        alpha=0.35,
        depthshade=False,
        label="Tetrahedral lattice",
    )

    ax.plot(
        coords[:, 0],
        coords[:, 1],
        coords[:, 2],
        "-",
        color=color,
        lw=3,
        alpha=0.9,
        label="Protein chain",
    )
    ax.scatter(
        coords[:, 0],
        coords[:, 1],
        coords[:, 2],
        s=70,
        color=color,
        edgecolor="black",
        alpha=0.9,
    )

    for i, (x, y, z) in enumerate(coords):
        ax.text(x, y, z, str(i), color="black", fontsize=8, ha="center", va="center")

    max_range = (coords.max(axis=0) - coords.min(axis=0)).max() / 2.0
    mid = coords.mean(axis=0)
    ax.set_xlim(mid[0] - max_range, mid[0] + max_range)
    ax.set_ylim(mid[1] - max_range, mid[1] + max_range)
    ax.set_zlim(mid[2] - max_range, mid[2] + max_range)

    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")
    ax.set_title("3D Protein Folding on Tetrahedral (Diamond) Lattice")

    ax.legend(loc="upper left")
    plt.tight_layout()
    plt.show()
