from __future__ import annotations

from typing import TYPE_CHECKING, Any, cast

import matplotlib.pyplot as plt
import numpy as np

from logger import get_logger

if TYPE_CHECKING:
    from mpl_toolkits.mplot3d import Axes3D
    from numpy.typing import NDArray
    from result_interpretation_utils import BeadPosition

logger = get_logger()


def visualize_3d(
    coords: list[BeadPosition],
    *,
    color: str = "blue",
    figsize: tuple[int, int] = (8, 8),
) -> None:
    if len(coords) == 0:
        logger.warning("visualize_3d received empty coords; nothing to plot.")
        return

    coords_arr: NDArray[np.float64] = np.array(
        [bp.position for bp in coords], dtype=float
    )

    fig = plt.figure(figsize=figsize)
    ax: Axes3D = fig.add_subplot(111, projection="3d")  # type: ignore[assignment]
    ax3d = cast(Any, ax)

    base_dirs: NDArray[np.float64] = np.array(
        [[1.0, 1.0, 1.0], [1.0, -1.0, -1.0], [-1.0, 1.0, -1.0], [-1.0, -1.0, 1.0]],
        dtype=float,
    )
    base_dirs = base_dirs / float(np.linalg.norm(base_dirs[0]))

    lattice_points: list[NDArray[np.float64]] = []
    n: int = 4
    for i in range(-n, n + 1):
        for j in range(-n, n + 1):
            for k in range(-n, n + 1):
                point: NDArray[np.float64] = (
                    i * base_dirs[0] + j * base_dirs[1] + k * base_dirs[2]
                )
                lattice_points.append(point)
    lattice_points_arr: NDArray[np.float64] = np.unique(
        np.round(np.vstack(lattice_points), 3), axis=0
    )

    center: NDArray[np.float64] = np.mean(coords_arr, axis=0)
    lattice_points_arr = lattice_points_arr + center

    ax3d.scatter(
        lattice_points_arr[:, 0],
        lattice_points_arr[:, 1],
        lattice_points_arr[:, 2],
        s=25,
        color="lightgray",
        alpha=0.35,
        depthshade=False,
        label="Tetrahedral lattice",
    )

    ax3d.plot(
        coords_arr[:, 0],
        coords_arr[:, 1],
        coords_arr[:, 2],
        linestyle="-",
        color=color,
        lw=3,
        alpha=0.9,
        label="Protein chain",
    )
    ax3d.scatter(
        coords_arr[:, 0],
        coords_arr[:, 1],
        coords_arr[:, 2],
        s=70,
        color=color,
        edgecolor="black",
        alpha=0.9,
    )

    for i, (x, y, z) in enumerate(coords_arr):
        ax3d.text(
            float(x),
            float(y),
            float(z),
            str(i),
            color="black",
            fontsize=8,
            ha="center",
            va="center",
        )

    max_range: float = float(
        (coords_arr.max(axis=0) - coords_arr.min(axis=0)).max() / 2.0
    )
    mid: NDArray[np.float64] = coords_arr.mean(axis=0)
    ax3d.set_xlim(float(mid[0] - max_range), float(mid[0] + max_range))
    ax3d.set_ylim(float(mid[1] - max_range), float(mid[1] + max_range))
    ax3d.set_zlim(float(mid[2] - max_range), float(mid[2] + max_range))

    ax3d.set_xlabel("X")
    ax3d.set_ylabel("Y")
    ax3d.set_zlabel("Z")
    ax3d.set_title("3D Protein Folding on Tetrahedral (Diamond) Lattice")

    ax3d.legend(loc="upper left")
    plt.tight_layout()
    plt.show()
