from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Any, cast

import matplotlib.pyplot as plt
import numpy as np

from constants import GIF_FILENAME
from logger import get_logger

if TYPE_CHECKING:
    from pathlib import Path

    from mpl_toolkits.mplot3d import Axes3D
    from numpy.typing import NDArray
    from result_interpretation_utils import BeadPosition

logger = get_logger()


def _compute_axes_limits(
    coords_arr: NDArray[np.float64],
) -> tuple[float, NDArray[np.float64]]:
    """Return (max_range, mid) for setting symmetric 3D limits."""
    max_range: float = float(
        (coords_arr.max(axis=0) - coords_arr.min(axis=0)).max() / 2.0
    )
    mid: NDArray[np.float64] = coords_arr.mean(axis=0)
    return max_range, mid


@dataclass(slots=True)
class AxesLimits:
    x: tuple[float, float]
    y: tuple[float, float]
    z: tuple[float, float]


@dataclass(slots=True)
class RotGifConfig:
    figsize: tuple[int, int]
    elev: float
    azim_start: float
    frames: int
    fps: int


@dataclass(slots=True)
class PlotScene:
    coords_arr: NDArray[np.float64]
    coords: list[BeadPosition]
    lattice_points_arr: NDArray[np.float64]
    bead_colors: NDArray[np.floating]
    mid: NDArray[np.float64]
    max_range: float


def _build_lattice_points(
    *,
    base_dirs: NDArray[np.float64],
    max_range: float,
    center: NDArray[np.float64],
    limits: AxesLimits,
) -> NDArray[np.float64]:
    """Generate tetrahedral lattice points centered at `center` and clip to axes limits."""
    lattice_points: list[NDArray[np.float64]] = []
    n: int = max(2, int(np.ceil(max_range)) + 2)
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
    lattice_points_arr = lattice_points_arr + center
    x_min, x_max = limits.x
    y_min, y_max = limits.y
    z_min, z_max = limits.z
    mask = (
        (lattice_points_arr[:, 0] >= x_min)
        & (lattice_points_arr[:, 0] <= x_max)
        & (lattice_points_arr[:, 1] >= y_min)
        & (lattice_points_arr[:, 1] <= y_max)
        & (lattice_points_arr[:, 2] >= z_min)
        & (lattice_points_arr[:, 2] <= z_max)
    )
    return lattice_points_arr[mask]


def _draw_chain(
    ax3d: Any, coords_arr: NDArray[np.float64], colors: NDArray[np.floating]
) -> None:
    for i in range(coords_arr.shape[0] - 1):
        ax3d.plot(
            coords_arr[i : i + 2, 0],
            coords_arr[i : i + 2, 1],
            coords_arr[i : i + 2, 2],
            linestyle="-",
            color=colors[i],
            lw=3,
            alpha=0.9,
            label="Protein chain" if i == 0 else None,
        )
    ax3d.scatter(
        coords_arr[:, 0],
        coords_arr[:, 1],
        coords_arr[:, 2],
        s=70,
        c=colors,
        edgecolor="black",
        alpha=0.95,
    )


def _annotate_beads(
    ax3d: Any,
    coords_arr: NDArray[np.float64],
    coords: list[BeadPosition],
    offset: float,
) -> None:
    for i, (x, y, z) in enumerate(coords_arr):
        symbol = coords[i].symbol
        ax3d.text(
            float(x + offset),
            float(y + offset),
            float(z + offset),
            symbol,
            color="black",
            fontsize=9,
            ha="center",
            va="center",
        )


def _save_rotating_gif(scene: PlotScene, cfg: RotGifConfig, *, dirpath: Path) -> None:
    try:
        from matplotlib import animation

        fig_anim = cast(Any, plt).figure(figsize=cfg.figsize)
        ax_anim: Axes3D = fig_anim.add_subplot(111, projection="3d")  # type: ignore[assignment]
        ax_anim3d = cast(Any, ax_anim)

        ax_anim3d.set_xlim(
            float(scene.mid[0] - scene.max_range), float(scene.mid[0] + scene.max_range)
        )
        ax_anim3d.set_ylim(
            float(scene.mid[1] - scene.max_range), float(scene.mid[1] + scene.max_range)
        )
        ax_anim3d.set_zlim(
            float(scene.mid[2] - scene.max_range), float(scene.mid[2] + scene.max_range)
        )
        if scene.lattice_points_arr.size > 0:
            ax_anim3d.scatter(
                scene.lattice_points_arr[:, 0],
                scene.lattice_points_arr[:, 1],
                scene.lattice_points_arr[:, 2],
                s=20,
                color="lightgray",
                alpha=0.3,
                depthshade=False,
                label="Tetrahedral lattice",
            )
        _draw_chain(ax_anim3d, scene.coords_arr, scene.bead_colors)
        _annotate_beads(
            ax_anim3d, scene.coords_arr, scene.coords, float(scene.max_range * 0.06)
        )
        ax_anim3d.set_xlabel("X")
        ax_anim3d.set_ylabel("Y")
        ax_anim3d.set_zlabel("Z")
        ax_anim3d.set_title("3D Protein Folding on Tetrahedral Lattice")
        ax_anim3d.view_init(elev=cfg.elev, azim=cfg.azim_start)
        plt.tight_layout()

        def _animate(frame: int) -> list[Any]:
            ax_anim3d.view_init(
                elev=cfg.elev,
                azim=cfg.azim_start + (360.0 * frame / cfg.frames),
            )
            return []

        anim = animation.FuncAnimation(
            fig_anim,
            _animate,
            frames=cfg.frames,
            interval=int(1000 / cfg.fps),
            blit=False,
        )

        path_str = str(dirpath / GIF_FILENAME)
        try:
            anim.save(path_str, writer="pillow", fps=cfg.fps)
        except Exception as e:
            logger.warning(f"Falling back to default writer for GIF due to: {e}")
            anim.save(path_str, fps=cfg.fps)
        logger.info(f"Saved rotating GIF to: {path_str}")
        plt.close(fig_anim)
    except Exception:
        logger.exception("Failed to create GIF")


def visualize_3d(
    coords: list[BeadPosition],
    *,
    dirpath: Path,
    figsize: tuple[int, int] = (8, 8),
    cmap: str = "viridis",
    show: bool = True,
) -> None:
    if len(coords) == 0:
        logger.warning("visualize_3d received empty coords; nothing to plot.")
        return

    coords_arr: NDArray[np.float64] = np.array(
        [bp.position for bp in coords], dtype=float
    )

    fig = cast(Any, plt).figure(figsize=figsize)
    ax: Axes3D = fig.add_subplot(111, projection="3d")  # type: ignore[assignment]
    ax3d = cast(Any, ax)

    base_dirs: NDArray[np.float64] = np.array(
        [[1.0, 1.0, 1.0], [1.0, -1.0, -1.0], [-1.0, 1.0, -1.0], [-1.0, -1.0, 1.0]],
        dtype=float,
    )
    base_dirs = base_dirs / float(np.linalg.norm(base_dirs[0]))

    max_range, mid = _compute_axes_limits(coords_arr)
    ax3d.set_xlim(float(mid[0] - max_range), float(mid[0] + max_range))
    ax3d.set_ylim(float(mid[1] - max_range), float(mid[1] + max_range))
    ax3d.set_zlim(float(mid[2] - max_range), float(mid[2] + max_range))

    limits = AxesLimits(
        x=cast(tuple[float, float], ax3d.get_xlim()),
        y=cast(tuple[float, float], ax3d.get_ylim()),
        z=cast(tuple[float, float], ax3d.get_zlim()),
    )
    lattice_points_arr = _build_lattice_points(
        base_dirs=base_dirs,
        max_range=max_range,
        center=coords_arr.mean(axis=0),
        limits=limits,
    )

    if lattice_points_arr.size > 0:
        ax3d.scatter(
            lattice_points_arr[:, 0],
            lattice_points_arr[:, 1],
            lattice_points_arr[:, 2],
            s=20,
            color="lightgray",
            alpha=0.3,
            depthshade=False,
            label="Tetrahedral lattice",
        )

    cmap_obj = plt.get_cmap(cmap)
    t = np.linspace(0.0, 1.0, coords_arr.shape[0])
    bead_colors = cmap_obj(t)

    _draw_chain(ax3d, coords_arr, bead_colors)
    _annotate_beads(ax3d, coords_arr, coords, float(max_range * 0.06))

    ax3d.set_xlabel("X")
    ax3d.set_ylabel("Y")
    ax3d.set_zlabel("Z")
    ax3d.set_title("3D Protein Folding on Tetrahedral Lattice")

    elev: float = 20.0
    azim_start: float = -60.0
    gif_frames: int = 120
    gif_fps: int = 20
    ax3d.view_init(elev=elev, azim=azim_start)
    ax3d.legend(loc="upper left")
    plt.tight_layout()

    scene = PlotScene(
        coords_arr=coords_arr,
        coords=coords,
        lattice_points_arr=lattice_points_arr,
        bead_colors=bead_colors,
        mid=mid,
        max_range=max_range,
    )
    cfg = RotGifConfig(
        figsize=figsize,
        elev=elev,
        azim_start=azim_start,
        frames=gif_frames,
        fps=gif_fps,
    )
    _save_rotating_gif(scene, cfg, dirpath=dirpath)

    if show:
        cast(Any, plt).show()
    else:
        plt.close(fig)
