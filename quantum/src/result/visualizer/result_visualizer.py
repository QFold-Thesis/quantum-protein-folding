# ruff: noqa: I001
from __future__ import annotations

from typing import TYPE_CHECKING, Any, cast

if TYPE_CHECKING:
    from pathlib import Path

import numpy as np
from matplotlib import animation
from matplotlib import pyplot as plt

from constants import GIF_FILENAME, PLOT2D_FILENAME, QUBITS_PER_TURN
from logger import get_logger
from result.models import AxesLimits, BeadPosition, PlotScene, RotGifConfig

logger = get_logger()


class ResultVisualizer:
    """
    Create 3D and 2D visualizations of a protein conformation.

    - 3D: static window (no rotation) plus a rotating GIF saved to disk.
    - 2D: simple diagram derived directly from the final turn sequence.

    Configuration is intentionally minimal (no customization): only dirpath is
    required; other plot parameters are sensible defaults.
    """

    def __init__(self, dirpath: Path) -> None:
        self.dirpath: Path = dirpath

        self._figsize3d: tuple[int, int] = (8, 8)
        self._cmap3d: str = "viridis"
        self._elev: float = 20.0
        self._azim_start: float = -60.0
        self._gif_frames: int = 120
        self._gif_fps: int = 20

    def generate_3d(
        self, *, coords: list[BeadPosition], show: bool = False
    ) -> Path | None:
        """Create a static 3D plot and save a rotating GIF. Returns GIF path or None."""
        if not coords:
            logger.warning("generate_3d received empty coords; nothing to plot.")
            return None

        coords_arr = np.array([bp.position for bp in coords], dtype=float)
        fig = cast("Any", plt).figure(figsize=self._figsize3d)
        ax = fig.add_subplot(111, projection="3d")
        ax3d: Any = ax

        # Symmetric limits
        max_range, mid = self._compute_axes_limits(coords_arr)
        ax3d.set_xlim(float(mid[0] - max_range), float(mid[0] + max_range))
        ax3d.set_ylim(float(mid[1] - max_range), float(mid[1] + max_range))
        ax3d.set_zlim(float(mid[2] - max_range), float(mid[2] + max_range))

        # Lattice (tetrahedral directions), clipped to axes
        base_dirs = np.array(
            [[1.0, 1.0, 1.0], [1.0, -1.0, -1.0], [-1.0, 1.0, -1.0], [-1.0, -1.0, 1.0]],
            dtype=float,
        )
        base_dirs = base_dirs / float(np.linalg.norm(base_dirs[0]))
        limits = AxesLimits(
            x=cast("tuple[float, float]", ax3d.get_xlim()),
            y=cast("tuple[float, float]", ax3d.get_ylim()),
            z=cast("tuple[float, float]", ax3d.get_zlim()),
        )
        lattice_points_arr = self._build_lattice_points(
            base_dirs=base_dirs,
            max_range=max_range,
            center=coords_arr.mean(axis=0),
            limits=limits,
        )
        if lattice_points_arr.size > 0:
            ax3d.scatter(
                lattice_points_arr[:, 0],
                lattice_points_arr[:, 1],
                zs=lattice_points_arr[:, 2],
                s=20,
                color="lightgray",
                alpha=0.3,
                depthshade=False,
                label="Tetrahedral lattice",
            )

        # Chain colors and drawing
        cmap_obj = plt.get_cmap(self._cmap3d)
        t = np.linspace(0.0, 1.0, coords_arr.shape[0])
        bead_colors = cmap_obj(t)
        self._draw_chain(ax3d, coords_arr, bead_colors)
        self._annotate_beads(
            ax3d, coords_arr, coords=coords, offset=float(max_range * 0.06)
        )

        ax3d.set_xlabel("X")
        ax3d.set_ylabel("Y")
        ax3d.set_zlabel("Z")
        ax3d.set_title("3D Protein Folding on Tetrahedral Lattice")
        ax3d.view_init(elev=self._elev, azim=self._azim_start)
        ax3d.legend(loc="upper left")
        plt.tight_layout()

        # Save rotating GIF on a separate figure
        scene = PlotScene(
            coords_arr=coords_arr,
            coords=coords,
            lattice_points_arr=lattice_points_arr,
            bead_colors=bead_colors,
            mid=mid,
            max_range=max_range,
        )
        cfg = RotGifConfig(
            figsize=self._figsize3d,
            elev=self._elev,
            azim_start=self._azim_start,
            frames=self._gif_frames,
            fps=self._gif_fps,
        )
        gif_path = self._save_rotating_gif(scene, cfg)

        if show:
            cast("Any", plt).show()
        else:
            plt.close(fig)

        return gif_path

    def generate_2d(self, turn_sequence: str) -> Path:
        """Create and save a simple 2D diagram built directly from the turn sequence."""
        bits_per = int(QUBITS_PER_TURN)
        if len(turn_sequence) % bits_per != 0:
            usable_len = len(turn_sequence) - (len(turn_sequence) % bits_per)
            turns_bits = [
                turn_sequence[i : i + bits_per] for i in range(0, usable_len, bits_per)
            ]
        else:
            turns_bits = [
                turn_sequence[i : i + bits_per]
                for i in range(0, len(turn_sequence), bits_per)
            ]

        try:
            turns = [int(b, 2) for b in turns_bits]
        except ValueError as exc:
            msg = "turn_sequence must contain only '0' and '1' characters in groups of QUBITS_PER_TURN"
            raise ValueError(msg) from exc

        # Map turns to 2D steps (right, up, left, down)
        dir_map = np.array(
            [[1.0, 0.0], [0.0, 1.0], [-1.0, 0.0], [0.0, -1.0]], dtype=float
        )
        positions = [np.array([0.0, 0.0], dtype=float)]
        for t in turns:
            idx = int(t) % dir_map.shape[0]
            positions.append(positions[-1] + dir_map[idx])
        pts = np.vstack(positions)

        # Plot
        fig = cast("Any", plt).figure(figsize=(8, 6))
        ax = fig.add_subplot(111)

        # Connecting line
        min_points_for_line = 2
        if pts.shape[0] >= min_points_for_line:
            ax.plot(
                pts[:, 0],
                pts[:, 1],
                linestyle="-",
                color="#444444",
                linewidth=2,
                zorder=1,
            )
        # Beads (large circles)
        cmap2d = plt.get_cmap("tab10")
        colors = [cmap2d(i % 10) for i in range(pts.shape[0])]
        ax.scatter(
            pts[:, 0],
            pts[:, 1],
            s=300,
            c=colors,
            edgecolors="black",
            linewidths=0.8,
            zorder=2,
        )
        # Labels inside circles
        for i, (x, y) in enumerate(pts):
            ax.text(
                float(x),
                float(y),
                str(i + 1),
                ha="center",
                va="center",
                fontsize=10,
                fontweight="bold",
                color="white",
                zorder=3,
            )

        # Fit view with padding to prevent cropping
        x_min, x_max = float(pts[:, 0].min()), float(pts[:, 0].max())
        y_min, y_max = float(pts[:, 1].min()), float(pts[:, 1].max())
        x_span = x_max - x_min if x_max > x_min else 1.0
        y_span = y_max - y_min if y_max > y_min else 1.0
        pad = 0.25 * max(x_span, y_span)
        ax.set_xlim(x_min - pad, x_max + pad)
        ax.set_ylim(y_min - pad, y_max + pad)
        ax.set_aspect("equal", adjustable="box")
        ax.axis("off")

        out_path = self.dirpath / PLOT2D_FILENAME
        try:
            fig.tight_layout()
            fig.savefig(out_path, dpi=200, bbox_inches="tight")
            logger.info("Saved 2D diagram to: %s", out_path)
        except Exception:
            logger.exception("Failed to save 2D diagram")
        finally:
            plt.close(fig)
        return out_path

    def generate_plots(
        self, *, coords: list[BeadPosition], turn_sequence: str, show_3d: bool = False
    ) -> dict[str, Path | None]:
        """Convenience method to create both outputs and return their paths."""
        gif_path = self.generate_3d(coords=coords, show=show_3d)
        png_path = self.generate_2d(turn_sequence)
        return {"3d_gif": gif_path, "2d_png": png_path}

    @staticmethod
    def _compute_axes_limits(coords_arr: np.ndarray) -> tuple[float, np.ndarray]:
        max_range: float = float(
            (coords_arr.max(axis=0) - coords_arr.min(axis=0)).max() / 2.0
        )
        mid: np.ndarray = coords_arr.mean(axis=0)
        return max_range, mid

    @staticmethod
    def _build_lattice_points(
        *,
        base_dirs: np.ndarray,
        max_range: float,
        center: np.ndarray,
        limits: AxesLimits,
    ) -> np.ndarray:
        lattice_points: list[np.ndarray] = []
        n: int = max(2, int(np.ceil(max_range)) + 2)
        for i in range(-n, n + 1):
            for j in range(-n, n + 1):
                for k in range(-n, n + 1):
                    point: np.ndarray = (
                        i * base_dirs[0] + j * base_dirs[1] + k * base_dirs[2]
                    )
                    lattice_points.append(point)
        lattice_points_arr: np.ndarray = np.unique(
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

    @staticmethod
    def _draw_chain(ax3d: Any, coords_arr: np.ndarray, colors: np.ndarray) -> None:
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

    @staticmethod
    def _annotate_beads(
        ax3d: Any, coords_arr: np.ndarray, coords: list[BeadPosition], offset: float
    ) -> None:
        for i, (x, y, z) in enumerate(coords_arr):
            ax3d.text(
                float(x + offset),
                float(y + offset),
                float(z + offset),
                coords[i].symbol,
                color="black",
                fontsize=9,
                ha="center",
                va="center",
            )

    def _save_rotating_gif(self, scene: PlotScene, cfg: RotGifConfig) -> Path | None:
        """Render GIF on a separate figure so the shown window stays static."""
        try:
            fig_anim = cast("Any", plt).figure(figsize=cfg.figsize)
            ax_anim = fig_anim.add_subplot(111, projection="3d")
            ax_anim3d: Any = ax_anim

            ax_anim3d.set_xlim(
                float(scene.mid[0] - scene.max_range),
                float(scene.mid[0] + scene.max_range),
            )
            ax_anim3d.set_ylim(
                float(scene.mid[1] - scene.max_range),
                float(scene.mid[1] + scene.max_range),
            )
            ax_anim3d.set_zlim(
                float(scene.mid[2] - scene.max_range),
                float(scene.mid[2] + scene.max_range),
            )
            if scene.lattice_points_arr.size > 0:
                ax_anim3d.scatter(
                    scene.lattice_points_arr[:, 0],
                    scene.lattice_points_arr[:, 1],
                    zs=scene.lattice_points_arr[:, 2],
                    s=20,
                    color="lightgray",
                    alpha=0.3,
                    depthshade=False,
                    label="Tetrahedral lattice",
                )
            self._draw_chain(ax_anim3d, scene.coords_arr, scene.bead_colors)
            self._annotate_beads(
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
                    elev=cfg.elev, azim=cfg.azim_start + (360.0 * frame / cfg.frames)
                )
                return []

            anim = animation.FuncAnimation(
                fig_anim,
                _animate,
                frames=cfg.frames,
                interval=int(1000 / cfg.fps),
                blit=False,
            )

            path = self.dirpath / GIF_FILENAME
            path_str = str(path)
            try:
                anim.save(path_str, writer="pillow", fps=cfg.fps)
            except Exception as e:
                logger.warning("Falling back to default writer for GIF due to: %s", e)
                anim.save(path_str, fps=cfg.fps)
            logger.info("Saved rotating GIF to: %s", path_str)
            plt.close(fig_anim)
        except Exception:
            logger.exception("Failed to create GIF")
            return None
        else:
            return path
