"""Utilities for visualizing protein folding results.

This module provides the `ResultVisualizer` class, which processes
the interpreted results of protein folding simulations and generates 3D and 2D visualizations of
the folded protein structure.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from constants import (
    CONFORMATION_ENCODING,
    FCC_EDGE_LENGTH,
    FLAT_VISUALIZATION_FILENAME,
    GIF_VISUALIZATION_FILENAME,
    HTML_VISUALIZATION_FILENAME,
    INTERACTION_TYPE,
    QUBITS_PER_TURN,
    TETRAHEDRAL_LATTICE_PADDING,
)
from logger.logger import get_logger

if TYPE_CHECKING:
    from pathlib import Path

    from matplotlib.colors import Colormap
    from matplotlib.figure import Figure
    from mpl_toolkits.mplot3d.axes3d import Axes3D
    from numpy.typing import NDArray

    from enums import TurnDirection
    from result.models import BeadPosition

logger = get_logger()


class ResultVisualizer:
    """Visualizes and processes interpreted results of quantum protein folding simulations."""

    def __init__(
        self,
        dirpath: Path,
        turn_sequence: list[TurnDirection],
        coordinates_3d: list[BeadPosition],
        main_main_contacts_detected: dict[int, int],
    ) -> None:
        """Initialize the ResultVisualizer with contacts data from the main chain, directory path and decoded turn sequence and 3D coordinates.

        Args:
            dirpath (Path): Directory path for saving output files.
            turn_sequence (list[TurnDirection]): Decoded turn sequence.
            coordinates_3d (list[BeadPosition]): 3D coordinates of the protein structure.
            main_main_contacts_detected (dict[int, int]): Detected contacts between main chain beads.

        """
        self._dirpath: Path = dirpath
        self._turn_sequence: list[TurnDirection] = turn_sequence
        self._coordinates_3d: list[BeadPosition] = coordinates_3d
        self._main_main_contacts_detected: dict[int, int] = main_main_contacts_detected

    def visualize_3d(self, filename: str = HTML_VISUALIZATION_FILENAME) -> None:
        """Generate interactive 3D visualization of the resulting conformation in the .html file format.

        Args:
            filename (str): The name of the file to save the interactive .html visualization. Defaults to HTML_VISUALIZATION_FILENAME.

        """
        import plotly.graph_objects as go
        from matplotlib import cm

        logger.debug(
            "Generating interactive 3D HTML visualization of the conformation..."
        )

        coords: NDArray[np.float64] = np.array(
            [(b.x, b.y, b.z) for b in self._coordinates_3d]
        )
        symbols: list[str] = [b.symbol for b in self._coordinates_3d]

        cmap: Colormap = cm.get_cmap("hsv", len(coords))
        colors: list[str] = [
            f"rgb({int(r * 255)}, {int(g * 255)}, {int(b * 255)})"
            for r, g, b, _ in cmap(np.linspace(0, 1, len(coords)))
        ]

        lattice_points: NDArray[np.float64] = self._generate_lattice_points(coords)
        fig: go.Figure = go.Figure()

        # Add tetrahedral lattice points
        fig.add_trace(
            go.Scatter3d(
                x=lattice_points[:, 0],
                y=lattice_points[:, 1],
                z=lattice_points[:, 2],
                mode="markers",
                marker={"size": 5, "color": "lightgray", "opacity": 0.45},
                name="Tetrahedral lattice points",
                legendgroup="lattice",
                showlegend=True,
                hoverinfo="skip",
            )
        )

        # Add contact lines
        for i, j in self._main_main_contacts_detected.items():
            x_vals = [coords[i, 0], coords[j, 0]]
            y_vals = [coords[i, 1], coords[j, 1]]
            z_vals = [coords[i, 2], coords[j, 2]]
            fig.add_trace(
                go.Scatter3d(
                    x=x_vals,
                    y=y_vals,
                    z=z_vals,
                    mode="lines",
                    line={"color": "purple", "width": 8, "dash": "dot"},
                    name=f"Contact indicator {i}-{j}",
                    showlegend=True,
                    hoverinfo="text",
                    hovertext=f"<b>Contact between beads {i} and {j}</b>",
                )
            )

        # Add protein bonds (turns)
        fig.add_trace(
            go.Scatter3d(
                x=coords[:, 0],
                y=coords[:, 1],
                z=coords[:, 2],
                mode="lines",
                line={"color": "black", "width": 5},
                name="Protein bonds",
                legendgroup="bonds",
                showlegend=True,
                hoverinfo="skip",
            )
        )

        # Add labels to the bonds
        for i in range(len(coords) - 1):
            x_mid = (coords[i, 0] + coords[i + 1, 0]) / 2
            y_mid = (coords[i, 1] + coords[i + 1, 1]) / 2
            z_mid = (coords[i, 2] + coords[i + 1, 2]) / 2
            turn_type: str = str(self._turn_sequence[i])

            fig.add_trace(
                go.Scatter3d(
                    x=[x_mid],
                    y=[y_mid],
                    z=[z_mid],
                    mode="markers+text",
                    marker={"size": 20, "color": "lightgray", "opacity": 0.6},
                    text=[turn_type],
                    textposition="middle center",
                    textfont={"family": "Arial", "size": 14, "color": "black"},
                    showlegend=(i == 0),
                    legendgroup="turns",
                    name="Turn indicators",
                    hoverinfo="text",
                    hovertext=f"<b>Turn between Beads {i} and {i + 1} - {turn_type}<b>",
                )
            )

        # Add protein beads
        for i, (sym, (x, y, z), color) in enumerate(
            zip(symbols, coords, colors, strict=True)
        ):
            text_color: str = self._get_text_color(color)
            fig.add_trace(
                go.Scatter3d(
                    x=[x],
                    y=[y],
                    z=[z],
                    mode="markers+text",
                    marker={
                        "size": 25,
                        "color": color,
                        "line": {"width": 1, "color": "black"},
                    },
                    text=f"{sym} ({i})",
                    textposition="middle center",
                    textfont={"family": "Arial", "size": 16, "color": text_color},
                    name=f"{sym} (Index {i})",
                    legendgroup=f"{sym}_{i}",
                    showlegend=True,
                    hoverinfo="text",
                    hovertext=f"<b>Bead {sym} (Index: {i})</b><br>Position: ({x:.2f}, {y:.2f}, {z:.2f})",
                )
            )

        # Add overlay features
        fig.update_layout(
            title=f"3D Protein Folding Visualization for main chain sequence: {''.join(symbols)}<br><br>Encoding: {CONFORMATION_ENCODING.name} (Qubits per turn: {QUBITS_PER_TURN})<br>Interaction model: {INTERACTION_TYPE.name}",
            scene={
                "xaxis_title": "X",
                "yaxis_title": "Y",
                "zaxis_title": "Z",
                "aspectmode": "data",
                "xaxis": {"showbackground": False},
                "yaxis": {"showbackground": False},
                "zaxis": {"showbackground": False},
            },
            legend_title_text="Legend",
            margin={"l": 0, "r": 0, "b": 0, "t": 40},
            scene_camera={"projection": {"type": "orthographic"}},
        )

        html_path: Path = self._dirpath / filename
        fig.write_html(html_path, include_plotlyjs=True, auto_open=False)
        logger.info("Interactive 3D HTML visualization saved to: %s", html_path)

    @staticmethod
    def _get_text_color(rgb_color: str, brightness_threshhold: float = 0.5) -> str:
        """Helper function to choose black or white text depending on the brightness of the RGB color (using luminance formula for RGB).

        Args:
            rgb_color (str): The RGB color string in the format "rgb(r, g, b)".
            brightness_threshhold (float): The brightness threshold for determining text color.

        """
        import re

        r, g, b = map(int, re.findall(r"\d+", rgb_color))

        rgb_color_values: int = 255

        r, g, b = r / rgb_color_values, g / rgb_color_values, b / rgb_color_values

        luminance_red_factor: float = 0.299
        luminance_green_factor: float = 0.587
        luminance_blue_factor: float = 0.114

        brightness: float = (
            luminance_red_factor * r
            + luminance_green_factor * g
            + luminance_blue_factor * b
        )
        return "white" if brightness < brightness_threshhold else "black"

    def generate_3d_gif(self, filename: str = GIF_VISUALIZATION_FILENAME) -> None:
        """Generate interactive simplified 3D visualization of the resulting conformation as a rotating plot in the .html file format.

        Args:
            filename (str): The name of the file to save the rotating .gif plot. Defaults to GIF_VISUALIZATION_FILENAME constant.

        """
        import matplotlib.pyplot as plt
        import numpy as np
        from matplotlib import animation

        logger.debug("Generating 3D rotating GIF visualization of the conformation...")

        coords: NDArray[np.float64] = np.array(
            [(b.x, b.y, b.z) for b in self._coordinates_3d]
        )
        symbols: list[str] = [b.symbol for b in self._coordinates_3d]

        cmap: Colormap = plt.get_cmap("hsv", len(coords))
        colors: list[tuple[float, float, float, float]] = [
            cmap(i) for i in range(len(coords))
        ]

        lattice: NDArray[np.float64] = self._generate_lattice_points(coords)

        fig: Figure = plt.figure(figsize=(12, 10))
        ax: Axes3D = fig.add_subplot(111, projection="3d")
        ax.set_facecolor("white")

        ax.scatter(
            lattice[:, 0],
            lattice[:, 1],
            lattice[:, 2],
            c="gray",
            alpha=0.2,
            s=5,
            label="Lattice points",
        )

        ax.plot(
            coords[:, 0],
            coords[:, 1],
            coords[:, 2],
            c="black",
            lw=2,
            label="Protein bonds",
        )

        scatter_handles: list[Axes3D.scatter] = []
        for i, (x, y, z, sym, color) in enumerate(
            zip(coords[:, 0], coords[:, 1], coords[:, 2], symbols, colors, strict=True)
        ):
            sc: Axes3D.scatter = ax.scatter(
                x, y, z, c=[color], s=90, edgecolors="black", label=f"{sym} (Index {i})"
            )
            scatter_handles.append(sc)

        ax.set_title(
            f"3D Protein Folding Visualization for main chain sequence: {''.join(symbols)}",
            fontsize=14,
        )
        ax.set_xlabel("X")
        ax.set_ylabel("Y")
        ax.set_zlabel("Z")
        ax.legend(loc="upper left", bbox_to_anchor=(1.05, 1))

        def __update_animation(angle: float) -> list[Axes3D.scatter]:
            """Helper function handler for FuncAnimation."""
            ax.view_init(elev=30, azim=angle)
            return scatter_handles

        angles: NDArray[np.float64] = np.linspace(0, 360, 72, endpoint=False)
        ani = animation.FuncAnimation(
            fig, __update_animation, frames=angles, blit=False
        )

        gif_path: Path = self._dirpath / filename
        writer = animation.PillowWriter(fps=10)
        ani.save(gif_path, writer=writer)

        plt.close(fig)
        logger.info("3D rotating GIF visualization saved to: %s", gif_path)

    def visualize_2d(self, filename: str = FLAT_VISUALIZATION_FILENAME) -> None:
        """Generate flat, 2D visualization of the resulting conformation in the .png file format.

        This plot is a 2D projection (top-down view) of the 3D coordinates.

        Args:
            filename (str): The name of the file to save the static .png visualization. Defaults to FLAT_VISUALIZATION_FILENAME.

        """
        import matplotlib.pyplot as plt
        import numpy as np
        from matplotlib.lines import Line2D

        logger.debug("Generating 2D flat visualization of the conformation...")

        symbols: list[str] = [b.symbol for b in self._coordinates_3d]
        contacts: dict[int, int] = self._main_main_contacts_detected

        coords_3d: NDArray[np.float64] = np.array(
            [(b.x, b.y, b.z) for b in self._coordinates_3d]
        )
        x_coords: NDArray[np.float64] = coords_3d[:, 0]
        y_coords: NDArray[np.float64] = coords_3d[:, 1]
        z_coords: NDArray[np.float64] = coords_3d[:, 2]

        fig: Figure = plt.figure(figsize=(10, 8))
        ax: Axes3D = fig.add_subplot(111, projection="3d")

        node_color: str = "tab:blue"
        node_size: int = 1000

        ax.plot(
            x_coords,
            y_coords,
            z_coords,
            color=node_color,
            lw=2.5,
            zorder=1,
        )

        for i, j in contacts.items():
            x_vals: list[float] = [coords_3d[i, 0], coords_3d[j, 0]]
            y_vals: list[float] = [coords_3d[i, 1], coords_3d[j, 1]]
            z_vals: list[float] = [coords_3d[i, 2], coords_3d[j, 2]]

            ax.plot(
                x_vals,
                y_vals,
                z_vals,
                color="black",
                lw=2.0,
                linestyle="--",
                zorder=1,
            )

        ax.scatter(
            x_coords,
            y_coords,
            z_coords,
            s=node_size,
            color=node_color,
            edgecolors="black",
            linewidth=0.5,
            zorder=2,
        )

        for i, (x, y, z) in enumerate(coords_3d):
            ax.text(
                x,
                y,
                z,
                symbols[i],
                color="white",
                ha="center",
                va="center",
                fontsize=10,
                fontweight="bold",
                zorder=3,
            )

        v: NDArray[np.int64] = np.array([-6, -4, -1])

        azim: np.float64 = np.rad2deg(np.arctan2(v[0], -v[1]))

        elev: np.float64 = np.rad2deg(np.arctan2(v[2], np.linalg.norm(v[:2])))

        ax.view_init(elev=elev, azim=azim)

        legend_elements: list[Line2D] = [
            Line2D(
                [0],
                [0],
                marker="o",
                color="w",
                markerfacecolor=node_color,
                markeredgecolor="black",
                markersize=12,
                label="Beads",
            )
        ]

        if contacts:
            legend_elements.append(
                Line2D(
                    [0],
                    [0],
                    color="black",
                    lw=2,
                    linestyle="--",
                    label="Interactions detected",
                )
            )

        ax.legend(
            handles=legend_elements,
            loc="lower center",
            bbox_to_anchor=(0.5, -0.05),
            ncol=2,
            fontsize=12,
            frameon=False,
        )
        ax.set_title(
            f"2D Protein Folding Visualization for main chain sequence: {''.join(symbols)}",
            fontsize=14,
        )

        ax.axis("off")
        plt.tight_layout()

        filepath: Path = self._dirpath / filename
        plt.savefig(filepath, format="png", bbox_inches="tight", dpi=150)
        plt.close(fig)
        logger.info("2D flat visualization saved to: %s", filepath)

    def _generate_lattice_points(
        self, coords: NDArray[np.float64], padding: int = TETRAHEDRAL_LATTICE_PADDING
    ) -> NDArray[np.float64]:
        """Generate tetrahedral lattice points in range that matches the coordinates of the conformation (by taking the min/max values from the coordinates and padding them accordingly).

        Args:
            coords (NDArray[np.float64]): The coordinates of the conformation.
            padding (int): Value that determines the padding around the coordinates (how big should the space be). Defaults to TETRAHEDRAL_LATTICE_PADDING.

        """
        min_x, max_x = coords[:, 0].min(), coords[:, 0].max()
        min_y, max_y = coords[:, 1].min(), coords[:, 1].max()
        min_z, max_z = coords[:, 2].min(), coords[:, 2].max()
        padding: int = 1

        x_range: NDArray[np.float64] = np.arange(
            min_x - padding * FCC_EDGE_LENGTH,
            max_x + padding * FCC_EDGE_LENGTH,
            FCC_EDGE_LENGTH,
        )
        y_range: NDArray[np.float64] = np.arange(
            min_y - padding * FCC_EDGE_LENGTH,
            max_y + padding * FCC_EDGE_LENGTH,
            FCC_EDGE_LENGTH,
        )
        z_range: NDArray[np.float64] = np.arange(
            min_z - padding * FCC_EDGE_LENGTH,
            max_z + padding * FCC_EDGE_LENGTH,
            FCC_EDGE_LENGTH,
        )

        xv, yv, zv = np.meshgrid(x_range, y_range, z_range, indexing="ij")
        return np.vstack([xv.ravel(), yv.ravel(), zv.ravel()]).T
