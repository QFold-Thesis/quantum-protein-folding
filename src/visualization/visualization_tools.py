"""Visualization tools for particle folding studies.

Creates plots and visualizations for energy landscapes, conformations,
sensitivity analyses, and genetic algorithm convergence.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from logger import get_logger

if TYPE_CHECKING:
    from analysis.baseline_comparison import BaselineComparison
    from optimization.sequence_design_ga import SequenceDesignGA

logger = get_logger()


class VisualizationTools:
    """Tools for creating publication-quality visualizations.

    Methods cover: energy landscapes, conformations, parameter sensitivity,
    convergence plots, and comparative analyses.

    """

    @staticmethod
    def plot_energy_landscape_comparison(
        energy_no_particle: list[float],
        energy_with_particle: list[float],
        output_path: str = "energy_comparison.png",
    ) -> None:
        """Plot energy distributions with and without particle.

        Args:
            energy_no_particle: Energy values of configurations without particle.
            energy_with_particle: Energy values of configurations with particle.
            output_path: Where to save figure.

        """
        try:
            import matplotlib.pyplot as plt
            import numpy as np

            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

            # Histograms
            bins = 20
            ax1.hist(energy_no_particle, bins=bins, alpha=0.7, label="No particle", color="blue")
            ax1.hist(energy_with_particle, bins=bins, alpha=0.7, label="With particle", color="red")
            ax1.set_xlabel("Energy")
            ax1.set_ylabel("Frequency")
            ax1.set_title("Energy Distribution Comparison")
            ax1.legend()
            ax1.grid(alpha=0.3)

            # Box plot
            ax2.boxplot([energy_no_particle, energy_with_particle])
            ax2.set_xticklabels(["No particle", "With particle"])
            ax2.set_ylabel("Energy")
            ax2.set_title("Energy Statistics")
            ax2.grid(alpha=0.3)

            plt.tight_layout()
            plt.savefig(output_path, dpi=300)
            logger.info("Saved energy comparison plot to %s", output_path)
            plt.close()

        except ImportError:
            logger.warning("matplotlib not available for visualization")

    @staticmethod
    def plot_parameter_sensitivity(
        param_values: list[float],
        metrics: list[float],
        param_name: str = "Parameter",
        output_path: str = "sensitivity.png",
    ) -> None:
        """Plot parameter sensitivity curve.

        Args:
            param_values: Parameter values tested.
            metrics: Corresponding metric values.
            param_name: Name of parameter.
            output_path: Where to save figure.

        """
        try:
            import matplotlib.pyplot as plt

            fig, ax = plt.subplots(figsize=(10, 6))

            ax.plot(param_values, metrics, "o-", linewidth=2, markersize=8, color="darkblue")
            ax.set_xlabel(f"{param_name} Value")
            ax.set_ylabel("Metric Value")
            ax.set_title(f"Sensitivity to {param_name}")
            ax.grid(alpha=0.3)

            plt.tight_layout()
            plt.savefig(output_path, dpi=300)
            logger.info("Saved sensitivity plot to %s", output_path)
            plt.close()

        except ImportError:
            logger.warning("matplotlib not available for visualization")

    @staticmethod
    def plot_ga_convergence(
        ga: SequenceDesignGA,
        output_path: str = "ga_convergence.png",
    ) -> None:
        """Plot genetic algorithm convergence.

        Args:
            ga: SequenceDesignGA after run_evolution().
            output_path: Where to save figure.

        """
        try:
            import matplotlib.pyplot as plt

            generations, fitnesses = ga.get_convergence_data()

            fig, ax = plt.subplots(figsize=(10, 6))

            ax.plot(generations, fitnesses, "o-", linewidth=2, color="green")
            ax.set_xlabel("Generation")
            ax.set_ylabel("Best Fitness")
            ax.set_title("Genetic Algorithm Convergence")
            ax.grid(alpha=0.3)

            plt.tight_layout()
            plt.savefig(output_path, dpi=300)
            logger.info("Saved GA convergence plot to %s", output_path)
            plt.close()

        except ImportError:
            logger.warning("matplotlib not available for visualization")

    @staticmethod
    def plot_baseline_comparison_bar(
        comparisons: list[BaselineComparison],
        labels: list[str],
        output_path: str = "baseline_comparison.png",
    ) -> None:
        """Create bar plot comparing baseline energies.

        Args:
            comparisons: List of BaselineComparison objects.
            labels: Labels for each comparison.
            output_path: Where to save figure.

        """
        try:
            import matplotlib.pyplot as plt
            import numpy as np

            energies_no_particle = [c.energy_no_particle for c in comparisons]
            energies_with_particle = [c.energy_with_particle for c in comparisons]

            x = np.arange(len(labels))
            width = 0.35

            fig, ax = plt.subplots(figsize=(12, 6))

            bars1 = ax.bar(x - width / 2, energies_no_particle, width, label="No particle", color="blue")
            bars2 = ax.bar(x + width / 2, energies_with_particle, width, label="With particle", color="red")

            ax.set_xlabel("Condition")
            ax.set_ylabel("Energy")
            ax.set_title("Baseline Comparison")
            ax.set_xticks(x)
            ax.set_xticklabels(labels)
            ax.legend()
            ax.grid(alpha=0.3, axis="y")

            # Add value labels on bars
            for bars in [bars1, bars2]:
                for bar in bars:
                    height = bar.get_height()
                    ax.text(
                        bar.get_x() + bar.get_width() / 2,
                        height,
                        f"{height:.2f}",
                        ha="center",
                        va="bottom",
                        fontsize=9,
                    )

            plt.tight_layout()
            plt.savefig(output_path, dpi=300)
            logger.info("Saved baseline comparison plot to %s", output_path)
            plt.close()

        except ImportError:
            logger.warning("matplotlib not available for visualization")

    @staticmethod
    def plot_2d_conformation(
        configuration: dict,
        sequence: str = "",
        output_path: str = "conformation.png",
    ) -> None:
        """Plot 2D protein conformation.

        Args:
            configuration: Mapping of bead indices to (x, y, z) coordinates.
            sequence: Amino acid sequence (for coloring).
            output_path: Where to save figure.

        """
        try:
            import matplotlib.pyplot as plt
            import matplotlib.patches as patches

            fig, ax = plt.subplots(figsize=(10, 10))

            coords = [configuration.get(i, (0, 0, 0)) for i in range(len(configuration))]

            # Plot backbone
            xs = [c[0] for c in coords]
            ys = [c[1] for c in coords]
            ax.plot(xs, ys, "k-", linewidth=2, alpha=0.5)

            # Plot beads
            colors = ["red" if char in "FLI" else "blue" for char in sequence] if sequence else ["gray"] * len(coords)

            for i, (x, y) in enumerate(zip(xs, ys)):
                circle = patches.Circle((x, y), 0.3, color=colors[i] if i < len(colors) else "gray", ec="black")
                ax.add_patch(circle)
                ax.text(x, y, str(i), ha="center", va="center", fontsize=8, color="white", weight="bold")

            ax.set_aspect("equal")
            ax.set_xlabel("X")
            ax.set_ylabel("Y")
            ax.set_title("2D Protein Conformation")
            ax.grid(alpha=0.3)

            plt.tight_layout()
            plt.savefig(output_path, dpi=300)
            logger.info("Saved conformation plot to %s", output_path)
            plt.close()

        except ImportError:
            logger.warning("matplotlib not available for visualization")

    @staticmethod
    def plot_heatmap_2d_sweep(
        hh_values: list[float],
        pp_values: list[float],
        energy_matrix: list[list[float]],
        output_path: str = "sweep_heatmap.png",
    ) -> None:
        """Create heatmap of 2D parameter sweep.

        Args:
            hh_values: Particle-H interaction values.
            pp_values: Particle-P interaction values.
            energy_matrix: 2D matrix of energy values.
            output_path: Where to save figure.

        """
        try:
            import matplotlib.pyplot as plt
            import numpy as np

            fig, ax = plt.subplots(figsize=(10, 8))

            im = ax.imshow(
                energy_matrix,
                cmap="RdYlBu_r",
                aspect="auto",
                extent=[min(pp_values), max(pp_values), min(hh_values), max(hh_values)],
                origin="lower",
            )

            ax.set_xlabel("Particle-Polar Interaction")
            ax.set_ylabel("Particle-Hydrophobic Interaction")
            ax.set_title("Energy Landscape: 2D Parameter Sweep")

            cbar = plt.colorbar(im, ax=ax)
            cbar.set_label("Energy")

            plt.tight_layout()
            plt.savefig(output_path, dpi=300)
            logger.info("Saved heatmap to %s", output_path)
            plt.close()

        except ImportError:
            logger.warning("matplotlib not available for visualization")
