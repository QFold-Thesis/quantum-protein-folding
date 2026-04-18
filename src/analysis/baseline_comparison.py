"""Baseline comparison workflows for particle studies.

Provides utilities to compare protein folding with and without particles,
establishing baseline behavior for sensitivity analysis.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

from logger import get_logger

if TYPE_CHECKING:
    from interaction.interaction import Interaction
    from particle import Particle
    from protein import Protein

logger = get_logger()


@dataclass
class BaselineComparison:
    """Results comparing system with and without particle.

    Attributes:
        energy_no_particle (float): Minimum energy without particle.
        energy_with_particle (float): Minimum energy with particle.
        energy_delta (float): Difference (with - without).
        config_no_particle (dict): Reference configuration without particle.
        config_with_particle (dict): Configuration with particle present.
        stabilization_metric (float): How much particle stabilizes (-) or destabilizes (+).

    """

    energy_no_particle: float
    energy_with_particle: float
    config_no_particle: dict
    config_with_particle: dict

    @property
    def energy_delta(self) -> float:
        """Change in energy due to particle presence."""
        return self.energy_with_particle - self.energy_no_particle

    @property
    def is_stabilizing(self) -> bool:
        """Whether particle lowers energy (is stabilizing)."""
        return self.energy_delta < 0

    @property
    def stabilization_factor(self) -> float:
        """Ratio of baseline to particle-modified energy."""
        return (self.energy_no_particle + 1e-10) / (self.energy_with_particle + 1e-10)


class BaselineComparisonWorkflow:
    """Workflow for systematic baseline comparisons.

    Attributes:
        protein (Protein): Protein to study.
        interaction_base (Interaction): Interaction model without particle support.
        particle (Particle): Particle for study.

    """

    def __init__(
        self,
        protein: Protein,
        interaction_base: Interaction,
        particle: Particle,
    ) -> None:
        """Initialize comparison workflow.

        Args:
            protein (Protein): Protein sequence to fold.
            interaction_base (Interaction): Standard interaction model.
            particle (Particle): Particle to study effects of.

        """
        self.protein = protein
        self.interaction_base = interaction_base
        self.particle = particle
        self.comparisons: dict[str, BaselineComparison] = {}

        logger.info(
            "Initialized baseline comparison workflow for %s with %s particle variant",
            protein.main_chain.get_sequence(),
            particle.particle_type.value,
        )

    def run_baseline_comparison(self) -> BaselineComparison:
        """Run complete baseline vs particle comparison.

        Finds minimum energy configurations for both cases and compares.

        Returns:
            BaselineComparison: Results of comparison.

        """
        logger.info("Running baseline comparison...")

        from validation.brute_force_folding import BruteForceFolding

        # Without particle
        solver_baseline = BruteForceFolding(
            protein=self.protein,
            interaction=self.interaction_base,
            particle=None,
            dimension=2,
        )
        result_baseline = solver_baseline.solve()

        # With particle
        solver_particle = BruteForceFolding(
            protein=self.protein,
            interaction=self.interaction_base,
            particle=self.particle,
            dimension=2,
        )
        result_particle = solver_particle.solve()

        comparison = BaselineComparison(
            energy_no_particle=result_baseline["energy"],
            energy_with_particle=result_particle["energy"],
            config_no_particle=result_baseline["configuration"],
            config_with_particle=result_particle["configuration"],
        )

        logger.info(
            "Baseline energy: %f, With particle: %f (delta: %f)",
            comparison.energy_no_particle,
            comparison.energy_with_particle,
            comparison.energy_delta,
        )

        self.comparisons["default"] = comparison
        return comparison

    def run_multiple_particle_positions(
        self, positions: list[tuple[int, int, int]]
    ) -> dict[tuple[int, int, int], BaselineComparison]:
        """Compare baselines for particle at different positions.

        Args:
            positions (list): List of (x, y, z) positions to test.

        Returns:
            dict: Results keyed by position.

        """
        logger.info("Running comparison for %d particle positions...", len(positions))

        results = {}

        for pos in positions:
            # Update particle position
            if hasattr(self.particle, "position"):
                self.particle.position = pos

            comparison = self.run_baseline_comparison()
            results[pos] = comparison

            logger.debug("Position %s: delta=%f", pos, comparison.energy_delta)

        return results

    def analyze_stabilization_pattern(self) -> dict:
        """Analyze pattern of energy changes across positions.

        Returns:
            dict: Statistical analysis of stabilization effects.

        """
        if not self.comparisons:
            msg = "Must run comparisons first"
            raise ValueError(msg)

        deltas = [c.energy_delta for c in self.comparisons.values()]
        min_delta = min(deltas)
        max_delta = max(deltas)
        avg_delta = sum(deltas) / len(deltas)

        num_stabilizing = sum(1 for d in deltas if d < 0)
        num_destabilizing = sum(1 for d in deltas if d > 0)

        logger.info(
            "Stabilization analysis: avg_delta=%f, stabilizing=%d, destabilizing=%d",
            avg_delta,
            num_stabilizing,
            num_destabilizing,
        )

        return {
            "mean_energy_delta": avg_delta,
            "min_energy_delta": min_delta,
            "max_energy_delta": max_delta,
            "num_positions": len(deltas),
            "stabilizing_count": num_stabilizing,
            "destabilizing_count": num_destabilizing,
            "stabilization_fraction": num_stabilizing / len(deltas),
        }

    def get_report(self) -> str:
        """Generate text report of baseline comparisons.

        Returns:
            str: Formatted report of results.

        """
        report = "\n" + "=" * 60 + "\n"
        report += "BASELINE COMPARISON REPORT\n"
        report += "=" * 60 + "\n\n"

        report += f"Protein sequence: {self.protein.main_chain.get_sequence()}\n"
        report += f"Particle type: {self.particle.particle_type.value}\n"
        report += f"Comparisons performed: {len(self.comparisons)}\n\n"

        for name, comparison in self.comparisons.items():
            report += f"Comparison: {name}\n"
            report += f"  Energy (no particle): {comparison.energy_no_particle:.6f}\n"
            report += f"  Energy (with particle): {comparison.energy_with_particle:.6f}\n"
            report += f"  Energy delta: {comparison.energy_delta:.6f}\n"
            report += f"  Stabilizing: {comparison.is_stabilizing}\n"
            report += f"  Stabilization factor: {comparison.stabilization_factor:.4f}\n"
            report += "\n"

        report += "=" * 60 + "\n"
        return report
