"""Sensitivity analysis framework for particle model parameters.

Performs systematic parameter sweeps and bifurcation analysis to characterize
how particle interaction strengths and configurations affect folding behavior.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Callable

from logger import get_logger

if TYPE_CHECKING:
    from interaction.interaction import Interaction
    from particle import Particle
    from protein import Protein

logger = get_logger()


class SensitivityAnalysisFramework:
    """Framework for systematic parameter sensitivity studies.

    Attributes:
        protein (Protein): Protein being studied.
        interaction (Interaction): Interaction model.
        particle (Particle): Particle with tunable parameters.

    """

    def __init__(
        self,
        protein: Protein,
        interaction: Interaction,
        particle: Particle,
    ) -> None:
        """Initialize sensitivity analysis.

        Args:
            protein (Protein): Protein sequence.
            interaction (Interaction): Interaction model.
            particle (Particle): Particle for sensitivity analysis.

        """
        self.protein = protein
        self.interaction = interaction
        self.particle = particle
        self.scan_results: dict[str, dict] = {}

        logger.info(
            "Initialized sensitivity analysis for %s particle type",
            particle.particle_type.value,
        )

    def sweep_single_parameter(
        self,
        param_name: str,
        param_values: list[float],
        metric: Callable | None = None,
    ) -> dict:
        """Sweep a single particle/interaction parameter.

        Args:
            param_name (str): Parameter name to vary.
            param_values (list[float]): Values to test.
            metric (Callable, optional): Function to evaluate at each value.
                Should take (protein, interaction, particle) and return float.

        Returns:
            dict: Results mapping parameter values to metric results.

        """
        logger.info(
            "Starting sweep of parameter '%s' with %d values", param_name, len(param_values)
        )

        if metric is None:
            metric = self._default_metric

        results = {
            "parameter": param_name,
            "values": param_values,
            "metrics": [],
        }

        for value in param_values:
            # Update parameter
            self._set_parameter(param_name, value)

            # Evaluate metric
            try:
                metric_value = metric(self.protein, self.interaction, self.particle)
                results["metrics"].append(metric_value)
                logger.debug("Parameter %s = %f: metric = %f", param_name, value, metric_value)
            except Exception as e:
                logger.warning("Error evaluating at %s = %f: %s", param_name, value, e)
                results["metrics"].append(None)

        self.scan_results[param_name] = results
        logger.info("Parameter sweep complete for '%s'", param_name)

        return results

    def sweep_hp_interaction_strength(
        self, hh_values: list[float], pp_values: list[float]
    ) -> dict:
        """Sweep HP model particle interaction strengths.

        For HP model, sweep particle-hydrophobic and particle-polar interactions.

        Args:
            hh_values (list[float]): Particle-hydrophobic interaction energies to test.
            pp_values (list[float]): Particle-polar interaction energies to test.

        Returns:
            dict: 2D sweep results.

        """
        from interaction.hp_interaction_with_particle import HPInteractionWithParticle

        if not isinstance(self.interaction, HPInteractionWithParticle):
            msg = "Interaction model must be HPInteractionWithParticle"
            raise TypeError(msg)

        logger.info(
            "Sweeping HP particle interactions: %d HH values, %d PP values",
            len(hh_values),
            len(pp_values),
        )

        results = {
            "hh_values": hh_values,
            "pp_values": pp_values,
            "energy_matrix": [],
        }

        for hh in hh_values:
            row = []
            for pp in pp_values:
                self.interaction.set_particle_interaction(
                    particle_hh_energy=hh, particle_pp_energy=pp
                )

                metric = self._default_metric(
                    self.protein, self.interaction, self.particle
                )
                row.append(metric)

            results["energy_matrix"].append(row)

        self.scan_results["hp_2d_sweep"] = results
        logger.info("HP 2D sweep complete: %d x %d matrix", len(hh_values), len(pp_values))

        return results

    def identify_critical_parameters(self, threshold: float = 0.1) -> dict:
        """Identify parameters with large sensitivity.

        Analyzes scan results to find parameters where small changes
        cause significant metric changes.

        Args:
            threshold (float, optional): Fractional change threshold. Defaults to 0.1.

        Returns:
            dict: Critical parameters and their sensitivities.

        """
        logger.info("Identifying critical parameters (threshold=%f)", threshold)

        critical = {}

        for param_name, results in self.scan_results.items():
            if param_name.endswith("_2d_sweep"):
                continue

            metrics = [m for m in results["metrics"] if m is not None]
            if not metrics or len(metrics) < 2:
                continue

            # Calculate sensitivity: max slope in parameter space
            max_slope = 0.0
            for i in range(len(metrics) - 1):
                dm = abs(metrics[i + 1] - metrics[i])
                max_slope = max(max_slope, dm)

            # Normalize by metric magnitude
            avg_metric = sum(metrics) / len(metrics)
            if avg_metric != 0:
                sensitivity = max_slope / abs(avg_metric)
            else:
                sensitivity = max_slope

            if sensitivity > threshold:
                critical[param_name] = {
                    "sensitivity": sensitivity,
                    "max_slope": max_slope,
                }

        logger.info("Found %d critical parameters", len(critical))
        return critical

    def analyze_bifurcation_points(self) -> dict:
        """Analyze bifurcation points in parameter space.

        Detects where optimal configurations change qualitatively.

        Returns:
            dict: Bifurcation analysis across all swept parameters.

        """
        logger.info("Analyzing bifurcation points...")

        from validation.brute_force_folding import BruteForceFolding

        bifurcations = {}

        for param_name, results in self.scan_results.items():
            min_configs = []

            for value in results["values"]:
                self._set_parameter(param_name, value)

                solver = BruteForceFolding(
                    protein=self.protein,
                    interaction=self.interaction,
                    particle=self.particle,
                    dimension=2,
                )
                result = solver.solve()
                min_configs.append(result["configuration"])

            # Detect configuration transitions
            transitions = []
            for i in range(len(min_configs) - 1):
                if min_configs[i] != min_configs[i + 1]:
                    transitions.append(
                        {
                            "occurs_between": (results["values"][i], results["values"][i + 1]),
                            "from_config": min_configs[i],
                            "to_config": min_configs[i + 1],
                        }
                    )

            if transitions:
                bifurcations[param_name] = transitions
                logger.info("Found %d bifurcations for parameter %s", len(transitions), param_name)

        return bifurcations

    def _set_parameter(self, param_name: str, value: float) -> None:
        """Set a parameter in particle or interaction.

        Args:
            param_name (str): Parameter name.
            value (float): New value to set.

        """
        if hasattr(self.particle, param_name):
            setattr(self.particle, param_name, value)
        elif hasattr(self.interaction, param_name):
            setattr(self.interaction, param_name, value)
        else:
            logger.warning("Parameter %s not found in particle or interaction", param_name)

    def _default_metric(
        self, protein: Protein, interaction: Interaction, particle: Particle
    ) -> float:
        """Default metric: minimum energy found by brute force.

        Args:
            protein: Protein object.
            interaction: Interaction model.
            particle: Particle.

        Returns:
            float: Minimum energy.

        """
        from validation.brute_force_folding import BruteForceFolding

        solver = BruteForceFolding(
            protein=protein, interaction=interaction, particle=particle, dimension=2
        )
        result = solver.solve()
        return result["energy"]

    def get_report(self) -> str:
        """Generate text report of sensitivity analysis.

        Returns:
            str: Formatted sensitivity analysis report.

        """
        report = "\n" + "=" * 60 + "\n"
        report += "SENSITIVITY ANALYSIS REPORT\n"
        report += "=" * 60 + "\n\n"

        report += f"Protein: {self.protein.main_chain.get_sequence()}\n"
        report += f"Particle type: {self.particle.particle_type.value}\n"
        report += f"Parameters scanned: {len(self.scan_results)}\n\n"

        for param_name, results in self.scan_results.items():
            if param_name.endswith("_2d_sweep"):
                report += f"2D Sweep: {param_name}\n"
                report += f"  Range 1: {min(results['hh_values']):.4f} to {max(results['hh_values']):.4f}\n"
                report += f"  Range 2: {min(results['pp_values']):.4f} to {max(results['pp_values']):.4f}\n"
            else:
                report += f"Parameter: {param_name}\n"
                metrics = [m for m in results["metrics"] if m is not None]
                if metrics:
                    report += f"  Range: {results['values'][0]:.4f} to {results['values'][-1]:.4f}\n"
                    report += f"  Metric range: {min(metrics):.6f} to {max(metrics):.6f}\n"

            report += "\n"

        report += "=" * 60 + "\n"
        return report
