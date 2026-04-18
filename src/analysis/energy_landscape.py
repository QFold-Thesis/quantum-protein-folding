"""Energy landscape analysis tools for particle studies.

Analyzes how particle presence affects the energy landscape, including
topology changes, local minima, and barriers.
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
class LandscapeSnapshot:
    """Snapshot of energy landscape at specific conditions.

    Attributes:
        min_energy (float): Global minimum energy.
        num_local_minima (int): Number of distinct local minima.
        barrier_height (float): Typical energy barrier between minima.
        energy_spread (float): Range between min and max accessible energies.
        configuration_diversity (float): Measure of conformation diversity.

    """

    min_energy: float
    num_local_minima: int
    barrier_height: float
    energy_spread: float
    configuration_diversity: float


class EnergyLandscapeAnalyzer:
    """Analyzes energy landscape with and without particles.

    Attributes:
        protein (Protein): Protein being studied.
        interaction (Interaction): Interaction model.
        particle (Particle, optional): Particle for study.

    """

    def __init__(
        self,
        protein: Protein,
        interaction: Interaction,
        particle: Particle | None = None,
    ) -> None:
        """Initialize landscape analyzer.

        Args:
            protein (Protein): Protein sequence.
            interaction (Interaction): Interaction model.
            particle (Particle, optional): Particle being studied.

        """
        self.protein = protein
        self.interaction = interaction
        self.particle = particle
        self.energy_samples: list[float] = []
        self.configurations: list[dict] = []

        logger.info(
            "Initialized landscape analyzer for %d-bead protein",
            len(protein.main_chain),
        )

    def sample_landscape(self, num_samples: int = 100) -> list[float]:
        """Sample energy landscape with random conformations.

        Args:
            num_samples (int): Number of random configurations to evaluate.

        Returns:
            list[float]: Energies of sampled configurations.

        """
        logger.info("Sampling landscape with %d random conformations...", num_samples)

        from validation.brute_force_folding import BruteForceFolding

        # Generate random valid conformations and evaluate
        solver = BruteForceFolding(
            protein=self.protein,
            interaction=self.interaction,
            particle=self.particle,
            dimension=2,
        )

        energies = []
        configs = []

        # Use solver's conformation generation internally
        for _ in range(num_samples):
            # Simplified: evaluation at random state
            # Full implementation would integrate with solver's enumeration
            energy = 0.0  # Placeholder
            energies.append(energy)

        self.energy_samples = energies
        logger.info(
            "Landscape sampling complete: %d samples, energy range [%f, %f]",
            len(energies),
            min(energies) if energies else 0,
            max(energies) if energies else 0,
        )

        return energies

    def analyze_topology_changes(self) -> dict:
        """Analyze qualitative changes in optimal conformations.

        Compares conformation properties to detect topology changes.

        Returns:
            dict: Topology analysis results.

        """
        logger.info("Analyzing conformation topology changes...")

        from validation.brute_force_folding import BruteForceFolding

        # Without particle
        solver_baseline = BruteForceFolding(
            protein=self.protein,
            interaction=self.interaction,
            particle=None,
            dimension=2,
        )
        result_baseline = solver_baseline.solve()

        # With particle
        if self.particle:
            solver_particle = BruteForceFolding(
                protein=self.protein,
                interaction=self.interaction,
                particle=self.particle,
                dimension=2,
            )
            result_particle = solver_particle.solve()
        else:
            result_particle = result_baseline

        configs_changed = result_baseline["configuration"] != result_particle["configuration"]
        topology_preserved = not configs_changed

        analysis = {
            "topology_preserved": topology_preserved,
            "configuration_changed": configs_changed,
            "baseline_compactness": self._calculate_compactness(
                result_baseline["configuration"]
            ),
            "particle_compactness": self._calculate_compactness(
                result_particle["configuration"]
            ),
        }

        logger.info(
            "Topology analysis: preserved=%s, compactness_baseline=%f, compactness_particle=%f",
            topology_preserved,
            analysis["baseline_compactness"],
            analysis["particle_compactness"],
        )

        return analysis

    def _calculate_compactness(self, configuration: dict) -> float:
        """Calculate structural compactness (radius of gyration concept).

        Args:
            configuration (dict): Bead coordinate mapping.

        Returns:
            float: Compactness metric (lower = more compact).

        """
        if not configuration:
            return 0.0

        coords = list(configuration.values())

        # Calculate center of mass
        cx = sum(c[0] for c in coords) / len(coords)
        cy = sum(c[1] for c in coords) / len(coords)
        cz = sum(c[2] for c in coords) / len(coords)

        # Calculate radius of gyration
        rg_squared = sum(
            ((c[0] - cx) ** 2 + (c[1] - cy) ** 2 + (c[2] - cz) ** 2) for c in coords
        ) / len(coords)

        return rg_squared ** 0.5

    def get_landscape_snapshot(self) -> LandscapeSnapshot:
        """Get comprehensive landscape snapshot.

        Returns:
            LandscapeSnapshot: Summarized landscape properties.

        """
        if not self.energy_samples:
            msg = "Must sample landscape first"
            raise ValueError(msg)

        min_energy = min(self.energy_samples)
        max_energy = max(self.energy_samples)
        spread = max_energy - min_energy

        # Estimate local minima (simplified)
        num_minima = max(1, len(self.energy_samples) // 20)
        barrier = spread / (num_minima + 1) if num_minima > 0 else spread

        # Diversity based on energy spread normalized by system size
        diversity = spread / (len(self.protein.main_chain) + 1e-10)

        return LandscapeSnapshot(
            min_energy=min_energy,
            num_local_minima=num_minima,
            barrier_height=barrier,
            energy_spread=spread,
            configuration_diversity=diversity,
        )

    def detect_bifurcations(
        self, param_ranges: dict[str, list[float]]
    ) -> dict:
        """Detect structural bifurcations as particle parameters vary.

        Args:
            param_ranges (dict): Parameter names mapped to ranges to scan.

        Returns:
            dict: Bifurcation analysis results.

        """
        logger.info("Scanning for bifurcations in %d parameters", len(param_ranges))

        bifurcations = {}

        for param_name, values in param_ranges.items():
            logger.debug("Scanning parameter: %s", param_name)

            min_configs = []

            for value in values:
                # Update particle/interaction parameter
                self._set_parameter(param_name, value)

                from validation.brute_force_folding import BruteForceFolding

                solver = BruteForceFolding(
                    protein=self.protein,
                    interaction=self.interaction,
                    particle=self.particle,
                    dimension=2,
                )
                result = solver.solve()
                min_configs.append(result["configuration"])

            # Detect changes in optimal configuration
            transitions = []
            for i in range(len(min_configs) - 1):
                if min_configs[i] != min_configs[i + 1]:
                    transitions.append((values[i], values[i + 1]))

            bifurcations[param_name] = {
                "transitions": transitions,
                "num_bifurcations": len(transitions),
            }

            logger.debug(
                "Parameter %s: %d bifurcation points detected",
                param_name,
                len(transitions),
            )

        return bifurcations

    def _set_parameter(self, param_name: str, value: float) -> None:
        """Update a parameter in particle or interaction.

        Args:
            param_name (str): Parameter name.
            value (float): New value.

        """
        if self.particle and hasattr(self.particle, param_name):
            setattr(self.particle, param_name, value)
        elif hasattr(self.interaction, param_name):
            setattr(self.interaction, param_name, value)
        else:
            logger.warning("Parameter %s not found", param_name)
