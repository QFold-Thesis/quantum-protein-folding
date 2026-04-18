"""Validation framework for particle folding implementations.

Compares results between different solvers (brute force, quantum, classical)
to validate correctness of particle extensions.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from logger import get_logger

if TYPE_CHECKING:
    from interaction.interaction import Interaction
    from particle import Particle
    from protein import Protein

logger = get_logger()


class ParticleValidationFramework:
    """Framework for validating particle folding implementations.

    Attributes:
        protein (Protein): Test protein.
        interaction (Interaction): Interaction model.
        particle (Particle, optional): Particle for testing.

    """

    def __init__(
        self,
        protein: Protein,
        interaction: Interaction,
        particle: Particle | None = None,
    ) -> None:
        """Initialize validation framework.

       Args:
            protein (Protein): Protein to validate.
            interaction (Interaction): Interaction model.
            particle (Particle, optional): Particle for testing.

        """
        self.protein = protein
        self.interaction = interaction
        self.particle = particle
        self.results = {}

        logger.info(
            "Initialized validation framework for %d-bead protein %s",
            len(protein.main_chain),
            "with particle" if particle else "without particle",
        )

    def validate_with_brute_force(self) -> dict:
        """Run brute force solver and store reference results.

        Returns:
            dict: Reference solution from exhaustive search.

        """
        from validation.brute_force_folding import BruteForceFolding

        logger.info("Starting brute force validation...")

        solver = BruteForceFolding(
            protein=self.protein,
            interaction=self.interaction,
            particle=self.particle,
            dimension=2,
        )

        reference_solution = solver.solve()
        self.results["brute_force"] = reference_solution

        logger.info(
            "Brute force reference: energy = %f, checked %d configurations",
            reference_solution["energy"],
            reference_solution["num_evaluated"],
        )

        return reference_solution

    def compare_energies(
        self, test_solution: dict, tolerance: float = 1e-6
    ) -> dict:
        """Compare test solution energy with brute force reference.

        Args:
            test_solution (dict): Solution from another solver.
            tolerance (float, optional): Energy difference tolerance. Defaults to 1e-6.

        Returns:
            dict: Comparison results including energy difference and pass/fail.

        Raises:
            ValueError: If no reference solution available.

        """
        if "brute_force" not in self.results:
            msg = "Must run validate_with_brute_force() first"
            raise ValueError(msg)

        reference_energy = self.results["brute_force"]["energy"]
        test_energy = test_solution["energy"]
        energy_diff = abs(reference_energy - test_energy)
        relative_diff = energy_diff / (abs(reference_energy) + 1e-10)

        passed = energy_diff <= tolerance

        logger.info(
            "Energy comparison: reference=%f, test=%f, diff=%f (%f%%)",
            reference_energy,
            test_energy,
            energy_diff,
            relative_diff * 100,
        )

        return {
            "reference_energy": reference_energy,
            "test_energy": test_energy,
            "absolute_difference": energy_diff,
            "relative_difference": relative_diff,
            "tolerance": tolerance,
            "passed": passed,
        }

    def validate_particle_effect(self) -> dict:
        """Analyze effect of particle on system energy.

        Compares energy with and without particle to quantify particle impact.

        Returns:
            dict: Analysis of particle's effect on energy landscape.

        """
        logger.info("Analyzing particle effect...")

        # Run without particle (baseline)
        from validation.brute_force_folding import BruteForceFolding

        solver_no_particle = BruteForceFolding(
            protein=self.protein,
            interaction=self.interaction,
            particle=None,
            dimension=2,
        )
        baseline_solution = solver_no_particle.solve()

        # Run with particle (if defined)
        if self.particle:
            solver_with_particle = BruteForceFolding(
                protein=self.protein,
                interaction=self.interaction,
                particle=self.particle,
                dimension=2,
            )
            particle_solution = solver_with_particle.solve()
        else:
            particle_solution = baseline_solution

        baseline_energy = baseline_solution["energy"]
        particle_energy = particle_solution["energy"]
        energy_difference = particle_energy - baseline_energy
        stabilization_factor = baseline_energy / (particle_energy + 1e-10)

        logger.info(
            "Particle effect: baseline=%f, with_particle=%f, delta=%f",
            baseline_energy,
            particle_energy,
            energy_difference,
        )

        return {
            "baseline_energy": baseline_energy,
            "energy_with_particle": particle_energy,
            "energy_difference": energy_difference,
            "is_stabilizing": energy_difference < 0,
            "stabilization_factor": stabilization_factor,
            "baseline_config": baseline_solution["configuration"],
            "particle_config": particle_solution["configuration"],
        }

    def validate_interaction_symmetry(self) -> dict:
        """Check if particle-residue interactions are symmetric where expected.

        Returns:
            dict: Symmetry validation results.

        """
        logger.info("Validating interaction symmetry...")

        if not self.particle:
            return {"skipped": True, "reason": "No particle defined"}

        errors = []
        main_chain = self.protein.main_chain
        symbols = set(bead.symbol for bead in main_chain.beads)

        for sym in symbols:
            # Check symmetry: interaction(A, L) should equal interaction(L, A)
            try:
                e1 = self.interaction.get_energy(sym, self.particle.symbol)
                e2 = self.interaction.get_energy(self.particle.symbol, sym)

                if abs(e1 - e2) > 1e-10:
                    errors.append(f"Asymmetry for {sym}-{self.particle.symbol}: {e1} vs {e2}")
            except Exception as e:
                errors.append(f"Error checking {sym}: {e}")

        logger.info("Symmetry check found %d errors", len(errors))

        return {
            "passed": len(errors) == 0,
            "errors": errors,
            "symbols_checked": len(symbols),
        }

    def get_validation_report(self) -> str:
        """Generate text report of validation results.

        Returns:
            str: Human-readable validation report.

        """
        report = "\n" + "=" * 60 + "\n"
        report += "PARTICLE FOLDING VALIDATION REPORT\n"
        report += "=" * 60 + "\n\n"

        if "brute_force" in self.results:
            bf = self.results["brute_force"]
            report += f"Brute Force Solution:\n"
            report += f"  Energy: {bf['energy']:.6f}\n"
            report += f"  Configurations evaluated: {bf['num_evaluated']}\n\n"

        report += "=" * 60 + "\n"
        return report
