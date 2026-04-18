"""Complete workflow example for particle extension research.

This module provides end-to-end example for:
1. Creating particle-extended systems
2. Validation against brute force reference
3. Baseline comparison (with vs without particle)
4. Sensitivity analysis
5. Genetic algorithm for sequence design
6. Result visualization and reporting
"""

from __future__ import annotations

from logger import get_logger

logger = get_logger()


def example_complete_workflow() -> None:
    """Run complete example workflow for particle studies.
    
    This demonstrates all major capabilities of the particle extension system.
    """
    
    logger.info("=" * 60)
    logger.info("PARTICLE EXTENSION COMPLETE WORKFLOW EXAMPLE")
    logger.info("=" * 60)
    
    # ===== Step 1: Setup =====
    logger.info("\n[1] SETUP: Creating particle and interaction models")
    
    from protein import Protein
    from constants import EMPTY_SIDECHAIN_PLACEHOLDER
    from interaction.hp_interaction_with_particle import HPInteractionWithParticle
    from particle.external_field_particle import ExternalFieldParticle
    
    # Create test protein (example sequence)
    main_chain_seq = "APRLRF"  # 6 beads - short for fast brute force
    side_chain_seq = EMPTY_SIDECHAIN_PLACEHOLDER * len(main_chain_seq)
    
    # Create interaction model with particle support
    hp_interaction = HPInteractionWithParticle(
        particle_symbol="L",
        particle_hh_energy=-1.0,  # Attractive to hydrophobic
        particle_pp_energy=0.0     # No interaction with polar
    )
    
    # Create protein
    protein = Protein(main_chain_seq, side_chain_seq, hp_interaction.valid_symbols)
    logger.info("[OK] Created protein sequence: %s", main_chain_seq)
    
    # Create particle (external field variant)
    particle = ExternalFieldParticle(
        position=(0, 0, 0),
        potential_type="distance",
        well_depth=-1.0,
        well_radius=2.0
    )
    logger.info("[OK] Created external field particle at origin")
    
    # ===== Step 2: Validation =====
    logger.info("\n[2] VALIDATION: Compare against brute force reference")
    
    from validation.validation_framework import ParticleValidationFramework
    
    validator = ParticleValidationFramework(protein, hp_interaction, particle)
    reference_solution = validator.validate_with_brute_force()
    
    logger.info("[OK] Reference solution: energy = %.4f", reference_solution["energy"])
    logger.info("  Configurations evaluated: %d", reference_solution["num_evaluated"])
    
    # ===== Step 3: Baseline Comparison =====
    logger.info("\n[3] BASELINE ANALYSIS: Compare with vs without particle")
    
    from analysis.baseline_comparison import BaselineComparisonWorkflow
    
    comparison_workflow = BaselineComparisonWorkflow(
        protein, hp_interaction, particle
    )
    comparison = comparison_workflow.run_baseline_comparison()
    
    logger.info("[OK] Without particle: energy = %.4f", comparison.energy_no_particle)
    logger.info("[OK] With particle:    energy = %.4f", comparison.energy_with_particle)
    logger.info("[OK] Energy delta: %.4f", comparison.energy_delta)
    logger.info("[OK] Stabilizing: %s", comparison.is_stabilizing)
    
    # ===== Step 4: Energy Landscape Analysis =====
    logger.info("\n[4] LANDSCAPE ANALYSIS: Study energy landscape properties")
    
    from analysis.energy_landscape import EnergyLandscapeAnalyzer
    
    analyzer = EnergyLandscapeAnalyzer(protein, hp_interaction, particle)
    
    # Analyze topology changes
    topology_analysis = analyzer.analyze_topology_changes()
    logger.info("[OK] Topology preserved: %s", topology_analysis["topology_preserved"])
    logger.info("[OK] Baseline compactness: %.4f", topology_analysis["baseline_compactness"])
    logger.info("[OK] With particle compactness: %.4f", topology_analysis["particle_compactness"])
    
    # ===== Step 5: Sensitivity Analysis =====
    logger.info("\n[5] SENSITIVITY: Parameter sweep analysis")
    
    from analysis.sensitivity_analysis import SensitivityAnalysisFramework
    
    sensitivity = SensitivityAnalysisFramework(protein, hp_interaction, particle)
    
    # Sweep particle interaction strength with hydrophobic residues
    hh_values = [-2.0, -1.5, -1.0, -0.5, 0.0]
    sweep_results = sensitivity.sweep_single_parameter(
        "particle_hh_energy",
        hh_values
    )
    
    logger.info("[OK] Swept particle-hydrophobic interaction from %.1f to %.1f",
                hh_values[0], hh_values[-1])
    logger.info("  Results: %s", [f"{m:.3f}" if m else "None" for m in sweep_results["metrics"]])
    
    # Find critical parameters
    critical_params = sensitivity.identify_critical_parameters(threshold=0.1)
    if critical_params:
        logger.info("[OK] Critical parameters found:")
        for param, info in critical_params.items():
            logger.info("  - %s: sensitivity = %.4f", param, info["sensitivity"])
    else:
        logger.info("[OK] No critical parameters identified (< threshold)")
    
    # ===== Step 6: Sequence Design =====
    logger.info("\n[6] SEQUENCE DESIGN: Genetic algorithm optimization")
    
    from optimization.sequence_design_ga import (
        SequenceDesignGA,
        objective_stabilization_effect
    )
    
    ga = SequenceDesignGA(
        interaction=hp_interaction,
        particle=particle,
        objective=objective_stabilization_effect,
        sequence_length=6,
        population_size=20,
        generations=30  # Small for this example
    )
    
    best_individuals = ga.run_evolution(verbose=False)
    best_sequence = ga.get_best_sequence()
    
    logger.info("[OK] GA evolution complete")
    logger.info("[OK] Best sequence designed: %s", best_sequence)
    logger.info("[OK] Fitness improvement from gen 0 to %d: %.3f%%",
                len(best_individuals)-1,
                ((best_individuals[-1].fitness - best_individuals[0].fitness) 
                 / abs(best_individuals[0].fitness)) * 100 if best_individuals[0].fitness != 0 else 0)
    
    # ===== Step 7: Visualization and Reports =====
    logger.info("\n[7] VISUALIZATION: Generate publication-ready figures")
    
    from visualization.visualization_tools import VisualizationTools
    
    # Generate convergence plot
    VisualizationTools.plot_ga_convergence(ga, output_path="output/results/ga_convergence.png")
    logger.info("[OK] Saved GA convergence plot")
    
    # Generate baseline comparison bar plot
    from analysis.baseline_comparison import BaselineComparison
    baseline_data = [comparison]
    labels = ["Baseline comparison"]
    VisualizationTools.plot_baseline_comparison_bar(
        baseline_data, labels, output_path="output/results/baseline_comparison.png"
    )
    logger.info("[OK] Saved baseline comparison plot")
    
    # ===== Step 8: Generate Reports =====
    logger.info("\n[8] REPORTS: Comprehensive analysis reports")
    
    report = comparison_workflow.get_report()
    logger.info(report)
    
    sensitivity_report = sensitivity.get_report()
    logger.info(sensitivity_report)
    
    # ===== Summary =====
    logger.info("\n" + "=" * 60)
    logger.info("WORKFLOW COMPLETE - SUMMARY")
    logger.info("=" * 60)
    logger.info("[OK] Validation: Reference solution established")
    logger.info("[OK] Baseline: Energy stabilization quantified")
    logger.info("[OK] Landscape: Topology and compactness analyzed")
    logger.info("[OK] Sensitivity: Critical parameters identified")
    logger.info("[OK] Design: Optimized sequences discovered via GA")
    logger.info("[OK] Visualization: Publication-quality figures generated")
    logger.info("=" * 60 + "\n")


def example_dynamic_bead_particle() -> None:
    """Example using dynamic bead variant instead of external field.
    
    Demonstrates the more complex particle-as-bead approach.
    """
    
    logger.info("=" * 60)
    logger.info("DYNAMIC BEAD PARTICLE EXAMPLE")
    logger.info("=" * 60)
    
    from protein import Protein
    from constants import EMPTY_SIDECHAIN_PLACEHOLDER
    from interaction.mj_interaction_with_particle import MJInteractionWithParticle
    from particle.dynamic_bead_particle import DynamicBeadParticle
    from analysis.baseline_comparison import BaselineComparisonWorkflow
    
    # Create MJ interaction model with particle
    mj_interaction = MJInteractionWithParticle(particle_symbol="L")
    
    # Define particle-residue interactions
    # (typically would load from matrix file)
    particle_interactions = {
        "A": -0.5, "R": 0.2, "N": 0.1, "D": 0.3, "C": -0.8,
        "E": 0.4, "Q": 0.1, "G": 0.0, "H": -0.2, "I": -0.9,
        "K": 0.3, "L": -0.8, "M": -0.7, "F": -0.9, "P": 0.1,
        "S": 0.0, "T": 0.0, "W": -0.8, "Y": -0.6, "V": -0.8,
    }
    mj_interaction.set_particle_interactions(particle_interactions)
    
    # Create protein
    main_chain = "APRLRF"
    protein = Protein(main_chain, "------", mj_interaction.valid_symbols)
    
    # Create dynamic bead particle
    particle = DynamicBeadParticle(
        position=(0, 0, 0),
        symbol="L",
        fixed=False  # Can move
    )
    particle.set_interaction_potentials(particle_interactions)
    
    logger.info("[OK] Created dynamic bead particle")
    logger.info("  Position: %s", particle.position)
    logger.info("  Fixed: %s", particle.fixed)
    
    # Run baseline comparison
    workflow = BaselineComparisonWorkflow(protein, mj_interaction, particle)
    comparison = workflow.run_baseline_comparison()
    
    logger.info("[OK] Baseline analysis complete")
    logger.info("  Energy delta: %.4f", comparison.energy_delta)
    logger.info("  Stabilization factor: %.4f", comparison.stabilization_factor)


if __name__ == "__main__":
    # Run the complete workflow example
    example_complete_workflow()
    
    # Optionally, run dynamic bead example
    # example_dynamic_bead_particle()
