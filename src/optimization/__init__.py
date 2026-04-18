"""Optimization module for sequence design.

Provides genetic algorithms and other optimization methods for designing
sequences with target behaviors relative to particles.
"""

from optimization.sequence_design_ga import (
    Individual,
    SequenceDesignGA,
    objective_maximize_encirclement,
    objective_minimize_energy,
    objective_stabilization_effect,
)

__all__ = [
    "Individual",
    "SequenceDesignGA",
    "objective_maximize_encirclement",
    "objective_minimize_energy",
    "objective_stabilization_effect",
]
