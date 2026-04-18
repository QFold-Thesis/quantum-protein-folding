"""Particle module for extending protein folding with external particle interactions.

This module provides abstractions for modeling external particles that interact
with protein chains in lattice/diamond representations, supporting two variants:
1. External field: Static interaction potential
2. Dynamic bead: Particle with its own lattice coordinates and mobility
"""

from particle.particle import Particle
from particle.external_field_particle import ExternalFieldParticle
from particle.dynamic_bead_particle import DynamicBeadParticle

__all__ = [
    "Particle",
    "ExternalFieldParticle",
    "DynamicBeadParticle",
]
