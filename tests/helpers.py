"""Shared helpers for the test-suite."""

from __future__ import annotations

from functools import cache

import numpy as np

from analysis.exact_solver import to_diagonal
from builder import HamiltonianBuilder
from constants import EMPTY_SIDECHAIN_PLACEHOLDER
from utils.setup_utils import setup_folding_system


def diagonal_of(operator):
    """Return the diagonal of a Z-type operator as a real array."""
    return to_diagonal(operator)


@cache
def build_system(sequence: str):
    """Set up the folding system for a sequence, with empty side chains.

    Cached because the contact and distance maps are expensive to rebuild and
    are treated as read-only by every caller.
    """
    return setup_folding_system(
        main_chain=sequence,
        side_chain=EMPTY_SIDECHAIN_PLACEHOLDER * len(sequence),
    )


def build_hamiltonian(sequence: str, **kwargs):
    """Build a Hamiltonian for a sequence, forwarding optional terms."""
    protein, interaction, contact_map, distance_map = build_system(sequence)
    builder = HamiltonianBuilder(
        protein=protein,
        interaction=interaction,
        distance_map=distance_map,
        contact_map=contact_map,
        **kwargs,
    )
    return builder, builder.sum_hamiltonians()


def assert_operators_equal(left, right):
    """Assert two diagonal operators agree on every basis state."""
    assert left.num_qubits == right.num_qubits
    assert np.allclose(to_diagonal(left), to_diagonal(right))
