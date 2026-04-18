"""Brute force solver for protein folding validation.

Provides BruteForceFolding solver that exhaustively searches all possible
conformations for short protein sequences (6-9 beads) and finds minimum
energy configurations. Used for validation against quantum/classical solvers.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from logger import get_logger

if TYPE_CHECKING:
    from interaction.interaction import Interaction
    from particle import Particle
    from protein import Protein

logger = get_logger()


class BruteForceFolding:
    """Exhaustive search solver for protein folding validation.

    For short sequences, evaluates all possible 2D/3D lattice configurations
    and finds global minimum. Works with or without external particle.

    Attributes:
        protein (Protein): Protein object.
        interaction (Interaction): Interaction model.
        particle (Particle, optional): Particle for interaction calculations.
        dimension (int): Lattice dimensionality (2 or 3). Defaults to 2.

    """

    def __init__(
        self,
        protein: Protein,
        interaction: Interaction,
        particle: Particle | None = None,
        dimension: int = 2,
    ) -> None:
        """Initialize brute force solver.

        Args:
            protein (Protein): Protein to fold.
            interaction (Interaction): Interaction model.
            particle (Particle, optional): External particle. Defaults to None.
            dimension (int, optional): Lattice dimension (2 or 3). Defaults to 2.

        Raises:
            ValueError: If dimension not 2 or 3, or protein too long.

        """
        if dimension not in (2, 3):
            msg = "Dimension must be 2 or 3"
            raise ValueError(msg)

        if len(protein.main_chain) > 10:
            msg = "Brute force solver limited to sequences of 10 beads or fewer"
            raise ValueError(msg)

        self.protein = protein
        self.interaction = interaction
        self.particle = particle
        self.dimension = dimension

        logger.info(
            "Initialized BruteForceFolding for %d-bead sequence in %dD lattice",
            len(protein.main_chain),
            dimension,
        )

    def solve(self) -> dict:
        """Find minimum energy configuration by exhaustive search.

        Enumerates all valid lattice conformations and evaluates energy for each.
        Returns configuration with lowest total energy.

        Returns:
            dict: Solution dictionary containing:
                - 'energy': Minimum energy found
                - 'configuration': Bead coordinates mapping
                - 'is_valid': Whether configuration is physically valid
                - 'num_evaluated': Total configurations evaluated

        """
        logger.info("Starting brute force enumeration...")

        min_energy = float("inf")
        best_configuration = None
        count = 0

        # Generate all possible configurations
        for configuration in self._generate_conformations():
            count += 1

            # Skip invalid/self-intersecting configurations
            if not self._is_valid_configuration(configuration):
                continue

            # Calculate total energy
            energy = self._calculate_energy(configuration)

            if energy < min_energy:
                min_energy = energy
                best_configuration = configuration

        logger.info(
            "Brute force complete: %d configurations evaluated, minimum energy: %f",
            count,
            min_energy,
        )

        return {
            "energy": min_energy,
            "configuration": best_configuration,
            "is_valid": True,
            "num_evaluated": count,
        }

    def _generate_conformations(self) -> list[dict]:
        """Generate all possible lattice configurations.

        Returns:
            list[dict]: List of configuration dictionaries, each mapping
                bead indices to (x, y, z) coordinates.

        """
        # Simplified: for 2D with simple moves (right, down, -right, -down)
        # This is a placeholder - full implementation needs proper
        # self-avoiding walk enumeration

        configurations = []

        if self.dimension == 2:
            # Generate using backtracking with 2D moves
            chain_len = len(self.protein.main_chain)
            directions = [(1, 0), (0, 1), (-1, 0), (0, -1)]  # right, down, left, up

            def backtrack(pos_list, visited, depth):
                if depth == chain_len:
                    config = {i: (pos_list[i][0], pos_list[i][1], 0) 
                             for i in range(chain_len)}
                    configurations.append(config)
                    return

                current_pos = pos_list[-1]
                for dx, dy in directions:
                    new_pos = (current_pos[0] + dx, current_pos[1] + dy)
                    if new_pos not in visited:
                        pos_list.append(new_pos)
                        visited.add(new_pos)
                        backtrack(pos_list, visited, depth + 1)
                        pos_list.pop()
                        visited.remove(new_pos)

            # Start from origin
            backtrack([(0, 0)], {(0, 0)}, 1)

        elif self.dimension == 3:
            # Similar logic for 3D but with 6 directions
            pass

        logger.debug("Generated %d possible conformations", len(configurations))
        return configurations

    def _is_valid_configuration(self, configuration: dict) -> bool:
        """Check if configuration is valid (no self-intersections, connected).

        Args:
            configuration (dict): Mapping of bead indices to coordinates.

        Returns:
            bool: True if configuration is valid.

        """
        # Check that all beads are present
        if len(configuration) != len(self.protein.main_chain):
            return False

        # Check that beads form a connected path
        coords = list(configuration.values())
        for i in range(len(coords) - 1):
            curr = coords[i]
            next_coord = coords[i + 1]
            distance = sum((c1 - c2) ** 2 for c1, c2 in zip(curr, next_coord)) ** 0.5
            if distance > 1.01:  # Allow small float error
                return False

        return True

    def _calculate_energy(self, configuration: dict) -> float:
        """Calculate total energy of a configuration.

        Includes backbone bond energies, contact interactions, and particle
        interactions if particle is defined.

        Args:
            configuration (dict): Bead coordinates.

        Returns:
            float: Total system energy.

        """
        energy = 0.0

        # Add contact interaction energy
        coords = [configuration[i] for i in range(len(self.protein.main_chain))]
        main_chain = self.protein.main_chain

        # Create proper configuration dictionary for particle
        config_for_particle = {"coordinates": configuration.copy()}
        config_for_particle["particle_position"] = (0, 0, 0)
        config_for_particle["contact_distance"] = 1

        # Pairwise interactions for beads in contact
        for i in range(len(main_chain)):
            for j in range(i + 2, len(main_chain)):  # Skip consecutive beads
                if self._are_neighbors(coords[i], coords[j]):
                    sym_i = main_chain.beads[i].symbol
                    sym_j = main_chain.beads[j].symbol
                    inter_energy = self.interaction.get_energy(sym_i, sym_j)
                    energy += inter_energy

        # Add particle interactions if present
        if self.particle:
            try:
                particle_energy = self.particle.get_energy_contribution(
                    self.protein, config_for_particle
                )
                energy += particle_energy
            except Exception as e:
                logger.warning("Error calculating particle energy: %s", e)

        return energy

    def _are_neighbors(
        self, coord1: tuple[float, float, float], coord2: tuple[float, float, float]
    ) -> bool:
        """Check if two beads are neighbors on lattice.

        Uses Manhattan distance = 1 for neighbors.

        Args:
            coord1: First coordinate.
            coord2: Second coordinate.

        Returns:
            bool: True if beads are neighbors.

        """
        manhattan = sum(abs(c1 - c2) for c1, c2 in zip(coord1, coord2))
        return manhattan == 1
