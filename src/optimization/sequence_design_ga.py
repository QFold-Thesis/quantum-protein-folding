"""Genetic algorithm for inverse sequence design problem.

Finds protein sequences that exhibit target behaviors relative to a particle.
Behaviors include: maximizing particle encirclement, generating repulsive forces,
or achieving specific contact patterns.
"""

from __future__ import annotations

import random
from typing import TYPE_CHECKING, Callable

from logger import get_logger

if TYPE_CHECKING:
    from interaction.interaction import Interaction
    from particle import Particle
    from protein import Protein

logger = get_logger()


class Individual:
    """Represents a candidate sequence in genetic algorithm.

    Attributes:
        sequence (str): Amino acid sequence.
        fitness (float): Fitness metric for this sequence (higher is better).
        energy (float, optional): System energy for this sequence.

    """

    def __init__(self, sequence: str, valid_symbols: set[str]):
        """Initialize individual with sequence.

        Args:
            sequence (str): Amino acid sequence.
            valid_symbols (set[str]): Valid symbols for this evolution.

        """
        self.sequence = sequence
        self.valid_symbols = valid_symbols
        self.fitness: float | None = None
        self.energy: float | None = None

    def mutate(self, mutation_rate: float = 0.1) -> Individual:
        """Create mutated copy of this individual.

        Args:
            mutation_rate (float): Probability per position to mutate.

        Returns:
            Individual: New mutated individual.

        """
        new_sequence = list(self.sequence)

        for i in range(len(new_sequence)):
            if random.random() < mutation_rate:
                new_sequence[i] = random.choice(list(self.valid_symbols))

        return Individual("".join(new_sequence), self.valid_symbols)

    def crossover(self, other: Individual) -> Individual:
        """Create offspring via crossover with another individual.

        Args:
            other (Individual): Other parent.

        Returns:
            Individual: Offspring individual.

        """
        if len(self.sequence) != len(other.sequence):
            msg = "Parents must have equal sequence length"
            raise ValueError(msg)

        # Single-point crossover
        crossover_point = random.randint(1, len(self.sequence) - 1)
        new_sequence = self.sequence[:crossover_point] + other.sequence[crossover_point:]

        return Individual(new_sequence, self.valid_symbols)

    def __repr__(self) -> str:
        """Return string representation."""
        return f"Individual(seq={self.sequence}, fitness={self.fitness:.4f})"


class SequenceDesignGA:
    """Genetic algorithm for sequence design.

    Evolves protein sequences to optimize target behaviors relative to particle.

    Attributes:
        interaction (Interaction): Interaction model.
        particle (Particle): Target particle.
        objective (Callable): Fitness function to optimize.
        population_size (int): Number of individuals per generation.
        generations (int): Number of generations to evolve.

    """

    def __init__(
        self,
        interaction: Interaction,
        particle: Particle,
        objective: Callable,
        sequence_length: int = 8,
        population_size: int = 50,
        generations: int = 100,
    ) -> None:
        """Initialize genetic algorithm.

        Args:
            interaction (Interaction): Interaction model to use.
            particle (Particle): Particle for sequence design.
            objective (Callable): Fitness function (sequence, interaction, particle) -> float.
            sequence_length (int, optional): Length of sequences to design. Defaults to 8.
            population_size (int, optional): Population size. Defaults to 50.
            generations (int, optional): Number of generations. Defaults to 100.

        """
        self.interaction = interaction
        self.particle = particle
        self.objective = objective
        self.sequence_length = sequence_length
        self.population_size = population_size
        self.generations = generations
        self.valid_symbols = interaction.valid_symbols.copy()
        self.valid_symbols.discard(particle.symbol)  # Don't use particle symbol

        self.population: list[Individual] = []
        self.best_individuals: list[Individual] = []
        self.generation: int = 0

        logger.info(
            "Initialized SequenceDesignGA: length=%d, pop=%d, gen=%d, symbols=%d",
            sequence_length,
            population_size,
            generations,
            len(self.valid_symbols),
        )

    def initialize_population(self, include_random: bool = True) -> None:
        """Create initial population.

        Args:
            include_random (bool, optional): Start with random sequences. Defaults to True.

        """
        logger.info("Initializing population of %d individuals", self.population_size)

        self.population = []

        if include_random:
            for _ in range(self.population_size):
                seq = "".join(random.choice(list(self.valid_symbols)) for _ in range(self.sequence_length))
                self.population.append(Individual(seq, self.valid_symbols))
        else:
            # Could initialize with known good sequences
            pass

        logger.debug("Initial population created")

    def evaluate_fitness(self, individual: Individual) -> float:
        """Evaluate fitness of individual sequence.

        Args:
            individual (Individual): Individual to evaluate.

        Returns:
            float: Fitness value (higher is better).

        """
        # Create temporary protein for evaluation
        from protein import Protein
        from constants import EMPTY_SIDECHAIN_PLACEHOLDER

        try:
            side_chain = EMPTY_SIDECHAIN_PLACEHOLDER * len(individual.sequence)
            protein = Protein(individual.sequence, side_chain, self.interaction.valid_symbols)

            # Evaluate fitness using objective function
            fitness = self.objective(protein, self.interaction, self.particle)
            individual.fitness = fitness
            return fitness

        except Exception as e:
            logger.warning("Error evaluating sequence %s: %s", individual.sequence, e)
            individual.fitness = 0.0
            return 0.0

    def run_evolution(self, verbose: bool = True) -> list[Individual]:
        """Run genetic algorithm evolution.

        Args:
            verbose (bool, optional): Print progress. Defaults to True.

        Returns:
            list[Individual]: Best individuals from each generation.

        """
        logger.info("Starting evolution for %d generations", self.generations)

        self.initialize_population()
        self.best_individuals = []

        for gen in range(self.generations):
            self.generation = gen

            # Evaluate fitness
            for individual in self.population:
                if individual.fitness is None:
                    self.evaluate_fitness(individual)

            # Sort by fitness (descending)
            self.population.sort(key=lambda x: x.fitness, reverse=True)
            best = self.population[0]
            self.best_individuals.append(best)

            if verbose and gen % 10 == 0:
                logger.info(
                    "Gen %d: best_fitness=%.4f, seq=%s",
                    gen,
                    best.fitness,
                    best.sequence,
                )

            # Selection and reproduction
            new_population = []

            # Elitism: keep best individuals
            elite_size = max(2, self.population_size // 10)
            new_population.extend(self.population[:elite_size])

            # Generate offspring
            while len(new_population) < self.population_size:
                # Tournament selection
                parent1 = self._tournament_selection(tournament_size=3)
                parent2 = self._tournament_selection(tournament_size=3)

                # Crossover and mutation
                offspring = parent1.crossover(parent2)
                offspring = offspring.mutate(mutation_rate=0.15)

                new_population.append(offspring)

            self.population = new_population

        # Ensure final population is evaluated
        for individual in self.population:
            if individual.fitness is None:
                self.evaluate_fitness(individual)

        logger.info("Evolution complete")
        return self.best_individuals

    def _tournament_selection(self, tournament_size: int = 3) -> Individual:
        """Select individual via tournament selection.

        Args:
            tournament_size (int): Size of tournament.

        Returns:
            Individual: Selected individual.

        """
        tournament = [random.choice(self.population) for _ in range(tournament_size)]
        # Evaluate fitness if needed
        for ind in tournament:
            if ind.fitness is None:
                self.evaluate_fitness(ind)
        return max(tournament, key=lambda x: x.fitness)

    def get_best_sequence(self) -> str:
        """Get best sequence found.

        Returns:
            str: Best amino acid sequence.

        """
        if not self.best_individuals:
            msg = "Must run evolution first"
            raise ValueError(msg)

        # Ensure all individuals are evaluated
        for individual in self.population:
            if individual.fitness is None:
                self.evaluate_fitness(individual)

        # Return best from final generation
        best = max(self.population, key=lambda x: x.fitness)
        return best.sequence

    def get_convergence_data(self) -> tuple[list[float], list[float]]:
        """Get convergence history.

        Returns:
            tuple: (generations, best_fitness_values)

        """
        fitnesses = [ind.fitness for ind in self.best_individuals]
        generations = list(range(len(fitnesses)))
        return generations, fitnesses


def objective_maximize_encirclement(
    protein: Protein, interaction: Interaction, particle: Particle
) -> float:
    """Objective: maximize number of particle-chain contacts.

    Higher contacts = better encirclement.

    Args:
        protein (Protein): Candidate protein.
        interaction (Interaction): Interaction model.
        particle (Particle): Target particle.

    Returns:
        float: Encirclement score (number of contacts).

    """
    from validation.brute_force_folding import BruteForceFolding

    solver = BruteForceFolding(protein, interaction, particle, dimension=2)
    result = solver.solve()
    config = result["configuration"]

    # Count contacts with particle (simplified)
    contacts = 0
    if config and "particle_position" in config:
        # Placeholder: count beads near particle
        contacts = sum(1 for i in range(len(protein.main_chain)))

    return float(contacts)


def objective_minimize_energy(
    protein: Protein, interaction: Interaction, particle: Particle
) -> float:
    """Objective: minimize total system energy (maximize stability).

    More negative energy = better.

    Args:
        protein (Protein): Candidate protein.
        interaction (Interaction): Interaction model.
        particle (Particle): Target particle.

    Returns:
        float: Negative of minimum energy (for maximization).

    """
    from validation.brute_force_folding import BruteForceFolding

    solver = BruteForceFolding(protein, interaction, particle, dimension=2)
    result = solver.solve()
    return -result["energy"]  # Negate to maximize


def objective_stabilization_effect(
    protein: Protein, interaction: Interaction, particle: Particle
) -> float:
    """Objective: maximize energy stabilization from particle presence.

    Function returns the energy difference: E_no_particle - E_with_particle.
    Positive values = stabilization.

    Args:
        protein (Protein): Candidate protein.
        interaction (Interaction): Interaction model.
        particle (Particle): Target particle.

    Returns:
        float: Stabilization energy (positive = stabilizing).

    """
    from validation.brute_force_folding import BruteForceFolding

    # Without particle
    solver_baseline = BruteForceFolding(protein, interaction, None, dimension=2)
    result_baseline = solver_baseline.solve()

    # With particle
    solver_particle = BruteForceFolding(protein, interaction, particle, dimension=2)
    result_particle = solver_particle.solve()

    stabilization = result_baseline["energy"] - result_particle["energy"]
    return stabilization
