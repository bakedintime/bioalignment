#author: Murodese
#link: github.com/Murodese

import random
import logging
import sys, os
import pprint

logger = logging.getLogger(__name__)
logger.addHandler(logging.StreamHandler(sys.stdout))
logger.setLevel(logging.INFO)

debug = True


class DNASegment(object):
    def __init__(self, contig_sequence, filename, bufferPlaces):
        self.sequence = contig_sequence
        self.filename = filename
        self.fitness = -1
        self.bufferPlaces = bufferPlaces

    def is_fitness_set(self):
        if self.fitness == -1:
            return False
        else:
            return True

    def get_just_filename(self):
        return os.path.basename(self.filename.split('.')[0])

    def get_sequence(self):
        return self.sequence

    def set_fitness(self, f):
        self.fitness = f

    def get_fitness(self):
        return self.fitness

    def __repr__(self):
        return self.filename

class DNA(object):
    """Python GA DNA sequence, containing genes and all associated
    functions"""
    def __init__(self, generation_function, evaluation_function,
        mutation_function, gene_length = 0, dataset = [], placeholder = None):
        self.dataset = dataset
        self.placeholder = placeholder

        if dataset:
            self.size = len(dataset)
        else:
            self.size = gene_length

        self.fitness = 0

        self.genes = []

        self.generation_function = generation_function
        self.evaluation_function = evaluation_function
        self.mutation_function = mutation_function

        self.setup()

    def setup(self):
        """Create the DNA sequence and randomise genes"""
        self.genes = []
        for index in xrange(0, self.size):
            if self.dataset:
                #if self.dataset[index] == self.placeholder:
                self.genes.insert(index, self.dataset[index])
            else:
                self.genes.insert(index, self.generation_function(self))

    def evaluate(self):
        """Fitness function. Will automatically call the evaluation function
        given during initialisation. Fitness values should be from 0 to 1,
        the higher the better."""
        self.fitness = self.evaluation_function(self)

    def crossover(self, partner):
        """Generates a child from the crossover of this object and a
        parent. Chooses a random midpoint to split genes."""
        child = DNA(self.generation_function, self.evaluation_function,
                    self.mutation_function, self.size, self.dataset,
                    self.placeholder)

        midpoint = random.choice(xrange(0, self.size))

        for index in xrange(0, self.size):
            if index > midpoint:
                child.genes[index] = self.genes[index]
            else:
                child.genes[index] = partner.genes[index]

        return child

    def mutate(self, rate=0.01):
        """Mutates genes randomly, given a rate and a mutation function"""
        for index in xrange(0, self.size):
            if random.random() < rate:
                if self.dataset:
                    if self.genes[index] == self.dataset[index]:
                        continue
                else:
                    self.genes[index] = self.mutation_function(self)

    def __repr__(self):
        return '[' + ','.join(self.genes) + ']'

class Population(object):
    """Holds a population of DNA, as well as functions for iterating over
    them."""
    def __init__(self, size):
        self.size = size

        self.population = []
        self.mating_pool = []

    def setup(self, generation_function, evaluation_function,
        mutation_function, gene_length = 0, dataset = [], placeholder = None):
        """Sets up the population, including generating pools of DNA.
        Set either gene_length (for a totally randomised dataset using
        generation_function), or import a dataset. Placeholders represent
        mutable values."""

        self.population = []
        for index in xrange(0, self.size):
            self.population.insert(index,
                DNA(generation_function, evaluation_function,
                    mutation_function, gene_length, dataset, placeholder)
                )

        if debug:
            logger.debug('Population: ' + str(pprint.pformat(self.population)))

    def evaluate(self, target_fitness = 1.0):
        """Evaluates the fitness of each DNA strand."""

        total_fitness = 0
        best_fitness = 0
        best_dna = None

        for dna in self.population:
            dna.evaluate()

            total_fitness += dna.fitness
            if dna.fitness > best_fitness:
                best_fitness = dna.fitness
                best_dna = dna

            if dna.fitness >= target_fitness:
                return dna

        logger.info('Average fitness: ' + str(total_fitness / self.size))
        logger.info('Best fitness:' + str(best_fitness))
        logger.info('Best DNA: ' + ''.join(best_dna.genes))

        return None

    def spawn(self):
        """Spawns mating pools based on DNA fitness, giving preferential
        treatment to higher-performing genes, but still giving fat, slow
        dudes a chance."""

        self.mating_pool = []
        for dna in self.population:
            proportion = int((dna.fitness / 1.0) * 100)
            for index in xrange(0, proportion):
                self.mating_pool.insert(index, dna)

        if debug:
            logger.debug('Mating pool: ' + str(pprint.pformat(self.mating_pool)))

    def mate(self):
        """Go through the population and force-mate random pairs. Because
        the spawning pools were based on fitness, this still prefers better
        genes, but will allow smaller chances of pairing with idiots."""

        for index in xrange(0, self.size):
            first = random.choice(self.mating_pool)
            second = random.choice(self.mating_pool)
            while first is second:
                # make sure that we're not trying to masturbate
                # ...first time I think I've ever written that...
                second = random.choice(self.mating_pool)

            child = first.crossover(second)
            child.mutate()

            self.population[index] = child

    def stats(self):
        """Generate various statistics on the population. Namely, who
        won the award for "best gene", and the average fitness. Useful,
        because it provides feedback for the user and gives us some idea
        of how far through the process we are."""

        total_fitness = 0
        best_performer = ''
        highest_fitness = 0
        for dna in self.population:
            total_fitness += dna.fitness

            if dna.fitness > highest_fitness:
                highest_fitness = dna.fitness
                best_performer = dna.fitness

        return best_performer, total_fitness / self.size

    def run(self, target_fitness = 1.0):
        """Runs the genetic algorithm until target_fitness is reached."""

        generation = 1
        while 1:
            dna = self.evaluate(target_fitness)
            if dna:
                return (generation, dna)
            self.spawn()
            self.mate()

            generation += 1