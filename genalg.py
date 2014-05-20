# -*- coding: cp1252 -*-
# Basado en https://github.com/Murodese/Python-Genetic
import random
import logging
import sys
import pprint

logger = logging.getLogger(__name__)
logger.addHandler(logging.StreamHandler(sys.stdout))
logger.setLevel(logging.INFO)

debug = False

class DNA(object):
    """Clase con los genes y funciones"""
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
        """Crear secuencia y aleatorizar"""
        self.genes = []
        for index in xrange(0, self.size):
            if self.dataset:
                if self.dataset[index] == self.placeholder:
                    self.genes.insert(index, self.dataset[index])
            else:
                self.genes.insert(index, self.generation_function(self))

    def evaluate(self):
        """Evaluar en base a la funcion de fitness."""
        self.fitness = self.evaluation_function(self)

    def crossover(self, partner):
        """Genera hijos a través de midpoint aleatorio."""
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
        """Mutacion aleatoria"""
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
    """Poblacion que consiste en un conjunto de DNAs."""
    def __init__(self, size):
        self.size = size

        self.population = []
        self.mating_pool = []

    def setup(self, generation_function, evaluation_function,
        mutation_function, gene_length = 0, dataset = [], placeholder = None):


        self.population = []
        for index in xrange(0, self.size):
            self.population.insert(index,
                DNA(generation_function, evaluation_function,
                    mutation_function, gene_length, dataset, placeholder)
                )

        if debug:
            logger.debug('Population: ' + str(pprint.pformat(self.population)))

    def evaluate(self, target_fitness = 1.0):
        """Evalúa el fitness de cada DNA."""

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
        """Creación de generación."""

        self.mating_pool = []
        for dna in self.population:
            proportion = int((dna.fitness / 1.0) * 100)
            for index in xrange(0, proportion):
                self.mating_pool.insert(index, dna)

        if debug:
            logger.debug('Mating pool: ' + str(pprint.pformat(self.mating_pool)))

    def mate(self):
        """Cruzamiento. La probabilidad de cruce es más alta para genes con mayor fitness."""

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
        """Generación de estadísticas del proceso."""

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
        """Corre el algoritmo hasta alcanzar el fitness deseado."""

        generation = 1
        while 1:
            dna = self.evaluate(target_fitness)
            if dna:
                return (generation, dna)
            self.spawn()
            self.mate()

            generation += 1
