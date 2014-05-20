import glob, os
import pymongo
from GenAlgorithm import DNA, Population
from BioMatrix import BioMatrix
from DistributedWorker import compute_two_sequences
from needleman_wunsch import FastaReader

class MainWorker(object):
    def __init__(self, overwrite=False):
        self.fasta_files = glob.glob(('fasta_data/*.fasta'))
        # Inicializar la matriz y subirla a mongodb
        self.biom = BioMatrix()
        self.biom.init_contigs(len(self.fasta_files))
        if overwrite:
            self.biom.init_matrix()
            m = self.biom.get_matrix()
            self.biom.post_to_db(m)
        else:
            self.biom.get_from_db()
        self.fastaReader = FastaReader()
        self.sequences = []

    def load_sequences(self):
        for f1 in self.fasta_files:
            seq = self.fastaReader.get_sequence(f1)
            self.sequences.append([seq, f1])

    def build_genetic_alg(self):
        #def gen_func(dna):
        #    """Used to initialise DNA. Should return a single gene - in this case,
        #    a letter"""
        #    return random.choice(string.letters + ' ,.?')

        def eval_func(dna):
            """Used to evaluate fitness of DNA. Should return a fitness value, 0-1
            where higher is better."""
            score = 0
            for index, gene in enumerate(dna.genes):
                if gene == target[index]:
                    score += 1
            fitness = score / float(len(target))
            return fitness

        def mut_func(dna):
            """Used to mutate a gene. Should return a single gene."""
            return random.choice(string.letters + ' ,.')

        target = "Is Kevin a giant faglord? Yes."

        self.population = Population(len(self.sequences))
        self.population.setup(
            #gen_func,
            eval_func,
            mut_func,
            gene_length = len(target),
            dataset = self.sequences
        )



    def check_if_finished(self):
        conn = pymongo.Connection("mongodb://bioa:bioa@oceanic.mongohq.com:10091/bioalignment")
        db = conn.bioalignment
        messages = db["messages"]
        if messages.find().count() == 0:
            return True
        else:
            return False

    def delegate_jobs(self):
        for f1 in self.sequences:
            for f2 in self.sequences:
                if f1[1] != f2[1]:
                    seq1 = f1[0]
                    seq2 = f2[0]
                    filename1 = f1[1]
                    filename2 = f2[1]
                    result = compute_two_sequences.delay(seq1, seq2, os.path.basename(filename1), os.path.basename(filename2))
                    print result

    def get_results(self):
        print 'Matriz resultante:'
        print self.biom.get_from_db()
        print 'Segmentos alineados: '
        print self.biom.get_aligned_pairs()
        print 'Fitness (Ponderacion de alineacion): ', self.biom.get_score()

if __name__=="__main__":
    mw = MainWorker(overwrite=True)
    mw.load_sequences()
    #for i in range(0,5):
    mw.delegate_jobs()
    # Esperar hasta terminar los calculos
    while True:
        if (mw.check_if_finished()):
            mw.get_results()
            break

        #(generation, dna) = self.population.run(1.0)
        #print('Finished!')
        #print('Generation: ' + str(generation))
        #print('DNA: ' + ''.join(dna.genes))
        #mw.get_results()
    print 'ya'
