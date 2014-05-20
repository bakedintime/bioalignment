import glob, os, random
import pymongo
from GenAlgorithm import DNASegment, Population
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
        self.population = None

    def load_sequences(self):
        for f1 in self.fasta_files:
            seq = self.fastaReader.get_sequence(f1)
            self.sequences.append([seq, f1])

    def build_genetic_alg_description(self):
        def gen_func(dna):
            """Used to initialise DNA. Should return a single gene - in this case,
            a letter"""
            return [0, 1]

        def eval_func(dna):
            """Used to evaluate fitness of DNA. Should return a fitness value, 0-1
            where higher is better."""
            score = 0
            genes = dna.genes
            # matrix coordinates
            i = genes[0].get_just_filename()
            j = genes[1].get_just_filename()
            biom = BioMatrix()
            score = biom.get_coordinate_score(i, j)
            percentage = float(score*100)/len(genes[0].get_sequence())
            fitness = 1 - float(percentage/100)
            return fitness

        def mut_func(dna):
            """Used to mutate a gene. Should return a single gene."""
            index1 = random.randint(0, len(self.genes))
            gene = self.genes[index1]
            dna_seg1 = gene[0]
            dna_seg2 = gene[1]
            if '-' in dna_seg1.get_sequence():
                seq = dna_seg1.get_sequence.strip('-')
                diff = len(dna_seg1.get_sequence()) - len(dna_seg2.get_sequence())
                bufferPlaces1 = []
                for i in range(0, diff):
                    index = random.randint(0, len(seq))
                    seq = seq[:index] + '-' + seq[index:]
                    bufferPlaces1.append(index)
                dna_segment1 = DNASegment(seq, dna_seg1.filename, bufferPlaces1)
                new_gene = [dna_segment1, dna_seg2] 
            else:
                seq = dna_seg2.get_sequence.strip('-')
                diff = len(dna_seg1.get_sequence()) - len(dna_seg2.get_sequence())
                bufferPlaces2 = []
                for i in range(0, diff):
                    index = random.randint(0, len(seq))
                    seq = seq[:index] + '-' + seq[index:]
                    bufferPlaces2.append(index)
                dna_segment2 = DNASegment(seq, dna_seg2.filename, bufferPlaces2)
                new_gene = [dna_seg1, dna_segment2] 
            return new_gene

        combined_sequences = []

        for f1 in self.sequences:
            for f2 in self.sequences:
                if f1[1] != f2[1]:
                    seq1 = f1[0]
                    seq2 = f2[0]
                    diff = len(seq1) - len(seq2)
                    # insertar dashes de manera random
                    bufferPlaces1=[]
                    bufferPlaces2=[]
                    if diff > 0:
                        for i in range(0, diff):
                            index = random.randint(0, len(seq2))
                            seq2 = seq2[:index] + '-' + seq2[index:]
                            bufferPlaces2.append(index)
                    elif diff < 0:
                        for i in range(0, diff):
                            index = random.randint(0, len(seq1))
                            seq1 = seq1[:index] + '-' + seq1[index:]
                            bufferPlaces1.append(index)
                    filename1 = f1[1]
                    filename2 = f2[1]
                    dna_segment1= DNASegment(seq1, filename1, bufferPlaces1)
                    dna_segment2 = DNASegment(seq2, filename2, bufferPlaces2)
                    combined_sequences.append([dna_segment1, dna_segment2])


        self.population = Population(len(combined_sequences))
        self.population.setup(
            gen_func,
            eval_func,
            mut_func,
            # i**2 -i 
            dataset = combined_sequences
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
    # Carga a memoria los contigs
    # para que puedan ser modificados
    # por el algoritmo genetico
    mw.load_sequences()
    mw.build_genetic_alg_description()
    
    def recompute_fitness(mw):
        mw.delegate_jobs()
        # Esperar hasta terminar los calculos
        while True:
            if (mw.check_if_finished()):
                mw.get_results()
                break
        
    (generation, dna) = mw.population.run(recompute_fitness, mw, target_fitness=1.0)
    print('Finished!')
    print('Generation: ' + str(generation))
    #print('DNA: ' + str(dna.genes[0].bufferPlaces) + ', '+ str(dna.genes[1].bufferPlaces))
    print dna.genes
    print 'ya'
