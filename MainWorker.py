import glob, os
from BioMatrix import BioMatrix
from DistributedWorker import compute_two_sequences
from needleman_wunsch import FastaReader

class MainWorker(object):
    def __init__(self, overwrite=False):
        self.fasta_files = glob.glob(('fasta_data/*.fasta'))
        # Inicializar la matriz y subirla a mongodb
        self.biom = BioMatrix()
        self.biom.init_db()
        if overwrite:
            self.biom.init_matrix()
            m = self.biom.get_matrix()
            self.biom.post_to_db(m)
        else: 
            self.biom.get_from_db()
        self.biom.init_contigs(len(self.fasta_files))
        self.fastaReader = FastaReader()

    def delegate_jobs(self):
        # total len(self.fasta_files)**2 - len(self.fasta_files)
        for f1 in self.fasta_files:
            for f2 in self.fasta_files:
                if f1 != f2:
                    seq1 = self.fastaReader.get_sequence(f1)
                    seq2 = self.fastaReader.get_sequence(f2)
                    result = compute_two_sequences.delay(seq1, seq2, os.path.basename(f1), os.path.basename(f2))
                    print result

    def get_results(self):
        print self.biom.get_from_db()
        print self.biom.get_aligned_pairs()
        print self.biom.get_score()

if __name__=="__main__":
    mw = MainWorker(overwrite=False)
    mw.delegate_jobs()
    #mw.get_results()
    print 'ya'