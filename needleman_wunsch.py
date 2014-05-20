import nwalign as nw
import random

class FastaReader(object):
    def get_sequence(self, filename):
        # Read & Format
        #for file_name in fasta_files:
        f = open(filename, "r")
        flines = f.readlines()
        seq = ""
        for el in flines[1:]:
            seq += el.strip()
        f.close()
        return seq

class NWAligner(object):

    def align_sequences(self, seq1, seq2):
        result = nw.global_align(seq1, seq2, matrix='BLOSUM62.txt')
        return result

    def score_alignment(self, seq1, seq2, gp_e=-2, gp_o=-5):
        score = nw.score_alignment(seq1, seq2, matrix='BLOSUM62.txt', gap_extend=gp_e, gap_open=gp_o)
        return score
