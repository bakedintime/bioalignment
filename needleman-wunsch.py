import glob
import nwalign as nw

fasta_files = glob.glob(('fasta_data/*'))
seqs = []

# Read & Format
for file_name in fasta_files:
  f = open(file_name, "r")
  print file_name
  flines = f.readlines()
  seq = ""
  for el in flines[1:]:
    seq += el.strip()
  seqs.append(seq) # Place desired length of comparison here.
  f.close()

result = nw.global_align(seqs[0], seqs[1], matrix='BLOSUM62.txt')
score = nw.score_alignment(result[0], result[1], matrix='BLOSUM62.txt', gap_extend=-2, gap_open=-5)
print result
print score

