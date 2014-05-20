from celery import Celery
from BioMatrix import BioMatrix
from needleman_wunsch import NWAligner

app = Celery('tasks', broker='mongodb://bioa:bioa@oceanic.mongohq.com:10091/bioalignment')

@app.task
def compute_two_sequences(seq1, seq2, filename1, filename2):
    aligner = NWAligner()
    result = aligner.score_alignment(seq1, seq2)
    print 'Realizando la alineacion de secuencias correspondientes a los archivos: ['+filename1+', '+filename2+']'
    print 'Resultado: '+str(result)
    biom = BioMatrix()
    m = biom.get_from_db()
    x = filename1.split('.')[0]
    y = filename2.split('.')[0]
    if result < 0 or result == 0:
        m[int(x)-1, int(y)-1] = 2147483648
    else:
        m[int(x)-1, int(y)-1] = result
    biom.post_to_db(m)
    return 'Terminado'

