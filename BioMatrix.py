import pymongo
import numpy as np
import pickle
import bson

class BioMatrix(object):
    def __init__(self):
        self.biomatrix = None
        self.contigs = []
        self.total_contigs = 0
        self.score = 2147483648
        self.pairs = []

    def init_db(self):
        self.conn = pymongo.Connection("mongodb://bioa:bioa@oceanic.mongohq.com:10091/bioalignment")
        self.db = self.conn.bioalignment
        self.matrix = self.db["matrix"]

    def get_score(self):
        return self.score

    def get_matrix(self):
        return self.biomatrix

    def post_matrix(self, bmatrix):
        self.biomatrix = bmatrix

    def init_contigs(self, total_contigs):
        for i in range(0, total_contigs):
            self.contigs.append(0)
        self.total_contigs = total_contigs

    def init_matrix(self):
        self.biomatrix = np.zeros([self.total_contigs, self.total_contigs])
        for i in range(0, self.total_contigs):
            self.biomatrix[i,i] = 2147483648

    def get_connection(self):
        return self.conn

    def get_from_db(self):
        if not self.conn:
            print 'Primero debes ejecutar init_db()'
        else:
            matrix = self.matrix.find_one({"key":"main"})
            deserialized_matrix = pickle.loads(matrix['matrix'])
            self.biomatrix = deserialized_matrix
            return self.biomatrix

    def post_to_db(self, new_matrix):
        if not self.conn:
            print 'Primero debes ejecutar init_db()'
        else:
            array = self.matrix.find_one({"key":"main"})
            if array:
                array['matrix'] = bson.binary.Binary(pickle.dumps(new_matrix, protocol=2))
                self.matrix.save(array)
            else:
                self.matrix.insert({
                    'key':'main',
                    'matrix': bson.binary.Binary(pickle.dumps(new_matrix, protocol=2))
                })

    # Matrix traversal algorithms

    def get_min_alignment(self, row, idx):
        # indice del menor valor por fila
        m = np.zeros_like(self.biomatrix)
        m[:, idx] = 1
        masked = np.ma.masked_array(self.biomatrix, m)
        value = np.argmin(masked, axis=1)[row]
        return value

    def get_aligned_pairs(self):
        self.score = 0
        for i in range(0, self.total_contigs):
            j = self.get_min_alignment(i, np.array([self.contigs.index(x) for x in self.contigs if (x == 1 or x == i)]))
            if self.contigs[j] != 0:
                print 'El segmento ',i ,' no se ha podido emparejar o ya fue emparejado.'
            else:
                self.contigs[i] = 1
                self.contigs[j] = 1
                self.pairs.append([i,j])
                self.score += self.biomatrix[i,j]
                print 'El segmento', i, 'se empareja con el segmento ', j
        return self.contigs


"""
Ejemplo random de la agrupacion de contigs
y su ponderacion a traves de una matriz

if __name__=='__main__':
    biom = BioMatrix()
    biom.init_contigs(4)
    biom.init_matrix()
    m =  biom.get_matrix()
    m[0,1] = 34
    m[1,0] = 34
    m[0,2] = 12
    m[2,0] = 12
    m[0,3] = 5
    m[3,0] = 5
    m[1,2] = 14
    m[2,1] = 14
    m[1,3] = 7
    m[3,1] = 7
    m[2,3] = 50
    m[3,2] = 50
    biom.post_matrix(m)
    print biom.get_matrix()
    print biom.get_aligned_pairs()
    print biom.get_score()
"""
