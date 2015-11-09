import sys
import numpy as np
import random
import math
from Bio.PDB import *


class FileIO:
    def __init__(self, inFileName, outFileName):
        self.inFileName = inFileName
        self.outFileName = outFileName

    def readFile(self):
        inFile = open(self.inFileName)
        text = inFile.read()
        return text

    def writeFile(self, text):
        outFile = open(self.outFileName, 'w')
        outFile.write(text)


class NpFileIO(FileIO):
    def readFile(self):
        data = np.genfromtxt(self.inFileName, skip_header=1)
        return data


class ProcePDB:
    radii = {'N': 1.625,'CA': 1.9,'C': 1.875,'O': 1.48,'CB': 1.952,'CG': 1.952,'CG1': 1.952,'CG2': 1.955,
             'CD': 1.952,'CD1': 1.952,'CD2': 1.875,'CE': 1.952,'CE1': 1.875,'CE2': 1.875,'CE3': 1.875,
             'CZ': 1.875,'CZ2': 1.875,'CZ3': 1.875,'CH2': 1.875,'SG': 1.775,'SD': 1.775,'OG': 1.535,'OG1': 1.535,
             'OD1': 1.48,'OD2': 1.48,'OE1': 1.48,'OE2': 1.48,'ND1': 1.625,'ND2': 1.625,'NE': 1.625,'NE1': 1.625,
             'NE2': 1.625,'NZ': 1.625,'NH1': 1.625,'NH2': 1.625,'OH': 1.535
             }
    def __init__(self,fileName):
        pdbl = PDBList()
        #todo: what's this?
        pdbl.retrieve_pdb_file('1EN2') # arg is the pdb id

        parser = PDBParser()# will auto-correct obvious errors in pdb file
        #todo: filename not hard code
        self.structure = parser.get_structure('test','test.pdb')
        atom_list = []
        [atom_list.append(np.append(atom.get_coord(),self.radii[ atom.get_name()]).tolist()) for atom in self.structure.get_atoms()]
        self.atoms = np.asarray(atom_list)

    def get_atoms(self):
        return self.atoms


class MonteCarlo:
    def __init__(self, data):
        self.data = data

    def volume(self, N):
        atoms_data = self.data

        x_max = np.max(atoms_data[:, 0] + atoms_data[:, 3])
        x_min = np.min(atoms_data[:, 0] - atoms_data[:, 3])
        length = x_max - x_min
        # print x_min,x_max

        y_max = np.max(atoms_data[:, 1] + atoms_data[:, 3])
        y_min = np.min(atoms_data[:, 1] - atoms_data[:, 3])
        width = y_max - y_min
        # print y_min,y_max

        z_max = np.max(atoms_data[:, 2] + atoms_data[:, 3])
        z_min = np.min(atoms_data[:, 2] + atoms_data[:, 3])
        height = z_max - z_min
        # print z_min, z_max
        c = 0
        for i in range(N):
            random_x = random.uniform(x_min, x_max)
            random_y = random.uniform(y_min, y_max)
            random_z = random.uniform(z_min, z_max)
            dis = (random_x - atoms_data[:, 0]) ** 2 + (random_y - atoms_data[:, 1]) ** 2 + (random_z - atoms_data[:,
                                                                                                        2]) ** 2
            # print 'dis: ',dis
            d = dis - atoms_data[:, 3] ** 2
            # print 'd: ',d
            if len([x for x in d if x < 0]) > 0:
                c += 1
        # print c
        # print length*width*height
        c = float(c)
        N = float(N)
        v = c / N * length * width * height
        sd = length * width * height * math.sqrt((c / N - (c / N) ** 2) / N)
        # print v
        return [v, sd]

    def surface_area(self,N):
        atoms_data = self.data
        #todo: surface area:
        for i in xrange(N):
            phi = random.uniform(0,2*math.pi)
            costheta = random.uniform(-1,1)
            theta = math.acos(costheta)


if __name__ == '__main__':
    if len(sys.argv) == 2:
        inFile = sys.argv[1]
        outFile = 'result_' + inFile

        fileStream = NpFileIO(inFile, outFile)

        atoms_data = fileStream.readFile()
        #
        monteCarlo = MonteCarlo(atoms_data)
        #
        # result = map(monteCarlo.volume,[100,1000,10000])
        # print result

        pdbFile = ProcePDB(inFile)
        atoms = pdbFile.get_atoms()

        monteCarlo = MonteCarlo(atoms)
        result = map(monteCarlo.volume,[100,1000,10000])
        print result

    else:
        print 'command not correct!'
