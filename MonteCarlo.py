import sys
import numpy as np
import random
import math
from Bio.PDB import *
import time


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
        parser = PDBParser()# will auto-correct obvious errors in pdb file
        self.structure = parser.get_structure('test',fileName)
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
        z_min = np.min(atoms_data[:, 2] - atoms_data[:, 3])
        height = z_max - z_min
        # print z_min, z_max
        c = 0
        for i in xrange(N):
            random_x = random.uniform(x_min, x_max)
            random_y = random.uniform(y_min, y_max)
            random_z = random.uniform(z_min, z_max)
            dis = (random_x - atoms_data[:, 0]) ** 2 + (random_y - atoms_data[:, 1]) ** 2 + (random_z - atoms_data[:,
                                                                                                        2]) ** 2
            # print 'dis: ',dis
            d = dis - atoms_data[:, 3] ** 2
            # print 'd: ',d
            if min(d)<=0:
                c += 1
            # if len([x for x in d if x < 0]) > 0:
            #     c += 1

        c = float(c)
        N = float(N)
        # print c
        #
        # print length*width*height
        # print c/N
        v = c / N * length * width * height
        sd = length * width * height * math.sqrt((c / N - (c / N) ** 2) / N)
        # print v
        return [int(N), round(v,3), round(sd,3)]

    def surfatial_area(self,N):
        atoms_data = self.data
        surfatial_area = 0
        total_area = 0

        for i in xrange(len(atoms_data)):
            x0 = float(atoms_data[i][0])
            y0 = float(atoms_data[i][1])
            z0 = float(atoms_data[i][2])
            r  = float(atoms_data[i][3])
            c = 0
            for j in xrange(N):
                phi = random.uniform(0,360)
                costheta = random.uniform(-1,1)

                theta = math.acos(costheta)

                xp = x0 + r*math.sin(theta)*math.cos(phi)
                yp = y0 + r*math.sin(theta)*math.sin(phi)
                zp = z0 + r*costheta
                available = 1

                dis = (xp - atoms_data[:, 0]) ** 2 + (yp - atoms_data[:, 1]) ** 2 + (zp - atoms_data[:,2]) ** 2
                d = dis - atoms_data[:, 3] ** 2
                # d = np.delete(d,i,0)

                if min(d) < -0.000000001:
                    available = 0

                c += available

            surfatial_area += 4 * math.pi * r**2 * float(c)/float(N)
            total_area += 4 * math.pi * r**2

        return [int(N),round(surfatial_area,3)]


if __name__ == '__main__':
    if len(sys.argv) == 2:
        inFile = sys.argv[1]
        outFile = 'result_' + inFile

        pdbFile = ProcePDB(inFile)
        atoms = pdbFile.get_atoms()

        monteCarlo = MonteCarlo(atoms)

        t= time.time()
        volume = map(monteCarlo.volume,[10,100,1000,10000,100000])
        t=time.time()-t
        print 'results of volume:'
        print 'N'.rjust(7),
        for i in xrange(len(volume)):
            print str( volume[i][0]).rjust(12),
        print
        print 'volume'.rjust(7),
        for i in xrange(len(volume)):
            print str( volume[i][1]).rjust(12),
        print
        print 'SD'.rjust(7),
        for i in xrange(len(volume)):
            print str( volume[i][2]).rjust(12),
        print

        # print t
        print '========================================'
        t = time.time()
        surfatial = map(monteCarlo.surfatial_area,[5,10,20,50,100])
        t = time.time()-t
        print 'results of surfatial area: '
        print 'N'.rjust(9),
        for i in xrange(len(surfatial)):
            print str( surfatial[i][0]).rjust(12),
        print
        print 'srf_area'.rjust(9),
        for i in xrange(len(surfatial)):
            print str( surfatial[i][1]).rjust(12),
        print

        # print t
    else:
        print 'command not correct!'
        print 'python MonteCarlo.py <file_name>'
        print 'file_name needs end with .pdb'
