# Working with pdb files
# Adapted from Biopython docs
from Bio.PDB import *

# To fetch a file from the protein Data Bank
# File will be deposited off working directory in folder named from pdb id
# Ex: pdb file for 1EN2 is deposited as en/pdb1en2.ent:
pdbl = PDBList()
pdbl.retrieve_pdb_file('1EN2') # arg is the pdb id

# To work with a pdb file, begin by creating a parser:
parser = PDBParser()# will auto-correct obvious errors in pdb file
                    # if not desired use PERMISSIVE=0

# To read structure from file in current working directory
#  first arg is the name user-defined name assigned to structure. Ex:
structure = parser.get_structure('1IGT', '1IGT.pdb')
# structure = parser.get_structure('test','test.pdb')

# structure.header is a dictionary containing keys: name, head, deposition_date,
# release_date, structure_method, resolution, structure_reference (maps to a list of references),
# journal_reference, author and compound
resolution = structure.header['resolution']
keywords = structure.header['keywords']
print "resolution", resolution
print "keywords", keywords


# The Biopython PDB structure:
# A structure consists of models
# A model consists of chains
# A chain consists of residues
# A residue consists of atoms

# NB: with NMR there are several models, however,
#  with x-ray crystallography there is only model[0]

# for model in structure:
#     for chain in model:
#         for residue in chain:
#             for atom in residue:
#                 print atom.get_coord(), atom.get_name()  # type of coord is numpy.ndarray of float32
#                         # name of atom is used to get atomic radius from force field and other files
print "----------------------------------"
atom_list = []
for atom in structure.get_atoms():
    atom_list.append( atom.get_coord())
print atom_list
print "----------------------------------"


print
print
print  "************************************************"

# Ex: One way to iterate over all the chains in a model
for model in structure:
    for chain in model:
        print chain
print

# Ex: another way to iterate over all the chains in a model
for chain in structure.get_chains():
    print chain
print

## Ex: One way to iterate over all residues in a model
#for residue in model.get_residues():
#    print residue
#print

print  "************************************************"
# Drilling down into the structure for the entity you want:
print structure[0] # model
print structure[0]['A'] # add one more bracket operator for chain
res100 =  structure[0]['A'][100] # add one more bracket operator for residue
   # residue: use tuple to specify if there is an insertion code or if hetero-residue with same id exists
print "res100: ", res100
atom = structure[0]['A'][100]['CA']  # add one more bracket operator for atom
print atom
residue = atom.get_parent()
print residue
chain = residue.get_parent()
print chain


atoms = structure.get_atoms()
print atoms