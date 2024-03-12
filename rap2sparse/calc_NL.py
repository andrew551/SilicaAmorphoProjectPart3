from ase.atoms import Atoms
from ase.neighborlist import NeighborList, NewPrimitiveNeighborList
from ase.io.vasp import read_vasp
import pickle
import os
import sys
sys.path.insert(1, os.path.join(sys.path[0], '..')) # "hack" to add the parent directory to path
import makeconfig # this line needs to be after the sys.path.insert line
config = makeconfig.config()

POSCAR_file='relaxed_structure_starting_point.POSCAR'
output_ext='.NL'
FC2_cutoff = config['[fc2]_FC2_cutoff'] #used to be 6
Force_cutoff = config['[fc2]_Force_cutoff']

output_file=POSCAR_file+str(FC2_cutoff)+str(Force_cutoff)+output_ext
print(output_file)

print("Neighbourlist calculation starting.")

atoms = read_vasp(POSCAR_file)

NL_fc2= NeighborList(FC2_cutoff, skin=0.3, sorted=True, self_interaction=True, bothways=True,  primitive=NewPrimitiveNeighborList)
NL_fc2.update(atoms)

NL_force= NeighborList(Force_cutoff, skin=0.3, sorted=True, self_interaction=True, bothways=True,  primitive=NewPrimitiveNeighborList)
NL_force.update(atoms)

print("Calculation finished. Writeout starting.")

with open(output_file,"wb") as f:
    pickle.dump(NL_fc2,f,pickle.HIGHEST_PROTOCOL)
    pickle.dump(NL_force,f,pickle.HIGHEST_PROTOCOL)

print("Neighbourlist saving finished.")
