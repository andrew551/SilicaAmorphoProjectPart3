from ase.atoms import Atoms
from ase.neighborlist import NeighborList, NewPrimitiveNeighborList
from ase.io.vasp import read_vasp
import pickle
import numpy as np
import sys
import os
sys.path.insert(1, os.path.join(sys.path[0], '..')) # "hack" to add the parent directory to path
import makeconfig # this line needs to be after the sys.path.insert line
config = makeconfig.config()

POSCAR_file='POSCAR'
output_ext='.NLC'
Force_cutoff = config["[4_fc3]force_cutoff"]
FC3_cutoff = config["[4_fc3]fc3_cutoff"]

output_file=POSCAR_file+str(FC3_cutoff)+str(Force_cutoff)+output_ext
print(output_file)

print("Neighbourlist calculation starting.")

atoms = read_vasp(POSCAR_file)

natoms=len(atoms)
cell_parameters=atoms.cell.cellpar()
min_cell_length=min(cell_parameters[:3])

NL_fc3= NeighborList(FC3_cutoff, skin=0.3, sorted=True, self_interaction=True, bothways=True,  primitive=NewPrimitiveNeighborList)
NL_fc3.update(atoms)

if min_cell_length>2.1*(Force_cutoff+FC3_cutoff):
    print("Force NL calculation.")

    NL_force= NeighborList(Force_cutoff+FC3_cutoff, skin=0.3, sorted=True, self_interaction=True, bothways=True,  primitive=NewPrimitiveNeighborList)
    NL_force.update(atoms)
else:
    NL_force=None
    print("Force cutoff is too large. No force Neighbour list is calculated.")


print("Calculation finished. Writeout starting.")

compList=np.zeros((len(NL_fc3.get_neighbors(0)[0])*natoms*20,2))
ind=0
for i in range(natoms):
    indices=NL_fc3.get_neighbors(i)[0]
    for j in indices:
        if i<=j:
            compList[ind]=np.array([i,j])
            ind+=1

compList=compList[:ind,:]
print(ind, "pairs to be computed.")


with open(output_file,"wb") as f:
    pickle.dump(NL_fc3,f,pickle.HIGHEST_PROTOCOL)
    pickle.dump(NL_force,f,pickle.HIGHEST_PROTOCOL)
    pickle.dump(compList,f,pickle.HIGHEST_PROTOCOL)
    

print("Neighbourlist saving finished.")
