#!/usr/bin/env python

import sys, os
import time
import numpy        as     np
import h5py
# import phonopy
from ase.io.vasp import read_vasp
#from phono3py.phonon.velocity_operator import VelocityOperator 
#from ase.io                   import read
#from ase.units                import J, _hplanck  # conversion factors
#from ase.dft.kpoints          import bandpath

global start_time


save_single=False
withBORN=False
#def write__to_hdf5(filename, frequencies, qpoints, weights, velocity_operators):
#	compression="gzip"
#	with h5py.File(filename, "w") as w:
#		w.create_dataset("frequency", data=frequencies, compression=compression)
#		w.create_dataset("qpoint", data=qpoints)
#		w.create_dataset("weight", data=weights)
#		w.create_dataset("velocity_operator", data=velocity_operators, compression=compression)
#
#start_time = time.time()
#eV2Hz  = 1/ ( J * _hplanck )

# input 
primitive_filename       = "./POSCAR"
input_filename= 'fc2.hdf5'
supercell                = np.diag((1,1,1))
mesh                     = [1,1,1]
interpolation_mesh_size=mesh[0]


atoms = read_vasp(primitive_filename)

at_positions=atoms.get_scaled_positions()
nat=len(at_positions)
nat3=3*nat

print(nat3)

fin   =h5py.File(input_filename, 'r' )
force_constants=fin['fc2'][()]
fin.close()

if np.max(np.isnan(force_constants)):
  print("nan")

print('Reading finished, Writing starting.')

f = open('mat2d.dat', 'w')
for i in range(0, nat):
  for j in range(0, nat):
    for k in range(0, 3):
      for l in range(0, 3):
        f.write("%i %i %20.15f\n" % (i*3+k+1, j*3+l+1, force_constants[i][j][k][l]))
f.close()
 
print("done")