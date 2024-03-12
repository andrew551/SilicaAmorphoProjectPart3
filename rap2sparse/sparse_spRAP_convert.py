#!/usr/bin/env python

import sys, os
import time
import numpy        as     np
import h5py
from ase.io.vasp import read_vasp


# input 
primitive_filename       = "relaxed_structure_starting_point.POSCAR"

input_filename= 'fc2.hdf5'


atoms = read_vasp(primitive_filename)

at_positions=atoms.get_scaled_positions()
nat=len(at_positions)
nat3=3*nat

masses = np.array(atoms.get_masses())


print(nat3)

sys.stdout.flush()

with h5py.File(input_filename, 'r' ) as fin:
  fc2=fin['fc2'][()]
  #fc2=fin['fc2'][()]

#fc2=np.transpose(fc2,axes=(0,2,1,3))
print(fc2.shape)
inds=fc2[:,:4].astype(np.int64)

sys.stdout.flush()

dynmat=fc2[:,4]/np.sqrt((masses[inds[:,0]]*masses[inds[:,1]]))


print('Reading finished, Writing starting for %d values.'%(fc2.shape[0]))

sys.stdout.flush()



tosave=np.stack((inds[:,0]*3+inds[:,2]+1,inds[:,1]*3+inds[:,3]+1,dynmat),axis=1)

print(tosave.shape)

np.savetxt("dmat.dat",tosave,fmt='%d %d %20.15f')

print("done")