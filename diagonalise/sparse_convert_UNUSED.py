#!/usr/bin/env python

import sys, os
import time
import numpy        as     np
import h5py
from ase.io.vasp import read_vasp


# input 
primitive_filename      = '/mnt/scratch2/q13camb_scratch/adps2/output_folder_chik5001/relax_output/20240320212733/relaxed_structure_starting_point.POSCAR' #= "POSCAR_5184"

input_filename= 'fc2.hdf5'


atoms = read_vasp(primitive_filename)

at_positions=atoms.get_scaled_positions()
nat=len(at_positions)
nat3=3*nat

masses = np.array(atoms.get_masses())


print(nat3)

sys.stdout.flush()

with h5py.File(input_filename, 'r' ) as fin:
  force_constants=fin['fc2'][()]
  #force_constants=fin['force_constants'][()]

#force_constants=np.transpose(force_constants,axes=(0,2,1,3))
print(force_constants.shape)


#fc2_abs=np.abs(force_constants)
mask=force_constants!=0#fc2_abs>0#5e-2
fc2_num=np.sum(mask)



sys.stdout.flush()

dynmat=force_constants/np.sqrt((masses.reshape(-1,1)*masses.T.reshape(1,-1)).reshape(nat,nat,1,1))

#TODO the above line seems to be dodgy

print(np.sum(mask,axis=1)[0])

print('Reading finished, Writing starting for %d values.'%(fc2_num))

sys.stdout.flush()


i,j,k,l=np.meshgrid(np.arange(dynmat.shape[0]),np.arange(dynmat.shape[1]),np.arange(dynmat.shape[2]),np.arange(dynmat.shape[3]),indexing='ij')

tosave=np.stack((i[mask]*3+k[mask]+1,j[mask]*3+l[mask]+1,dynmat[mask]),axis=1)

print(tosave.shape)

np.savetxt("dmat.dat",tosave,fmt='%d %d %20.15f')

print("done")


