#!/usr/bin/env python

import sys, os
import time
import numpy as np
import h5py
# import phonopy
from phonopy import Phonopy
from phonopy import file_IO
from phonopy.interface.calculator import read_crystal_structure
from phonopy.units import VaspToEv, VaspToCm, PwscfToTHz, VaspToTHz

global start_time

save_single = False
withBORN = False

# input
#primitive_filename = "POSCAR"
primitive_filename      = '/mnt/scratch2/q13camb_scratch/adps2/output_folder_chik5001/relax_output/20240320212733/relaxed_structure_starting_point.POSCAR' #= "POSCAR_5184"

supercell = np.diag((1, 1, 1))
mesh = [1, 1, 1]
interpolation_mesh_size = mesh[0]

atoms, string = read_crystal_structure(filename=primitive_filename)

at_positions = atoms.get_scaled_positions()
nat = len(at_positions)
nat3 = 3 * nat

print(nat3)
# phonons = Phonopy(atoms, supercell)

force_constants = file_IO.read_force_constants_hdf5(filename="fc2.hdf5")
masses = atoms.get_masses()

# add 15 additional fake atoms
fake_nat = nat+15

f = open('dmat.dat', 'a')
for i in range(0, fake_nat):
    for j in range(0, fake_nat):
        for k in range(0, 3):
            for l in range(0, 3):
                if (i >= nat) or (j >= nat):
                    if (i == j) and (k == l):
                        value = 50 + (i-nat)*10 + k*2
                    else:
                        value = 0
                    f.write("%i %i %20.15f\n" % (i * 3 + k + 1, j * 3 + l + 1, value))
                else:
                    value = force_constants[i][j][k][l] / np.sqrt(masses[i] * masses[j])
                    f.write("%i %i %20.15f\n" % (i * 3 + k + 1, j * 3 + l + 1, value))

print("done")


