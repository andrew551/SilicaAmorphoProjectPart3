from mpi4py import MPI
from ase.atoms import Atoms
from ase.neighborlist import NeighborList, NewPrimitiveNeighborList
from ase.calculators.lammpslib import LAMMPSlib
from ase.io.vasp import read_vasp
import numpy as np
from sys import stdout
import h5py 
from os import makedirs,system
import pickle
from time import time
import lammps as lmp
import sys
import os
sys.path.insert(1, os.path.join(sys.path[0], '..')) # "hack" to add the parent directory to path
import makeconfig # this line needs to be after the sys.path.insert line
config = makeconfig.config()
import file_conversion

atoms = file_conversion.read_reg('relaxed_model/POSCAR')
#igenvectors = np.loadtxt('3_diagonalise/eigenvectors.dat')



displacement = 0.05


cmds = makeconfig.get_potential_command(config).split('\n')
atom_dict = makeconfig.get_atom_dict(config)
lammps_log_file = None

me = MPI.COMM_WORLD.Get_rank()
nprocs = MPI.COMM_WORLD.Get_size()
split=MPI.COMM_WORLD.Split(me,key=0)

print(me, nprocs, flush=True)

if me == 0:
    eigenvectors = np.loadtxt('3_diagonalise/eigenvectors.dat').T
    n = eigenvectors.shape[0]
    ev3 = eigenvectors.reshape((n, n//3, 3))
    print(eigenvectors[0, :10])
    print(eigenvectors[100, :10])
    print(ev3[0, :10, :])
    print(ev3[1000, :10, :])
    print('read eigenvectors of shape', eigenvectors.shape, flush=True)

else:
    eigenvectors = None
    ev3 = None

## TODO: eigenvectors scale by inverse square mass of each atom for proper movement

eigenvectors=MPI.COMM_WORLD.bcast(eigenvectors,root=0)
n = eigenvectors.shape[0]
ev3=MPI.COMM_WORLD.bcast(ev3,root=0)


lmp_test=lmp.lammps()
if lmp_test.has_mpi_support:
    comm_lammps=split
else:
    comm_lammps=None
lmp_test.close()
del lmp_test

lammps = LAMMPSlib(lmpcmds=cmds,log_file=lammps_log_file,atom_types=atom_dict, keep_alive='TRUE',comm=comm_lammps)
atoms.set_calculator(lammps)
MPI.COMM_WORLD.Barrier()

def compute_u0():
    # set calculators
    u = atoms.get_potential_energy()
    print('energy 0: ', u, flush=True)
    return u

def compute_element(i, j, k, l):
    s1 = -1 if k == 0 else 1
    s2 = -1 if l == 0 else 1
    atoms.positions += s1 * ev3[i, :, :] * displacement
    atoms.positions += s2 * ev3[j, :, :] * displacement
    u = atoms.get_potential_energy()
    atoms.positions -= s1 * ev3[i, :, :] * displacement
    atoms.positions -= s2 * ev3[j, :, :] * displacement
    
    print(f'u {i} {j} {k} {l} = ', u, flush=True)
    return u

def compute_element_single(i, k):
    s1 = -1 if k == 0 else 1
    atoms.positions += s1 * ev3[i, :, :] * displacement
    u = atoms.get_potential_energy()
    atoms.positions -= s1 * ev3[i, :, :] * displacement
    
    print(f'u{i} {k} = ', u, flush=True)
    return u

u0 = compute_u0()
compute_element(5, 5, 0, 0)

n=64


#elements = np.zeros(n, n, 2, 2)

##### shared memory output

comm = MPI.COMM_WORLD 
shape_out = (n, n, 2, 2)
size = np.prod(shape_out)
itemsize = MPI.DOUBLE.Get_size() 
if comm.Get_rank() == 0: 
    nbytes = size * itemsize 
else: 
    nbytes = 0

# on rank 0, create the shared block
# on rank 1 get a handle to it (known as a window in MPI speak)
win = MPI.Win.Allocate_shared(nbytes, itemsize, comm=comm) 

# create a numpy array whose data points to the shared mem
buf, itemsize = win.Shared_query(0) 
assert itemsize == MPI.DOUBLE.Get_size() 
elements = np.ndarray(buffer=buf, dtype='d', shape=shape_out) 


shape_out2 = (n, 2)
size2 = np.prod(shape_out)
itemsize2 = MPI.DOUBLE.Get_size() 
if comm.Get_rank() == 0: 
    nbytes2 = size2 * itemsize2
else: 
    nbytes2 = 0

# on rank 0, create the shared block
# on rank 1 get a handle to it (known as a window in MPI speak)
win2 = MPI.Win.Allocate_shared(nbytes2, itemsize2, comm=comm) 

# create a numpy array whose data points to the shared mem
buf2, itemsize2 = win2.Shared_query(0) 
assert itemsize2 == MPI.DOUBLE.Get_size() 
elements2 = np.ndarray(buffer=buf2, dtype='d', shape=shape_out2) 


######


if __name__ == '__main__':
    for i in range(n):
        if not i % nprocs == me:
            continue
        for j in range(n):
            if j > i:
                continue
            #if not j 
            for k in range(2):
                for l in range(2):
                    elements[i, j, k, l] = compute_element(i, j, k, l)
        for k in range(2):
            elements2[i, k] = compute_element_single(i, k)
MPI.COMM_WORLD.Barrier()

def restart_lammps(LAMMPSLibObject):
    if LAMMPSLibObject.started:
        LAMMPSLibObject.lmp.command("clear")
    # hope there's no other state to be reset
    LAMMPSLibObject.started = True
    LAMMPSLibObject.initialized = True
    LAMMPSLibObject.previous_atoms_numbers = []

def close_lammps(lammps):
    if lammps.lmp is not None:
        lammps.lmp.close()

close_lammps(lammps)

if me == 0:
    with open('10_quartic/elements.npy', 'wb') as f:
        np.save(f, elements)
    with open('10_quartic/elements2.npy', 'wb') as f:
        np.save(f, elements2)
    with open('10_quartic/u0.txt', 'w') as f:
        f.write(f'{u0}\n')
    #np.savetxt('10_quartic/elements.txt', elements)
    #np.savetxt('10_quartic/elements.txt', elements)

#np.savetxt(f'10_quartic/elements{me}.txt', elements)

MPI.Finalize()