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
import matplotlib.pyplot as plt
#atoms = file_conversion.read_reg('relaxed_model/POSCAR')
#igenvectors = np.loadtxt('3_diagonalise/eigenvectors.dat')

config['ZBL_flag'] = True
p0 = np.array([[10, 0, 0], [-0.5, np.sqrt(3)/2, 0], [-0.5, -np.sqrt(3)/2, 0]])
atoms = Atoms('Si3', positions=p0, cell=[50,50,50], pbc=[1, 1, 1])

file_conversion.write_file(atoms, 'temp_test_atoms.vasp', out_type='vasp')

cmds = makeconfig.get_potential_command(config).split('\n')
atom_dict = makeconfig.get_atom_dict(config)
lammps_log_file = None

me = MPI.COMM_WORLD.Get_rank()
nprocs = MPI.COMM_WORLD.Get_size()
split=MPI.COMM_WORLD.Split(me,key=0)

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

def compute_u0(r):
    # set calculators
    atoms.set_positions(p0*r)
    u = atoms.get_potential_energy()
    print('energy 0: ', u, flush=True)
    return u


######

#offset = 15000
offset = 0
if __name__ == '__main__':
    species = ['Si3', 'SiO2', 'Si2O', 'O3']
    for specie in species:
        atoms.set_chemical_symbols(specie)
        rs = np.linspace(0.277, 6, 50)
        us = [compute_u0(r) for r in rs]
        plt.plot(rs, us, label = specie)
    plt.legend()
    plt.ylabel("u / eV")
    plt.xlabel("r / Angstroms")
    plt.title(f"ACE three-body test zbl-- {config['potential']}")
    plt.savefig("SiSi_threebody10rzbl.png", dpi=400)
    plt.clf()
    for specie in species:
        atoms.set_chemical_symbols(specie)
        rs = np.logspace(-2, 0, 100)
        us = [compute_u0(r) for r in rs]
        plt.plot(rs, us, label = specie)
    plt.legend()
    plt.loglog()
    plt.ylabel("u / eV")
    plt.xlabel("r / Angstroms")
    plt.title(f"ACE two-body test close zbl-- {config['potential']}")
    plt.savefig("twothreebody_10rclosezbl.png", dpi=400)

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


MPI.Finalize()