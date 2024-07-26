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
import math


atoms = Atoms('Si2', positions=[[-1, 0, 0], [1, 0, 0]], cell=[50,50,50], pbc=[1, 1, 1])

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

def get_lammps_for_potential():
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
    return lammps

lammps = LAMMPSlib(lmpcmds=cmds,log_file=lammps_log_file,atom_types=atom_dict, keep_alive='TRUE',comm=comm_lammps)
atoms.calc = lammps
MPI.COMM_WORLD.Barrier()

def compute_u0(r):
    # set calculators
    atoms.set_positions([[r/2, 0, 0], [-r/2, 0, 0]])
    try:
        u = atoms.get_potential_energy()
    except Exception:
        u = math.nan
    print('energy 0: ', u, flush=True)
    return u

def compute_bks(r, specie):
    if specie == 'Si2':
        qq = 2.4**2
        A = 0
        B = 0
        C = 0
        D = 3423200
    elif specie == 'SiO':
        qq = -1.2*2.4
        A = 18003.757
        B = 4.87318
        C = 133.5381
        D = 29
    elif specie == 'O2':
        qq = 1.2 ** 2
        A = 1388.773
        B = 2.76
        C = 175.0
        D = 113
    D = 0
    return qq * 14.3997 / r + A * np.exp(-B*r) - C / r**6 + D / r**24

# from Carre, et al paper: New effective potential for Silica    (2008)
def compute_CHIK(r, specie):
    if specie == 'Si2':
        qq = 1.910418**2
        A = 3150.462646
        B = 2.851451
        C = 626.751953
        D = 3423200
    elif specie == 'SiO':
        qq = -1.910418**2 / 2
        A = 27029.419922
        B = 5.158606
        C = 148.099091
        D = 29
    elif specie == 'O2':
        qq = 1.910418**2 / 4
        A = 659.595398
        B = 2.590066
        C = 26.836679
        D = 113
    else:
        raise Exception("invalid specie")
    return qq * 14.3997 / r + A * np.exp(-B*r) - C / r**6 + D / r**24

def test_all_potentials_dimer():
    #potentials = ['CHIK', 'ACE_Chuck', 'ACE_Deringher', 'ACE_KAMIL_5', 'ACE_12_0.5_0.0003', 'ACE_18_N2750', 'GAP_SiO2']
    potentials = ['CHIK', 'ACE_Chuck', 'ACE_CHEAP_CHUCK_1', 'ACE_CHEAP_CHUCK_2', 'ACE_CHEAP_CHUCK_3', 'ACE_CHEAP_CHUCK_4']
    #potentials = ['ACE_REP1', 'ACE_REP2', 'ACE_REP10', 'ACE_REP100', 'ACE_REP1000', 'ACE_REP10000', 'ACE_REP100000', 'ACE_REP1000000', 'ACE_CHUCK_3']
    potentials = ['CHIK', 'ACE_Chuck', 'ACE_CHEAP_CHUCK_2']
    potentials = ['CHIK', 'ACE_Chuck', 'ACE_REP_100_REF', 'ACE_REP_0.01_REF', 'ACE_KAMIL_5']
    specie = 'SiO'
    atoms.set_chemical_symbols(specie)
    rs = np.linspace(0.01, 6, 500)
    us_chik = [compute_CHIK(r, specie) for r in rs]
    results = {'CHIK': us_chik}
    fig, ax = plt.subplots()
    for potential in potentials:
        if potential == 'CHIK':
            continue
        config['potential'] = potential
        lammps = get_lammps_for_potential()
        atoms.calc = lammps
        MPI.COMM_WORLD.Barrier()
        results[potential] = [compute_u0(r) for r in rs]
        close_lammps(lammps)
    for potential, result in results.items():
        linestyle = ':' if potential == 'CHIK' else '-'
        if potential == 'ACE_KAMIL_5' or potential=='ACE_REP1' or potential == 'ACE_CHEAP_CHUCK_4':
            linestyle = '--'
        if potential == 'ACE_REP1':
            potential = 'ACE_REP0.01'
        if potential=='ACE_REP2':
            potential = 'ACE_REP0.5'
        plt.plot(rs, result, label=potential, linestyle=linestyle)
    plt.ylim((-20, 50))
    plt.legend()
    plt.ylabel("u / eV")
    plt.xlabel("r / Å")
    plt.title("SiO Dimer potentials comparison")
    plt.savefig("SiO_Dimer_comparison_ACE_new_REF_LINEAR.png", dpi=400)
    plt.ylim((1e-2, 1e8))
    plt.xlim((0, 1.2))
    ax.set_yscale('log')
    plt.savefig("SiO_Dimer_comparison_ACE_new_REF_LINEAR_LOG_CROPPED.png", dpi=400)

def compute_cristobalite_energy():
    atoms = file_conversion.read_reg('/mnt/scratch2/q13camb_scratch/adps2/work/cristabolite_alpha_444/relaxed_model/POSCAR')
    original_cell = atoms.get_cell()
    strains = np.linspace(0.9, 1.1, 50)
    config['potential'] = 'ACE_12_0.5_0.0003'
    lammps = get_lammps_for_potential()
    atoms.calc = lammps
    u = atoms.get_potential_energy()
    print("energy is", u)
    
    energies = []
    for scale_factor in strains:
        atoms.set_cell(original_cell * scale_factor, scale_atoms=True) # apply an affine distortion
        energies.append(atoms.get_potential_energy())
    plt.plot(strains * 1.65, energies)
    plt.savefig("cristobalite_energy_test2.png", dpi=400)
    
    close_lammps(lammps)
######

#offset = 15000
offset = 0
if __name__ == '__main__':
    test_all_potentials_dimer()
    #compute_cristobalite_energy()
    '''
    species = ['Si2', 'O2', 'SiO']
    for specie in species:
        atoms.set_chemical_symbols(specie)
        rs = np.linspace(0.01, 4, 50)
        us_bks = [compute_CHIK(r, specie) for r in rs]
        plt.plot(rs, us_bks, label = (specie+'(CHIK)'), linestyle = ':')
    plt.ylim((-100, 200))
    plt.legend()
    plt.ylabel("u / eV")
    plt.xlabel("r / Angstroms")
    plt.title(f"CHIK two-body test -- {config['potential']}")
    plt.savefig("twobody_chik.png", dpi=400)
    plt.clf()
    species = ['Si2', 'O2', 'SiO']
    for specie in species:
        atoms.set_chemical_symbols(specie)
        rs = np.linspace(0.01, 2, 50)
        us = [compute_u0(r) for r in rs]
        us_bks = [compute_CHIK(r, specie) for r in rs]
        plt.plot(rs, us, label = specie)
        #plt.plot(rs, us_bks, label = (specie+'(CHIK)'), linestyle = ':')
    plt.ylim((-100, 200))
    plt.legend()
    plt.ylabel("u / eV")
    plt.xlabel("r / Angstroms")
    plt.title(f"ACE two-body test -- {config['potential']}")
    plt.savefig("twobody_large_scale.png", dpi=400)
    plt.ylim((-20, 60))
    plt.savefig("twobody_small_scale.png", dpi=400)
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
    plt.title(f"ACE two-body test close-- {config['potential']}")
    plt.savefig("twobody_close.png", dpi=400)
    plt.clf()
    for specie in ['SiO']:
        atoms.set_chemical_symbols(specie)
        rs = np.linspace(0.01, 5, 50)
        us = [compute_u0(r) for r in rs]
        us_chik = [compute_CHIK(r, specie) for r in rs]
        plt.plot(rs, us, label = specie)
        plt.plot(rs, us_chik, label = (specie+'(CHIK)'), linestyle = ':')
    plt.ylim((-20, 50))
    plt.legend()
    plt.ylabel("u / eV")
    plt.xlabel("r / Angstroms")
    plt.title(f"ACE two-body test -- {config['potential']}")
    plt.savefig(f"twobody_SiO_CHIK_{config['potential']}_comp.png", dpi=400)
    plt.loglog()
    plt.ylim(top=None)
    plt.savefig(f"twobody_SiO_CHIK_{config['potential']}_comp_loglog.png", dpi=400)
    plt.clf()
    '''

MPI.COMM_WORLD.Barrier()



close_lammps(lammps)


MPI.Finalize()