import numpy as np
from ase.io import read, write
from ase.io.vasp import read_vasp
from ase import Atoms
from ase.constraints import FixSymmetry
from ase.filters import UnitCellFilter, ExpCellFilter, StrainFilter,FrechetCellFilter
from ase.optimize import BFGS, FIRE, MDMin, GPMin
from ase.spacegroup.symmetrize import check_symmetry
import makeconfig
from mpi4py import MPI
from ase.calculators.lammpslib import LAMMPSlib
import lammps as lmp
import file_conversion

convert_ase_to_bar=1.602176634e-19/1e-30/1e5



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

def run_joint_relax(atoms,calculator,fmax=1e-4,allow_tilt=False,Optimizer=BFGS,alpha=None):
    
    atoms.calc=calculator
    
    no_tilt_mask=[True,True,True,False,False,False]
    vc_mask=None
    if not allow_tilt:
        vc_mask=no_tilt_mask
    
    print("Initial Energy", atoms.get_potential_energy()," ev")
    print("Initial Stress",atoms.get_stress()*convert_ase_to_bar," bar")
 
    print("Initial symmetry at precision 1e-6")
    check_symmetry(atoms, 1.0e-6, verbose=True)
    check_symmetry(atoms, 1.0e-3, verbose=True)
 
    atoms.set_constraint(FixSymmetry(atoms))
 
    vc_filter=StrainFilter(atoms,mask=vc_mask)
    exp_filter=ExpCellFilter(atoms,mask=no_tilt_mask)
    if Optimizer.__name__ in ["BFGS",'SciPyFminCG']:
        dyn_cell = Optimizer(vc_filter,alpha=alpha)
        dyn_atoms_only = Optimizer(atoms,alpha=alpha)
        dyn_total=Optimizer(exp_filter,alpha=alpha)
    else:
        dyn_cell = Optimizer(vc_filter)
        dyn_atoms_only = Optimizer(atoms)
        dyn_total=Optimizer(exp_filter)
    
 
 
    # Run a optimisation for atomic positions
    # with every step rescaling the cell to minimise stress
    dyn_total.run(fmax=1e-3,steps=200)
    
    for _ in dyn_atoms_only.irun(fmax=fmax,steps=500):
        dyn_cell.run(fmax=fmax,steps=1000)
    dyn_atoms_only.run(fmax=fmax,steps=500)
 
    print("After keeping symmetry VC/FC relax Energy", atoms.get_potential_energy()," ev")
    print("After keeping symmetry VC/FC relax Stress",atoms.get_stress()*convert_ase_to_bar," bar")
 
    # We print out the initial symmetry groups at two different precision levels
    print("After keeping symmetry VC/FC relax symmetry at precision 1e-5")
    check_symmetry(atoms, 1.0e-5, verbose=True)
 
    # delete constrainsts and run a optimisation for atomic positions
    # with every step rescaling the cell to minimise stress
    atoms_symmetry=atoms.copy()
 
    del atoms.constraints
    for _ in dyn_atoms_only.irun(fmax=fmax,steps=200):
        dyn_cell.run(fmax=fmax,steps=500)
    dyn_atoms_only.run(fmax=fmax,steps=200)
 
    print("Final Energy", atoms.get_potential_energy()," ev")
    print("Final Stress",atoms.get_stress()*convert_ase_to_bar," bar")
 
    print("Final symmetry at precision 1e-4")
    check_symmetry(atoms, 1.0e-4, verbose=True)
    print("Final symmetry at precision 1e-5")
    check_symmetry(atoms, 1.0e-5, verbose=True)
 
    # compare symmetries
    from ase.utils.structure_comparator import SymmetryEquivalenceCheck
    comp = SymmetryEquivalenceCheck()
    if not comp.compare(atoms, atoms_symmetry):
        atoms=atoms_symmetry
        print("\nSYMMETRY IS NOT KEPT AFTER RELAXATION, USING SYMMETRYIC STRUCTURE\n")
 
    atoms_write=Atoms(symbols=atoms.symbols,positions=atoms.positions,cell=atoms.cell,pbc=atoms.pbc)
    spc = config['[1b_relax_sym]supercell_result']
    write("1b_relax_sym/POSCAR_relaxed_unitcell",atoms_write,format='vasp')
    atoms_write = atoms_write * spc
    write("1b_relax_sym/POSCAR_relaxed",atoms_write,format='vasp')
    return atoms

config = makeconfig.config()

if __name__ == '__main__':
    atoms = file_conversion.read_reg(config['input_struct'])

    run_joint_relax(atoms, calculator=get_lammps_for_potential())