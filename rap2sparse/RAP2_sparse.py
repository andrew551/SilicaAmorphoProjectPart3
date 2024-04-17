from mpi4py import MPI
from ase.atoms import Atoms
from ase.neighborlist import NeighborList, NewPrimitiveNeighborList
from ase.calculators.lammpslib import LAMMPSlib
from ase.io.vasp import read_vasp
import numpy as np
from time import time
from os import makedirs,system
from sys import stdout
import h5py 
import pickle
import lammps as lmp
import os
import sys
sys.path.insert(1, os.path.join(sys.path[0], '..')) # "hack" to add the parent directory to path
import makeconfig # this line needs to be after the sys.path.insert line
config = makeconfig.config()
me = MPI.COMM_WORLD.Get_rank()
nprocs = MPI.COMM_WORLD.Get_size()
split=MPI.COMM_WORLD.Split(me,key=0)



FC2_cutoff = config['[fc2]_FC2_cutoff']
Force_cutoff = config['[fc2]_Force_cutoff']
displacement=0.0212
POSCAR_file='POSCAR'
lammps_log_file=None#"test.log"
#potential='ACE'
chunks=None #(24, 24, 3, 3)
NeighFile=POSCAR_file+str(FC2_cutoff)+str(Force_cutoff)+".NL"
folder="Results"
output_file=folder+"/fc2_%d.hdf5"%(me)


#write out type
WRITE_OUT_TYPE=0 
## 0 - file per process
## 1 - single shared file

##Checkpointing
CHECKPOINT=0
## 0 - default - no checkpointing


## Initialize cell and lammps calculators, please see https://wiki.fysik.dtu.dk/ase/ase/calculators/lammps.html
## Suggest have lammps python API installed first on the supercomputer, see README
if config['material'] == 'SiO2':
    atom_dict={'O':1, 'Si':2}
elif config['material'] == 'Si':
    atom_dict={'Si':1}
else:
    raise Exception(f"unsupported material: {config['material']}")

cmds = makeconfig.get_potential_command(config).split('\n')

'''
if potential=='ACE':
    cmds = makeconfig.get_potential_command(config).split('\n')
    path_ACE_potential = config['path_ACE_potential']
    cmds = [    'pair_style hybrid/overlay pace table spline 6000',
                'pair_coeff * * pace %s/SiO2-4_24-20-16-12.yace O Si'%path_ACE_potential,
                'pair_coeff 1 1 table %s/SiO2-4_24-20-16-12_pairpot.table O_O'%path_ACE_potential,
                'pair_coeff 1 2 table %s/SiO2-4_24-20-16-12_pairpot.table O_Si'%path_ACE_potential,
                'pair_coeff 2 2 table %s/SiO2-4_24-20-16-12_pairpot.table Si_Si'%path_ACE_potential
    ]    
if potential=='GAP':
    raise Exception("unimplemented GAP!")
    path_GAP_potential = '/work/e89/e89/bp443/GAP'
    cmds = [	'pair_style quip',
	            'pair_coeff * * %s/sio2_potential_data/potential/silica_gap.xml "Potential xml_label=GAP_2021_4_19_120_7_32_55_336" 8 14'%(path_GAP_potential)]
'''

if me == 0:
    print("Initialising LAMMPS.")
    stdout.flush()


lmp_test=lmp.lammps()
if lmp_test.has_mpi_support:
    comm_lammps=split
else:
    comm_lammps=None
lmp_test.close()
del lmp_test

lammps = LAMMPSlib(lmpcmds=cmds,log_file=lammps_log_file,atom_types=atom_dict, keep_alive='TRUE',comm=comm_lammps)


# generate neighbor list
start_time = time()
atoms = read_vasp(POSCAR_file)
natoms=len(atoms)
cell_parameters=atoms.cell.cellpar()
min_cell_length=min(cell_parameters[:3])

# get numbers to compute
if me < natoms%nprocs:
    ind_nums=natoms//nprocs+1
    ind_start=me*(natoms//nprocs+1)
else:
    ind_nums=natoms//nprocs
    ind_start=me*(natoms//nprocs)+natoms%nprocs

my_ind=np.arange(ind_start,ind_start+ind_nums)





#validate chunks
if not chunks is None:
    if chunks[0]>ind_nums:
        chunks=(ind_nums,chunks[1],3,3)


## Lammps all commputation function
def f1(i):
    print(i, 'out of', len(atoms), 'environments')
    stdout.flush()
    # neighbour list is sorted, so have to find its indice for displacement
    self_index = i
    atoms_local=atoms.copy()
    atoms_local.set_calculator(lammps)

    # forcep is force on positive displacement, forcen is force negative displacement
    forcep=np.zeros((3,natoms,3), dtype=np.float64)
    forcen=np.zeros((3,natoms,3), dtype=np.float64)
    dynmat=np.zeros((3,natoms,3), dtype=np.float64)
    
    # Generate 6 displacements for each axis and calculate with lammps
    atoms_local[self_index].position += [displacement,0,0]
    forcep[0]=atoms_local.get_forces()
    atoms_local[self_index].position += [-2*displacement,0,0]
    forcen[0]=atoms_local.get_forces()

    atoms_local[self_index].position += [displacement,displacement,0]
    forcep[1]=atoms_local.get_forces()
    atoms_local[self_index].position += [0,-2*displacement,0]
    forcen[1]=atoms_local.get_forces()
    
    atoms_local[self_index].position += [0,displacement,displacement]
    forcep[2]=atoms_local.get_forces()
    atoms_local[self_index].position += [0,0,-2*displacement]
    forcen[2]=atoms_local.get_forces()

    # Calculate dynamic matrix, crude formula
    fc=(forcen-forcep)/(2*displacement)

    return np.transpose(fc,(1,0,2))

## Lammps regional commputation function
def f2_dense(i):
    print(i, 'out of', len(atoms), 'environments')
    stdout.flush()
    indices= NL_force.get_neighbors(i)[0]
    ind_fc2=NL_fc2.get_neighbors(i)[0]
    n=len(indices)
    # neighbour list is sorted, so have to find its indice for displacement
    self_index = np.searchsorted(indices,i)
    ind_fc2_internal=np.searchsorted(indices,ind_fc2)
    atoms_local=atoms[indices].copy()
    atoms_local.set_calculator(lammps)

    # forcep is force on positive displacement, forcen is force negative displacement
    forcep=np.zeros((3,n,3), dtype=np.float64)
    forcen=np.zeros((3,n,3), dtype=np.float64)

    # Generate 6 displacements for each axis and calculate with lammps
    atoms_local[self_index].position += [displacement,0,0]
    forcep[0]=atoms_local.get_forces()
    atoms_local[self_index].position += [-2*displacement,0,0]
    forcen[0]=atoms_local.get_forces()

    atoms_local[self_index].position += [displacement,displacement,0]
    forcep[1]=atoms_local.get_forces()
    atoms_local[self_index].position += [0,-2*displacement,0]
    forcen[1]=atoms_local.get_forces()

    atoms_local[self_index].position += [0,displacement,displacement]
    forcep[2]=atoms_local.get_forces()
    atoms_local[self_index].position += [0,0,-2*displacement]
    forcen[2]=atoms_local.get_forces()

    # Calculate fc2 matrix, crude formula
    
    fc=(forcen-forcep)/(2*displacement)
    fc=fc[:,ind_fc2_internal,:]


    return np.transpose(fc,(1,0,2)), ind_fc2

def f2_sparse(i):
    print(i, 'out of', len(atoms), 'environments')
    stdout.flush()
    indices= NL_force.get_neighbors(i)[0]
    ind_fc2=NL_fc2.get_neighbors(i)[0]
    n=len(indices)
    # neighbour list is sorted, so have to find its indice for displacement
    self_index = np.searchsorted(indices,i)
    ind_fc2_internal=np.searchsorted(indices,ind_fc2)
    atoms_local=atoms[indices].copy()
    atoms_local.set_calculator(lammps)

    # forcep is force on positive displacement, forcen is force negative displacement
    forcep=np.zeros((3,n,3), dtype=np.float64)
    forcen=np.zeros((3,n,3), dtype=np.float64)

    # Generate 6 displacements for each axis and calculate with lammps
    atoms_local[self_index].position += [displacement,0,0]
    forcep[0]=atoms_local.get_forces()
    atoms_local[self_index].position += [-2*displacement,0,0]
    forcen[0]=atoms_local.get_forces()

    atoms_local[self_index].position += [displacement,displacement,0]
    forcep[1]=atoms_local.get_forces()
    atoms_local[self_index].position += [0,-2*displacement,0]
    forcen[1]=atoms_local.get_forces()

    atoms_local[self_index].position += [0,displacement,displacement]
    forcep[2]=atoms_local.get_forces()
    atoms_local[self_index].position += [0,0,-2*displacement]
    forcen[2]=atoms_local.get_forces()

    # Calculate fc2 matrix, crude formula
    
    fc=(forcen-forcep)/(2*displacement)

    js,xs,ys= np.meshgrid(ind_fc2,[0,1,2],[0,1,2],indexing="ij")
    i_s=np.ones((len(ind_fc2),3,3))*i
    fc=fc[:,ind_fc2_internal,:].transpose((1,0,2))

    res=np.stack((i_s.flatten(),js.flatten(),xs.flatten(),ys.flatten(),fc.flatten()),axis=1)


    return res

def write_fpp(output_file,fc2,chunks,my_ind,me,ind_cp=0):
    with h5py.File(output_file, "w") as f:
        # Create a dataset with chunking, compression, and filters
        
        if chunks is None:
            dataset = f.create_dataset(
                "fc2",
                data=fc2,
                #compression="gzip",           # Use gzip compression
                #compression_opts=4            # Set compression level
            )
        else:
            dataset = f.create_dataset(
                "fc2",
                data=fc2,
                chunks=chunks  # Define chunk size
                #compression="gzip",           # Use gzip compression
                #compression_opts=4            # Set compression level
            )
        dataset = f.create_dataset(
            "indices",
            data=my_ind
        )
        dataset = f.create_dataset(
            "rank",
            data=me
        )
        dataset = f.create_dataset(
            "checkpoint",
            data=ind_cp
        )

def restart_lammps(LAMMPSLibObject):
    if LAMMPSLibObject.started:
        LAMMPSLibObject.lmp.command("clear")
    # hope there's no other state to be reset
    LAMMPSLibObject.started = True
    LAMMPSLibObject.initialized = False
    LAMMPSLibObject.previous_atoms_numbers = []

def close_lammps(lammps):
    if lammps.lmp is not None:
        lammps.lmp.close()

def print_once(args):
    if me==0:
        print(args)
        stdout.flush()


##----------------------
##      MAIN CODE
##----------------------

if __name__ == "__main__":

    if CHECKPOINT!=0:
        ind_cp=ind_nums//(CHECKPOINT+1)+ind_start
        cp_num=0
        print_once("Using %d checkpoints"%(CHECKPOINT))
    else:
        ind_cp=-1

    #Create output folder
    if me==0:
        makedirs(folder,exist_ok=True)
        system("rm ./"+folder+"/*")

    MPI.COMM_WORLD.Barrier()

    #Whole unit cell calculation
    if (min_cell_length<Force_cutoff*2.1):
        ## Define dynmat
        fc2=np.zeros((ind_nums,natoms,3,3))
        print_once('Whole unit cell is being calculated.')
        for i in my_ind:
            fc=f1(i)
            fc2[i-ind_start,:,:,:]=fc
            if ind_cp==i:
                write_fpp(output_file,fc2,chunks,my_ind,me,ind_cp)
                if cp_num<CHECKPOINT:
                    cp_num+=1
                restart_lammps(lammps)
            stdout.flush()
        close_lammps(lammps)
    
    #Regional calculation
    else:
        print_once("Regional acceleration protocol is active.\nLoading neighbour list from %s:"%(NeighFile))
        time_step_neigh1=time()
        if me == 0:
            with open(NeighFile, "rb") as f:
                NL_fc2=pickle.load(f)
                NL_force=pickle.load(f)
        else:
            NL_fc2=None
            NL_force=None
        NL_fc2=MPI.COMM_WORLD.bcast(NL_fc2,root=0)
        NL_force=MPI.COMM_WORLD.bcast(NL_force,root=0)


        time_step_neigh2=time()
        print_once("Neighbour list is loaded in %d seconds\nForce constants are being calculated."%(time_step_neigh2-time_step_neigh1))

        MPI.COMM_WORLD.Barrier()

        fc2=[]

        for i in my_ind:
            fc2.append(f2_sparse(i))
            
            if ind_cp==i:
                fc2_np=np.vstack(fc2)
                write_fpp(output_file,fc2_np,chunks,my_ind,me,ind_cp)
                if cp_num<CHECKPOINT:
                    cp_num+=1
                ind_cp=cp_num*ind_nums//(CHECKPOINT+1)+ind_start
                lammps.clean()
            stdout.flush()

        close_lammps(lammps)

        fc2=np.vstack(fc2)


    MPI.COMM_WORLD.Barrier()
    stdout.flush()


    step1_time=time()
    print_once("Calculation finished in %d seconds.\nStarting writeout."%(step1_time-start_time))



    if WRITE_OUT_TYPE==0:
        write_fpp(output_file,fc2,chunks,my_ind,me)
        
    if me==0:    
        end_time=time()
        print("Total runtime: %d seconds"%(end_time-start_time))



MPI.COMM_WORLD.Barrier()
MPI.Finalize()
