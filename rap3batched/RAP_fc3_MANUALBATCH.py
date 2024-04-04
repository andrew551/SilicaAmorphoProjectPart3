import os, sys
sys.setdlopenflags(os.RTLD_NOW | os.RTLD_GLOBAL)
 
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
sys.path.insert(1, os.path.join(sys.path[0], '..')) # "hack" to add the parent directory to path
import makeconfig # this line needs to be after the sys.path.insert line
config = makeconfig.config()


me = MPI.COMM_WORLD.Get_rank()
nprocs = MPI.COMM_WORLD.Get_size()
split=MPI.COMM_WORLD.Split(me,key=0)



## All user-defined parameters
supercell_size = (1,1,1)
Force_cutoff = 9
FC3_cutoff = 3.2
displacement=0.01 #ShengBTE uses 0.01, while phono3py is using 0.03, phonopy is using 0.0212 
POSCAR_file='POSCAR'
lammps_log_file='test2.log'
potential='GAP'
NeighFile=POSCAR_file+str(FC3_cutoff)+str(Force_cutoff)+".NLC"



#write out type
WRITE_OUT_TYPE=0 
## 0 - file per process
## 1 - single shared file

##Checkpointing
TOTAL_CHECKPOINT=6
ME_CHECKPOINT=name = int(sys.argv[1])
## 0 - default - no checkpointing
# Depends on how much tasks per node, ChatGPT suggest 1 node and match the no_of task with no_process, in this case .lsf file is natoms ntasks

## Initialize cell and lammps calculators, please see https://wiki.fysik.dtu.dk/ase/ase/calculators/lammps.html
## Suggest have lammps python API installed first on the supercomputer, see README
atom_dict={'O':1, 'Si':2}
if potential=='ACE':
    cmds = makeconfig.get_potential_command(config).split('\n')
if potential=='GAP':
    raise Exception("unimplemented GAP!")
    path_GAP_potential = '/mnt/scratch2/users/axu/GAP'
    cmds = [	'pair_style quip',
	            'pair_coeff * * %s/sio2_potential_data/potential/silica_gap.xml "Potential xml_label=GAP_2021_4_19_120_7_32_55_336" 8 14'%(path_GAP_potential)]

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

"""
# get numbers to compute
if me < natoms%nprocs:
    ind_nums=natoms//nprocs+1
    ind_start=me*(natoms//nprocs+1)
else:
    ind_nums=natoms//nprocs
    ind_start=me*(natoms//nprocs)+natoms%nprocs
"""

folder="results_CP"+str(ME_CHECKPOINT)+"_"+str(TOTAL_CHECKPOINT)
output_file=folder+"/fc3_%d.hdf5"%((me+ME_CHECKPOINT*nprocs))


########
# Neighbourlist read in
########
time_step_neigh1=time()


if me == 0:
    with open(NeighFile, "rb") as f:
        NL_fc3=pickle.load(f)
        NL_force=pickle.load(f)
        compList=pickle.load(f)
else:
    NL_fc3=None
    NL_force=None
    compList=None

NL_fc3=MPI.COMM_WORLD.bcast(NL_fc3,root=0)
NL_force=MPI.COMM_WORLD.bcast(NL_force,root=0)
compList=MPI.COMM_WORLD.bcast(compList,root=0)

time_step_neigh2=time()
if me==0:
    print("Neighbour list is loaded in %d seconds\nForce constants are being calculated."%(time_step_neigh2-time_step_neigh1))



cL_length=compList.shape[0]

"""
if me < nprocs-cL_length%nprocs:
    ind_nums=cL_length//nprocs
    ind_start=me*(cL_length//nprocs)
else:
    ind_nums=cL_length//nprocs+1
    ind_start=me*(cL_length//nprocs)+me+cL_length%nprocs-nprocs"""

if (me+ME_CHECKPOINT*nprocs) < cL_length%(nprocs*TOTAL_CHECKPOINT):
    ind_nums=cL_length//(nprocs*TOTAL_CHECKPOINT)+1
    ind_start=(me+ME_CHECKPOINT*nprocs)*(cL_length//(nprocs*TOTAL_CHECKPOINT)+1)
else:
    ind_nums=cL_length//(nprocs*TOTAL_CHECKPOINT)
    ind_start=(me+ME_CHECKPOINT*nprocs)*(cL_length//(nprocs*TOTAL_CHECKPOINT))+cL_length%(nprocs*TOTAL_CHECKPOINT)

my_ind=np.arange(ind_start,ind_start+ind_nums)


print('Neighbour list completed')
displacement_vectors=np.array([[displacement,0,0], [0,displacement,0], [0,0,displacement]])

## Lammps regional commputation function
def fc3_all(i):
    print('Displace', i, 'out of', natoms, 'atoms.')
    # Print the stored terminal output. 
    stdout.flush()
    # Find the intersection of the regions of the first atom and second atom, only need to compute this region due to cutoff
    atoms_local=atoms.copy()
    # set calculators
    atoms_local.set_calculator(lammps)
    # set neighbourus of the first atom
    indice2= NL_fc3.get_neighbors(i)[0]
    if(min_cell_length<FC3_cutoff*2):
        indice2=np.unique(indice2)
    # find the index of the central atom
    self_index_indice2 = np.searchsorted(indice2, i)
    # To prevent redundant calculations, only displace it ones, but need to save self energy, so that ijk and jik can be computed (different in this new definition)
    indice_2_required=indice2[self_index_indice2:]
    
    # Number of displace then calculations
    ntasks=len(indice_2_required)
    forces_pp=np.zeros((3,natoms,3),dtype=np.float64)
    forces_pn=np.zeros((3,natoms,3),dtype=np.float64)
    forces_np=np.zeros((3,natoms,3),dtype=np.float64)
    forces_nn=np.zeros((3,natoms,3),dtype=np.float64)


    #initialise fc3 COO matrix and it's offset
    F_offset=0
    F=np.zeros((54*((len(indice2))**2),7))

    for m in range(ntasks):
        # Find the neighbour of i
        j=indice_2_required[m]
        intersect_indice_3=[]
        
        
        if(not j==i):
            # Find the neighbour of j
            indice3=NL_fc3.get_neighbors(j)[0]
            if(min_cell_length<FC3_cutoff*2):
                indice3=np.unique(indice3)
            # This will be the region where k has distance smaller or equal to 3.1 from i and j. 
            intersect_indice_3 = np.intersect1d(indice2,indice3)
            natoms3=len(intersect_indice_3)
        else:
            # For the same atom case, the intersection is basically the same region
            intersect_indice_3=indice2
            natoms3=len(intersect_indice_3)

        M=np.zeros((3,3,natoms3,3), dtype=np.float64)

        # Generate 3 displacements for each axis and calculate with lammps
        for x in range(3):
            atoms_local[i].position += displacement_vectors[x]
            for y in range(3):
                atoms_local[j].position += displacement_vectors[y]
                forces_pp[y]=atoms_local.get_forces()
                atoms_local[j].position -= 2*displacement_vectors[y]
                forces_pn[y]=atoms_local.get_forces()
                atoms_local[j].position += displacement_vectors[y]
                # Calculate anharmonic matrix for 1 atom iteration    
            # Return the displaced first atom to its original position,
            atoms_local[i].position -= 2*displacement_vectors[x]
            for y in range(3):
                atoms_local[j].position += displacement_vectors[y]
                forces_np[y]=atoms_local.get_forces()
                atoms_local[j].position -= 2*displacement_vectors[y]
                forces_nn[y]=atoms_local.get_forces()
                atoms_local[j].position += displacement_vectors[y]
                for k_idx in range(natoms3):
                    k=intersect_indice_3[k_idx]
                    M[x,y,k_idx]=(-forces_pp[y,k]+forces_pn[y,k]+forces_np[y,k]-forces_nn[y,k])/(4*displacement**2)
            atoms_local[i].position += displacement_vectors[x]
        
        xs,ys,ks,zs= np.meshgrid([0,1,2],[0,1,2],intersect_indice_3,[0,1,2])
        i_s=np.ones((3,3,natoms3,3))*i
        js=np.ones((3,3,natoms3,3))*j
        Temp=np.stack((i_s.flatten(),xs.flatten(),js.flatten(),ys.flatten(),ks.flatten(),zs.flatten(),M.flatten()),axis=1)
        #if not diagonal, symmetrize:
        if i!=j:
            Temp=np.vstack((Temp,np.stack((js.flatten(),ys.flatten(),i_s.flatten(),xs.flatten(),ks.flatten(),zs.flatten(),M.flatten()),axis=1)))

        len_temp=np.size(Temp,axis=0)
        F[F_offset:F_offset+len_temp,:]=Temp
        #increase offset
        F_offset=F_offset+len_temp
        del Temp


    #return values up to calculated
    return F[:F_offset,:]

def fc3_compList_all(i,j):
    print('Displace', i," and ",j,' atoms out of', cL_length, 'pairs.')
    stdout.flush()
    # First create the region

    atoms_local=atoms.copy()
    atoms_local.set_calculator(lammps)
    natoms_local=len(atoms_local)


    forces_pp=np.zeros((3,natoms_local,3),dtype=np.float64)
    forces_pn=np.zeros((3,natoms_local,3),dtype=np.float64)
    forces_np=np.zeros((3,natoms_local,3),dtype=np.float64)
    forces_nn=np.zeros((3,natoms_local,3),dtype=np.float64)
    
    indice2= NL_fc3.get_neighbors(i)[0]

    
    if(not j==i):
        indice3=NL_fc3.get_neighbors(j)[0]
        intersect_indice_3 = np.intersect1d(indice2,indice3)           
    else:
        intersect_indice_3=indice2
            
    natoms3=len(intersect_indice_3)     
        
    M=np.zeros((3,3,natoms3,3), dtype=np.float64)
    
    # Generate 3 displacements for each axis and calculate with lammps
    for x in range(3):
        atoms_local[i].position += displacement_vectors[x]
        for y in range(3):
            atoms_local[j].position += displacement_vectors[y]
            forces_pp[y]=atoms_local.get_forces()
            atoms_local[j].position -= 2*displacement_vectors[y]
            forces_pn[y]=atoms_local.get_forces()
            atoms_local[j].position += displacement_vectors[y]
            # Calculate anharmonic matrix for 1 atom iteration    
        # Return the displaced first atom to its original position,
        atoms_local[i].position -= 2*displacement_vectors[x]
        for y in range(3):
            atoms_local[j].position += displacement_vectors[y]
            forces_np[y]=atoms_local.get_forces()
            atoms_local[j].position -= 2*displacement_vectors[y]
            forces_nn[y]=atoms_local.get_forces()
            atoms_local[j].position += displacement_vectors[y]
            for k_idx in range(natoms3):
                k_force_local=intersect_indice_3[k_idx]
                M[x,y,k_idx]=(-forces_pp[y, k_force_local]+forces_pn[y, k_force_local]+forces_np[y, k_force_local]-forces_nn[y, k_force_local])/(4*displacement**2)
        atoms_local[i].position += displacement_vectors[x]
            
        xs,ys,ks,zs= np.meshgrid([0,1,2],[0,1,2],intersect_indice_3,[0,1,2])
        i_s=np.ones((3,3,natoms3,3))*i
        js=np.ones((3,3,natoms3,3))*j
        F=np.stack((i_s.flatten(),xs.flatten(),js.flatten(),ys.flatten(),ks.flatten(),zs.flatten(),M.flatten()),axis=1)
        #if not diagonal, symmetrize:
        if i!=j:
            F=np.vstack((F,np.stack((js.flatten(),ys.flatten(),i_s.flatten(),xs.flatten(),ks.flatten(),zs.flatten(),M.flatten()),axis=1)))


    #return values up to calculated
    return F



def fc3_region(i,j):
    print('Displace', i," and ",j,' atoms out of', cL_length, 'pairs.')
    stdout.flush()
    # First create the region
    force_indices= np.intersect1d(NL_force.get_neighbors(i)[0],NL_force.get_neighbors(j)[0])

    i_force_local = np.searchsorted(force_indices, i)
    j_force_local = np.searchsorted(force_indices, j)
    atoms_local=atoms[force_indices].copy()
    atoms_local.set_calculator(lammps)
    natoms_local=len(atoms_local)


    forces_pp=np.zeros((3,natoms_local,3),dtype=np.float64)
    forces_pn=np.zeros((3,natoms_local,3),dtype=np.float64)
    forces_np=np.zeros((3,natoms_local,3),dtype=np.float64)
    forces_nn=np.zeros((3,natoms_local,3),dtype=np.float64)
    
    indice2= NL_fc3.get_neighbors(i)[0]

    
    if(not j==i):
        indice3=NL_fc3.get_neighbors(j)[0]
        intersect_indice_3 = np.intersect1d(indice2,indice3)           
    else:
        intersect_indice_3=indice2
            
    natoms3=len(intersect_indice_3) 

    local_forces_k_indices=np.searchsorted(force_indices, intersect_indice_3)       
        
    M=np.zeros((3,3,natoms3,3), dtype=np.float64)
    
    # Generate 3 displacements for each axis and calculate with lammps
    for x in range(3):
        atoms_local[i_force_local].position += displacement_vectors[x]
        for y in range(3):
            atoms_local[j_force_local].position += displacement_vectors[y]
            forces_pp[y]=atoms_local.get_forces()
            atoms_local[j_force_local].position -= 2*displacement_vectors[y]
            forces_pn[y]=atoms_local.get_forces()
            atoms_local[j_force_local].position += displacement_vectors[y]
            # Calculate anharmonic matrix for 1 atom iteration    
        # Return the displaced first atom to its original position,
        atoms_local[i_force_local].position -= 2*displacement_vectors[x]
        for y in range(3):
            atoms_local[j_force_local].position += displacement_vectors[y]
            forces_np[y]=atoms_local.get_forces()
            atoms_local[j_force_local].position -= 2*displacement_vectors[y]
            forces_nn[y]=atoms_local.get_forces()
            atoms_local[j_force_local].position += displacement_vectors[y]
            for k_idx in range(natoms3):
                k_force_local=local_forces_k_indices[k_idx]
                M[x,y,k_idx]=(-forces_pp[y, k_force_local]+forces_pn[y, k_force_local]+forces_np[y, k_force_local]-forces_nn[y, k_force_local])/(4*displacement**2)
        atoms_local[i_force_local].position += displacement_vectors[x]
            
        xs,ys,ks,zs= np.meshgrid([0,1,2],[0,1,2],intersect_indice_3,[0,1,2])
        i_s=np.ones((3,3,natoms3,3))*i
        js=np.ones((3,3,natoms3,3))*j
        F=np.stack((i_s.flatten(),xs.flatten(),js.flatten(),ys.flatten(),ks.flatten(),zs.flatten(),M.flatten()),axis=1)
        #if not diagonal, symmetrize:
        if i!=j:
            F=np.vstack((F,np.stack((js.flatten(),ys.flatten(),i_s.flatten(),xs.flatten(),ks.flatten(),zs.flatten(),M.flatten()),axis=1)))


    #return values up to calculated
    return F


def write_fpp(output_file,fc3,my_ind,me,me_cp=0):
    with h5py.File(output_file, "w") as f:
        # Create a dataset with chunking, compression, and filters
        
        
        dataset = f.create_dataset(
            "fc3",
            data=fc3,
            compression="gzip",           # Use gzip compression
            compression_opts=4            # Set compression level
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
            data=me_cp
        )

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

def print_once(args):
    if me==0:
        print(args)
        stdout.flush()

##----------------------
##      MAIN CODE
##----------------------


if __name__ == "__main__":

    if TOTAL_CHECKPOINT!=0:
        cp_start=ME_CHECKPOINT*ind_nums//(TOTAL_CHECKPOINT)+ind_start
        cp_num=ind_nums//(TOTAL_CHECKPOINT)
        #my_ind=np.arange(cp_start,cp_start+cp_nums)
        #print_once("Using %d checkpoints"%(ME_CHECKPOINT))
    else:
        ind_cp=-1

    #Create output folder
    if me==0:
        makedirs(folder,exist_ok=True)
        #system("rm ./"+folder+"/*")

    MPI.COMM_WORLD.Barrier()

    fc3=[]

    if (min_cell_length<(Force_cutoff+FC3_cutoff)*2.1):
        print_once('Whole unit cell is being calculated.')
        # Initialize lists to store the results
        
        for i in my_ind:
            res=fc3_compList_all(int(compList[i,0]),int(compList[i,1]))
            fc3.append(res)

            if ind_cp==i:
                write_fpp(output_file,fc3,my_ind,me,ind_cp)
                print('write_output, %s'%i)
                if cp_num<CHECKPOINT:
                    cp_num+=1
                restart_lammps(lammps)
            stdout.flush()
        close_lammps(lammps)
    
    #Regional calculation
    else:
        print_once("Regional acceleration protocol is active.\nLoading neighbour list from %s:"%(NeighFile))



        for i in my_ind:
            res=fc3_region(int(compList[i,0]),int(compList[i,1]))
            fc3.append(res)

            fc3_checkpoint=np.vstack(fc3)
            write_fpp(output_file,fc3_checkpoint,my_ind,me,ME_CHECKPOINT)
            #lammps.clean()

            stdout.flush()
        close_lammps(lammps)
    
    fc3=np.vstack(fc3)

    stdout.flush()


    step1_time=time()
    print_once("Calculation finished in %d seconds.\nStarting writeout."%(step1_time-start_time))

    if WRITE_OUT_TYPE==0:
        write_fpp(output_file,fc3,my_ind,me,ME_CHECKPOINT)
        
    if me==0:    
        end_time=time()
        print("Total runtime: %d seconds"%(end_time-start_time))

    stdout.flush()

MPI.COMM_WORLD.Barrier()
MPI.Finalize()