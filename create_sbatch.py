import os
import makeconfig

config = makeconfig.config()

# load correct modules for Kelvin
def get_load_modules():
	cmd = 'module purge;\n\
module load compilers/gcc/9.3.0;\n\
module load mpi/openmpi/4.1.1/gcc-9.3.0;\n\
module load libs/atlas/3.10.3/gcc-9.3.0;\n\
module load apps/cmake/3.25.1/gcc-9.3.0;\n\
module load apps/python3/3.10.5/gcc-9.3.0;\n'
	return cmd

def get_lammps_srun(filename_ctrl, output_file, num_MPI_proc=None):
	if num_MPI_proc is None:
		return 'srun %s -in %s > %s\n'%(config['path_lammps'], filename_ctrl, output_file)
	return 'srun -n %d %s -in %s > %s\n'%(num_MPI_proc, config['path_lammps'], filename_ctrl, output_file)

def run_lammps_direct(filename_ctrl, output_file, num_MPI_proc=None):
	cmd = get_load_modules()
	cmd += get_lammps_srun(filename_ctrl, output_file)
	print(cmd)
	os.system(cmd)


def run_lammps_sbatch(files_in, files_out, launch_name, ntasks, cpus_per_task, mem_per_cpu=5):
	text =f'#!/bin/bash --login\n\
#SBATCH --nodes=1\n\
#SBATCH --ntasks={ntasks}\n\
#SBATCH --time=1:00:00\n\
#SBATCH --mail-type=NONE\n\
#SBATCH --cpus-per-task={cpus_per_task}\n\
#SBATCH --mem-per-cpu={mem_per_cpu}G\n\
#SBATCH --partition=k2-epsrc\n'

	text += get_load_modules()


	# TODO: understand this line 2
	text += f'source {config["path_venv"]} \n\
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK \n\
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:{config["path_lammps"]}\n\
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:{config["path_lammps_src"]}\n'

	for fin, fout in zip(files_in, files_out):
		text += get_lammps_srun(fin, fout, num_MPI_proc=cpus_per_task)
	with open(launch_name,'w') as tmpf:
		tmpf.write(text)
	cmd=f'sbatch {launch_name}'
	os.system(cmd)
	print (cmd)