#!/bin/bash
#SBATCH --time=0-00:30:00
#SBATCH --job-name=compile
#####SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu 5G
####SBATCH --cpus-per-task=1
#SBATCH --partition=k2-epsrc-himem

#SBATCH --output %j.out
#SBATCH --error  %j.err

module purge
module load mpi/openmpi/4.1.1/gcc-9.3.0

echo "start compile"
gfortran -O2 -ffree-line-length-0 internal_strain.f90 -o int_strain.exe
mpif90 -O2 -ffree-line-length-0 Gruneisen_mpi_opt.f90 -o grun_mpi.exe
gfortran -O2 -ffree-line-length-0 postproc.f90 -o postproc.exe
gfortran -O2 -ffree-line-length-0 convert_d3.f90 -o convert_d3.exe

echo "end of compilation"