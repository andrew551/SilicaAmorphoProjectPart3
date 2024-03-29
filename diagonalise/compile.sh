module purge
module load GCC/11.3.0  OpenMPI/4.1.4
module load ELPA/2022.11.001 ScaLAPACK/2.2.0-fb

EPATH=/sulis/easybuild/software/ELPA/2022.11.001-foss-2022a

mpif90  -O3 test_real2.F90  -I${EPATH}/include/elpa_openmp-2022.11.001/modules \
-L${EPATH}/libs -lelpa_openmp  -lscalapack -o elpa_eigen.exe









