# RAP3 - MPI version

This is the MPI version of RAP 3, which allows better scalability without a master process. The pairs of atoms to displacce are listed before computation to get approximately the same runtime on all the processes.

## To use
- Modify the scripts to use your filenames, potentials and parameters.
- Run `calc_NL_compList.py` to calculate the NeighbourLists and the ComputeList.
- Run `RAP_fc3_compList.py` with srun/mpirun
- Collect outputs with `collect_RAP3.py`, which also includes acoustice sum rule (in very near future version).
- Output is in `fc3.hdf5`

**If you have any questions or bugs, let me (Balazs Pota, bp443@cam.ac.uk) know.**

### Dependecies:

`scipy`, `lammps`, `ase` and `h5py` python packages are used by the main code.
