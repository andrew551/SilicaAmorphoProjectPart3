Pre-requisites:
===============
lammps, with QUIP (for GAP)
python library requirements: ase, numypy, scipy, lammps
virtual enviroment

Installation
============
git clone https://github.com/andrew551/SilicaAmorphoProjectPart3.git

Compilation of fortran binaries:
================================
cd  path/to/SilicaAmorphoProjectPart3/gruneisen
adjust compile.sh (if needed)
sbatch compile.sh
move generated binaries to  path/to/SilicaAmorphoProjectPart3/bin/gruneisen

Running workflow
================

1) make working directory for structure analysis (let its name be "wrkdir")

2) Copy bash scripts:
    cd workdir
    cp path/to/SilicaAmorphoProjectPart3/demo_scripts_kelvin/* .

3) update files AMORPHO_PATH, LAMMPS_PATH, VENV_PATH
4) configure config.json
    ** input_struct_path: absolute path to vasp, xyz, lammps (or other) file containing atomic positions (the code automatically determines the input type -- if there is a bug for a certain type of input contact me)
    ** material: currently Si and SiO2 supported
    ** potentials: currently ACE_Chuck, ACE_Deringher, GAP_Si, GAP_SiO2 are supported -- to add a new type of potential the function makeconfig.get_potential_command() may be updated
    "need_anneal?":false,               ##is an anneal needed before relax?
    "need_relax?":true,                 ##is relaxation of the input model required?
    "[1_relax]keep_cuboidal":true,      ##constrain to cuboidal cell relaxation
    "[1_relax]thresholds": [1e-9,1e-10,1e-11,1e-12,1e-12,1e-12,1e-13,1e-13,1e-13,1e-14,1e-14,1e-14,1e-15,1e-15, 8e-16, 5e-16, 3e-16, 1e-16], ## thresholds for relax cycles
    "[fc2]_FC2_cutoff":12,              ## force constant cutoffs (in Angstrom)
    "[fc2]_Force_cutoff":12,
    "[4_fc3]fc3_cutoff": 5.0,
    "[4_fc3]force_cutoff": 12.0,
    "[4_fc3]displacement": 0.01        ## fc3 atom displacement (Angstroms)
5) update partition, ntasks, cpus_per_task, etc. as required

6) launch scripts:
(sbatch submit_anneal) - if anneal is required
sbatch submit_relax

sbatch submit_fc2

sbatch submit_NL3
sbatch submit_workfc3

sbatch submit_gruneisen_2
sbatch submit_gruneisen_postproc



