import makeconfig
import sys, os
import datetime
from pathlib import Path

import create_anneal_ctrl
import makeconfig
import create_sbatch

nproc = 1

if __name__ == '__main__':
    config = makeconfig.config()
    input_struct_path = '/users/asmith/grun_in/models24k/Coords_fixed'
    config['starttime'] = datetime.datetime.now().strftime("%Y%m%d%H%M%S")
    config['output_dir'] = config['output_dir_base'] / 'anneal_output' / config['starttime'] # the initial slash cannot be skipped!
    print('output_dir:', config['output_dir'])
    os.makedirs(config['output_dir'])
    command = create_anneal_ctrl.create_anneal_ctrl(config, input_struct_path)
    print(command)
    lammps_ctrl = str(config['output_dir'] / 'anneal1500K_silica.lammps_ctrl')
    with open(lammps_ctrl, 'w') as f:
        f.write(command)
    lammps_out = str(config['output_dir'] / 'anneal1500K_silicaL.out')
    #create_sbatch.run_lammps_sbatch([lammps_ctrl], [lammps_out], 'temp', nproc, 1)
    anneal_command = create_sbatch.get_lammps_srun(lammps_ctrl, lammps_out, nproc)
    print(anneal_command)
    os.system(anneal_command)
    # .//users/asmith/programs/lammps/build/lmp -in /users/asmith/grun_out/anneal_output/20240306190529/anneal1500K_silica.lammps_ctrl