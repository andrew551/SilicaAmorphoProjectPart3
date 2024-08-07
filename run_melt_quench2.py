import makeconfig
import sys, os
import datetime
from pathlib import Path

import create_melt_quench_ctrl
import makeconfig
import file_conversion
import json
nproc = 1

def prepare_output_folder(config):
    config['starttime'] = datetime.datetime.now().strftime("%Y%m%d%H%M%S")
    config['output_dir'] = config['output_dir_base'] / 'anneal_output' / config['starttime'] # note: we are using pathlib.Path object for the paths
    print('output_dir:', config['output_dir'])
    os.makedirs(config['output_dir'])

if __name__ == '__main__':
    config = makeconfig.config()
    nproc = config['NTASKS']
    if not config['OMP_NUM_THREADS'] == 1:
        raise Exception(f"expected OMP_NUM_THEADS = 1, got {config['OMP_NUM_THEADS']}")
    print('config data=', config, flush=True)
    #input_struct_path = Path('/mnt/scratch2/q13camb_scratch/adps2/input_folder2/models24k/Coords_3.dat')
    #input_struct_path = Path('/users/asmith/grun_in/model1536/POSCAR_1536')
    #input_struct_path = Path('/mnt/scratch2/q13camb_scratch/silica_plateau/chik5001/Coords_5001atoms_chik_min.dat')
    #prepare_output_folder(config)
    config['output_dir'] = Path('0_melt_quench').resolve()
    if config['need_anneal?']:
        input_struct_path = Path(config['input_struct'])
    else:
        raise Exception("config.json says anneal not needed!")
    
    regularised_input_path = config['output_dir'] / (input_struct_path.stem + '_regularised.dat')
    # convert the weird input format into normal lammps format
    file_conversion.convert_and_regularize_file(input_struct_path, regularised_input_path)
    # srun lammps anneal / relax
    command = create_melt_quench_ctrl.create_melt_quench_ctrl(regularised_input_path, config)
    print(command, flush=True)
    lammps_ctrl = str(config['output_dir'] / 'mq_silica.lammps_ctrl')
    with open(lammps_ctrl, 'w') as f:
        f.write(command)
    lammps_out = str(config['output_dir'] / 'mq_silicaL.out')
    anneal_command = f"srun -n {config['NTASKS']} {config['path_lammps']} -in {lammps_ctrl} > {lammps_out}"
    print(anneal_command, flush=True)
    meta_data = config.copy()
    meta_data['input_file'] = str(input_struct_path)
    meta_data['lammps.out'] = str(lammps_out)
    meta_data['jobtype'] = 'anneal'
    for k, v in meta_data.items():
        if isinstance(v, Path):
            meta_data[k] = str(v) # convert back to string for writing to txt
    with open(config['output_dir'] / 'metadata.txt', 'w', encoding='utf-8') as fp:           
        json.dump(meta_data, fp, indent=4)
    # run the lammps command
    os.system(anneal_command)
