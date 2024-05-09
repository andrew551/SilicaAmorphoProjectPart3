import makeconfig
import sys, os
import datetime
from pathlib import Path

import create_relax_ctrl
import makeconfig
import file_conversion
import json
import datetime

nproc = 1

def prepare_output_folder(config):
    config['starttime'] = datetime.datetime.now().strftime("%Y%m%d%H%M%S")
    config['output_dir'] = config['output_dir_base'] / 'relax_output' / config['starttime'] # note: we are using pathlib.Path object for the paths
    print('output_dir:', config['output_dir'])
    os.makedirs(config['output_dir'] / 'steps')

if __name__ == '__main__':
    t_start = datetime.datetime.now()
    print('start-time:', t_start)
    config = makeconfig.config()
    nproc = config['NTASKS']
    if not config['OMP_NUM_THREADS'] == 1:
        raise Exception(f"expected OMP_NUM_THEADS = 1, got {config['OMP_NUM_THEADS']}")
    config['output_dir'] = Path('1_relax').resolve()
    if config['need_anneal?']:
        input_struct_path = Path('0_anneal/Coords_ACE_cg_min.dat').resolve()
    else:
        input_struct_path = config['input_struct']
    print('config data=', config, flush=True)
    #input_struct_path = Path('/mnt/scratch2/q13camb_scratch/adps2/output_folder1/anneal_output/20240308155312/Coords_ACE_cg_min.dat') # 24k
    #input_struct_path = Path('/mnt/scratch2/q13camb_scratch/silica_plateau/chik5001/Coords_5001atoms_chik_min.dat')
    #input_struct_path  = Path('/mnt/scratch2/q13camb_scratch/silica_plateau/chik5001_q10/Coords_5001atoms_chik_min1_4000q10.dat')
    #input_struct_path = Path('/mnt/scratch2/q13camb_scratch/adps2/output_folder_chik5001/md_output/20240321184831/2_md300.dat')
    #input_struct_path = Path('/mnt/scratch2/q13camb_scratch/adps2/output_folder1/anneal_output/20240308185606/Coords_ACE_cg_min.dat') #1536
    #input_struct_path = Path('/mnt/scratch2/q13camb_scratch/silica_plateau/chik5001_q10_5200K/Coords_5001atoms_chik_min_5200q10.dat')
    #input_struct_path  =Path('/mnt/scratch2/q13camb_scratch/adps2/input_folder2/deringher5184/POSCAR_5184 1')
    #input_struct_path = Path('/mnt/scratch2/q13camb_scratch/adps2/quenched6000K/initial_model/Coords_5001atoms_chik_min_q11_6000K.dat')
    #prepare_output_folder(config)

    # iteratively relax the structure
    #thresholds = [1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8,1e-9,1e-10,1e-11,1e-12,1e-12,1e-12,1e-13,1e-13,1e-13,1e-14,1e-14,1e-14,1e-15,1e-15, 8e-16, 5e-16, 3e-16, 1e-16]
    #thresholds = [1e-3, 1e-4, 1e-6, 1e-8, 1e-10, 1e-11,1e-12,1e-13,1e-14,1e-14,1e-14,1e-14,1e-14,1e-14,1e-14,1e-14,1e-14,1e-14,1e-15,1e-15]
    #thresholds = [1e-10, 1e-11,1e-12,1e-13,1e-14,1e-14,1e-14,1e-14,1e-14,1e-14,1e-14,1e-14,1e-14,1e-14,1e-15,1e-15]
    #thresholds = [1e-9,1e-10,1e-11,1e-12,1e-12,1e-12,1e-13,1e-13,1e-13,1e-14,1e-14,1e-14,1e-15,1e-15, 8e-16, 5e-16, 3e-16, 1e-16]
    thresholds = config['[1_relax]thresholds']
    number_relaxations = len(thresholds)*2 + 1

    if number_relaxations % 2 == 0:
        raise Exception("number of relaxations must be odd (to get non-vc constrained last relax)")
    if (number_relaxations-2)//2 >= len(thresholds):
        raise Exception(f"not enough thresholds provided ({len(thresholds)}) for given number_relaxations ({number_relaxations})")
    meta_data = config.copy()
    meta_data['input_file'] = str(input_struct_path)
    meta_data['initial_density'] = file_conversion.get_density(file_conversion.read_reg(input_struct_path))
    meta_data['jobtype'] = 'relax'
    meta_data['keep_cuboidal'] = config['[1_relax]keep_cuboidal']
    #meta_data['thresholds'] = thresholds
    meta_data['number_relaxations'] = number_relaxations
    for k, v in meta_data.items():
        if isinstance(v, Path):
            meta_data[k] = str(v) # convert back to string for writing to txt
    with open(config['output_dir'] / 'metadata.txt', 'w', encoding='utf-8') as fp:           
        json.dump(meta_data, fp, indent=4)

    #####################################################
    #  relaxation procedure
    #######################################################

    # standardize POSCAR file and write to POSCAR_SYMMETRISED # skip this step (for glasses this is definately ok as they don't have symmetry)
    #standardize_vasp_and_write_POSCAR(input_vasp_file)

    # convert to LAMMPS
    n_relax = 1
    number_structure = 1

    lammps_structure_file = config['output_dir'] / 'steps' / ('struct_%d_relax_%d.lammps_struct' % (number_structure, n_relax))
    #convert_vasp_to_lammps('POSCAR_SYMMETRISED',lammps_structure_file)

    initial_struct_atoms = file_conversion.read_reg(input_struct_path)
    '''
    if '[1_relax]supercell' in config:
        supercell = config['[1_relax]supercell']
        print('using supercell', supercell)
        initial_struct_atoms *= tuple(supercell)
    else:
        print('no supercell provided - assume [1, 1, 1]')
    '''
    num_atoms_in_primitive_cell = len(initial_struct_atoms.get_atomic_numbers())
    print(f'natoms total: {num_atoms_in_primitive_cell}')
    file_conversion.write_file(initial_struct_atoms, lammps_structure_file, out_type='lammps', style=None)
    #file_conversion.convert_and_regularize_file(input_struct_path, lammps_structure_file)

    
    temp_dump = config['output_dir'] / 'steps' / 'final_conf.dump'
    for n_relax in range(1, number_relaxations):
        lammps_dump_file=config['output_dir'] / 'steps' / ('struct_%d_relax_%d.dump'%(number_structure,n_relax))
        lammps_ctrl_file=config['output_dir'] / 'steps' / ('struct_%d_relax_%d.lammps_ctrl'%(number_structure,n_relax))
        lammps_structure_file=config['output_dir'] / 'steps' / ('struct_%d_relax_%d.lammps_struct'%(number_structure, n_relax))
        lammps_out_file=config['output_dir'] / 'steps' / ('struct_%d_relax_%d.out'%(number_structure, n_relax))
        # do relax and vc-relax(volume changes) one after another
        if n_relax % 2 == 0:
            relax_type = 'relax'
        else:
            relax_type = 'vc-relax'

        filename_ctrl = create_relax_ctrl.create_relax_ctrl(lammps_structure_file, lammps_dump_file, lammps_ctrl_file, config,
                                            minimize=relax_type, threshold=thresholds[(n_relax - 1) // 2], keep_cuboidal=config['[1_relax]keep_cuboidal'])
        
        cmd = f"srun -n {config['NTASKS']} {config['path_lammps']} -in {filename_ctrl} > {lammps_out_file}" 
        print(cmd, flush=True)
        os.system(cmd)

        print('n_relax=', n_relax, ' type=%s' % (relax_type))
        next_lammps_structure_file=config['output_dir'] / 'steps' / ('struct_%d_relax_%d.lammps_struct'%(number_structure, n_relax+1))
        cmd = f'tail -n {num_atoms_in_primitive_cell+9} {lammps_dump_file} > {temp_dump}'
        print(cmd, flush=True)
        os.system(cmd)
        file_conversion.convert_and_regularize_file(temp_dump, next_lammps_structure_file)
    #convert final one to vasp
    file_conversion.convert_and_regularize_file(next_lammps_structure_file, config['output_dir'] / 'relaxed_structure_starting_point.POSCAR', out_type='vasp')
    t_end = datetime.datetime.now()

    meta_data['final_density (g/cm3)'] = file_conversion.get_density(file_conversion.read_reg(next_lammps_structure_file))
    with open(config['output_dir'] / 'metadata.txt', 'w', encoding='utf-8') as fp:           
        json.dump(meta_data, fp, indent=4)

    print('end-time:', t_end)
    print('elapsed time:', (t_end-t_start))
   