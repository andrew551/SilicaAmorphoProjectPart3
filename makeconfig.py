from pathlib import Path
import os
import json

_config = {
    'potential' : 'ACE_Chuck', # ACE_Deringher, or ACE_Chuck, or GAP_Si
    'path_ACE_potential_Chuck' : '/mnt/scratch2/q13camb_scratch/POTENTIALS/sio2/ACE/',
    'path_ACE_potential_Deringher' : '/mnt/scratch2/q13camb_scratch/adps2/ACE_POTENTIAL_SIOX/',
    'path_GAP_Si': '/mnt/scratch2/q13camb_scratch/adps2/GAP_POTENTIAL_SI/gp_iter6_sparse9k.xml',
    'GAP_Si_label': 'GAP_2017_6_17_60_4_3_56_165',
    'path_lammps' : '/mnt/userapps/q13camb_apps/lammps/build',
    'material' : 'SiO2',
    #'output_dir_base' : '/users/asmith/grun_out',
    #'output_dir_base' : '/mnt/scratch2/q13camb_scratch/adps2/output_folder1',
    #'output_dir_base' : '/mnt/scratch2/q13camb_scratch/adps2/output_folder_chik5001',
    'output_dir_base': '/mnt/scratch2/q13camb_scratch/adps2/quenched6000K/1_relax',
    #'path_venv': '/users/asmith/programs/kelvenv/bin/activate',
    '[fc2]_FC2_cutoff':12, # FC2 cutoff for FC2
    '[fc2]_Force_cutoff':12, # force cutoff for FC2 (what's the difference?)
    'keep_cuboidal':True,
    '[fc3_batch]n_batches': 6,
}

def is_path_like(x):
    if isinstance(x, str) and '/' in x:
        try:
            _ = Path(x)
            return True
        except Exception:
            return False
    return False

def config():
    conf_copy = _config.copy()
    # convert strings to Path objects
    for k, v in conf_copy.items():
        if is_path_like(v):
            conf_copy[k] = Path(v)
    try:
        with open('config.json', encoding='utf-8') as f:
            conf_copy.update(json.load(f))
    except Exception: # try parent folder
        with open('../config.json', encoding='utf-8') as f:
            conf_copy.update(json.load(f))
    x = os.getenv('OMP_NUM_THREADS', '1')
    conf_copy['OMP_NUM_THREADS'] = int(x) if x else 1
    y = os.getenv('NTASKS', '1')
    conf_copy['NTASKS'] = int(y) if y else 1
    conf_copy['AMORPHO_PATH'] = os.getenv('AMORPHO_PATH', '')
    conf_copy['path_lammps'] = os.getenv('LAMMPS_PATH', conf_copy['path_lammps']) + '/lmp'
    #print(f'DEBUG: x, y = {x}, {y}')
    return conf_copy

def get_potential_command(config):
    if config['potential'] == 'ACE_Deringher':
        return f'pair_style  pace\n\
pair_coeff  * * {config["path_ACE_potential_Deringher"]}/SiOx_potential.yace O Si\n'
    elif config['potential'] == 'ACE_Chuck':
        return f'pair_style      hybrid/overlay pace table spline 6000\n\
pair_coeff      * * pace {config["path_ACE_potential_Chuck"]}/SiO2-4_24-20-16-12.yace O Si\n\
pair_coeff      1 1 table {config["path_ACE_potential_Chuck"]}/SiO2-4_24-20-16-12_pairpot.table O_O \n\
pair_coeff      1 2 table {config["path_ACE_potential_Chuck"]}/SiO2-4_24-20-16-12_pairpot.table O_Si\n\
pair_coeff      2 2 table {config["path_ACE_potential_Chuck"]}/SiO2-4_24-20-16-12_pairpot.table Si_Si\n'
    elif config['potential'] == 'GAP_Si':
        # 	            'pair_coeff * * %s/sio2_potential_data/potential/silica_gap.xml "Potential xml_label=GAP_2021_4_19_120_7_32_55_336" 8 14'%(path_GAP_potential)]
        return f'pair_style quip\n\
pair_coeff * * {config["path_GAP_Si"]} \"Potential xml_label={config["GAP_Si_label"]}\" 14\n'

    else:
        raise Exception(f"invalid potential name {config['potential']}")

