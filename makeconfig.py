from pathlib import Path
import os
import json

_config = {
    'potential' : 'ACE_Chuck', # ACE_Deringher, or ACE_Chuck, or GAP_Si
    'path_ACE_potential_Chuck' : '/mnt/scratch2/q13camb_scratch/POTENTIALS/sio2/ACE/',
    'path_ACE_potential_Deringher' : '/mnt/scratch2/q13camb_scratch/adps2/ACE_POTENTIAL_SIOX/',
    'path_GAP_Si': '/mnt/scratch2/q13camb_scratch/adps2/GAP_POTENTIAL_SI/gp_iter6_sparse9k.xml',
    'path_GAP_SiO2': '/mnt/scratch2/q13camb_scratch/POTENTIALS/sio2/GAP/sio2_potential_data/potential/silica_gap.xml',
    'GAP_Si_label': 'GAP_2017_6_17_60_4_3_56_165',
    'GAP_SiO2_label': 'GAP_2021_4_19_120_7_32_55_336',
    'path_ACE_cheap_CHUCK_1' : '/mnt/scratch2/q13camb_scratch/adps2/ACE_NEW_CHUCK/1',
    'path_ACE_cheap_CHUCK_2' : '/mnt/scratch2/q13camb_scratch/adps2/ACE_NEW_CHUCK/2',
    'path_ACE_cheap_CHUCK_3' : '/mnt/scratch2/q13camb_scratch/adps2/ACE_NEW_CHUCK/3',
    'path_ACE_cheap_CHUCK_4' : '/mnt/scratch2/q13camb_scratch/adps2/ACE_NEW_CHUCK/4',

    'path_lammps' : '/mnt/userapps/q13camb_apps/lammps/build',
    'material' : 'SiO2',
    
    #'[1_relax]thresholds': [1e-9,1e-10,1e-11,1e-12,1e-12,1e-12,1e-13,1e-13,1e-13,1e-14,1e-14,1e-14,1e-15,1e-15, 8e-16, 5e-16, 3e-16, 1e-16],
    #'[1_relax]keep_cuboidal':True,
    #'[fc2]_FC2_cutoff':12, # FC2 cutoff for FC2
    #'[fc2]_Force_cutoff':12, # force cutoff for FC2 (what's the difference?)
}

def is_path_like(x):
    if isinstance(x, str) and '/' in x:
        try:
            _ = Path(x)
            return True
        except Exception:
            return False
    return False

'''
get the full config dictionary
combines enviromental variables, config.json file, and default parameters

returns: dictionary
'''
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
    conf_copy['path_lammps'] = os.getenv('LAMMPS_PATH', str(conf_copy['path_lammps'])) + '/lmp'
    #print(f'DEBUG: x, y = {x}, {y}')
    return conf_copy

'''
get the lammps command to run each supported potential
to add a new potential, simply update this function and anneal, relax, fc2, fc3 should all work
'''
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
        return f'pair_style quip\n\
pair_coeff * * {config["path_GAP_Si"]} \"Potential xml_label={config["GAP_Si_label"]}\" 14\n'
    elif config['potential'] == 'GAP_SiO2':
        return f'pair_style quip\n\
pair_coeff * * {config["path_GAP_SiO2"]} \"Potential xml_label={config["GAP_SiO2_label"]}\" 8 14\n'

    elif config['potential'] == 'ACE_CHEAP_CHUCK_1':
        return f'pair_style      hybrid/overlay pace table spline 6000\n\
pair_coeff      * * pace {config["path_ACE_cheap_CHUCK_1"]}/SiO2-3-12.yace O Si\n\
pair_coeff      1 1 table {config["path_ACE_cheap_CHUCK_1"]}/SiO2-3-12_pairpot.table O_O \n\
pair_coeff      1 2 table {config["path_ACE_cheap_CHUCK_1"]}/SiO2-3-12_pairpot.table O_Si\n\
pair_coeff      2 2 table {config["path_ACE_cheap_CHUCK_1"]}/SiO2-3-12_pairpot.table Si_Si\n' 
    
    elif config['potential'] == 'ACE_CHEAP_CHUCK_2':
        return f'pair_style      hybrid/overlay pace table spline 6000\n\
pair_coeff      * * pace {config["path_ACE_cheap_CHUCK_2"]}/SiO2-3-18.yace O Si\n\
pair_coeff      1 1 table {config["path_ACE_cheap_CHUCK_2"]}/SiO2-3-18_pairpot.table O_O \n\
pair_coeff      1 2 table {config["path_ACE_cheap_CHUCK_2"]}/SiO2-3-18_pairpot.table O_Si\n\
pair_coeff      2 2 table {config["path_ACE_cheap_CHUCK_2"]}/SiO2-3-18_pairpot.table Si_Si\n' 

    elif config['potential'] == 'ACE_CHEAP_CHUCK_3':
        return f'pair_style      hybrid/overlay pace table spline 6000\n\
pair_coeff      * * pace {config["path_ACE_cheap_CHUCK_3"]}/SiO2-4-12.yace O Si\n\
pair_coeff      1 1 table {config["path_ACE_cheap_CHUCK_3"]}/SiO2-4-12_pairpot.table O_O \n\
pair_coeff      1 2 table {config["path_ACE_cheap_CHUCK_3"]}/SiO2-4-12_pairpot.table O_Si\n\
pair_coeff      2 2 table {config["path_ACE_cheap_CHUCK_3"]}/SiO2-4-12_pairpot.table Si_Si\n' 

    elif config['potential'] == 'ACE_CHEAP_CHUCK_4':
        return f'pair_style      hybrid/overlay pace table spline 6000\n\
pair_coeff      * * pace {config["path_ACE_cheap_CHUCK_4"]}/SiO2-4_24-20-16-12.yace O Si\n\
pair_coeff      1 1 table {config["path_ACE_cheap_CHUCK_4"]}/SiO2-4_24-20-16-12_pairpot.table O_O \n\
pair_coeff      1 2 table {config["path_ACE_cheap_CHUCK_4"]}/SiO2-4_24-20-16-12_pairpot.table O_Si\n\
pair_coeff      2 2 table {config["path_ACE_cheap_CHUCK_4"]}/SiO2-4_24-20-16-12_pairpot.table Si_Si\n' 

    else:
        raise Exception(f"invalid potential name {config['potential']}")


def get_atom_dict(config):
    if config['material'] == 'SiO2':
        return {'O':1, 'Si':2}
    elif config['material'] == 'Si':
        return {'Si':1}
    raise Exception(f"unsupported material: {config['material']}")
