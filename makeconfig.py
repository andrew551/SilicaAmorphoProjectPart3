from pathlib import Path
import os

_config = {
    'path_ACE_potential' : '/mnt/scratch2/q13camb_scratch/POTENTIALS/sio2/ACE/',
    'path_lammps' : '/mnt/userapps/q13camb_apps/lammps/build/lmp',
    'material' : 'SiO2',
    'output_dir_base' : '/users/asmith/grun_out',
    #'output_dir_base' : '/mnt/scratch2/q13camb_scratch/adps2/output_folder1',
    'path_venv': '/users/asmith/programs/kelvenv/bin/activate',
    '[fc2]_FC2_cutoff':12, # FC2 cutoff for FC2
    '[fc2]_Force_cutoff':12, # force cutoff for FC" (what's the difference?)
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

    x = os.getenv('OMP_NUM_THREADS', '1')
    conf_copy['OMP_NUM_THREADS'] = int(x) if x else 1
    y = os.getenv('NTASKS', '1')
    conf_copy['NTASKS'] = int(y) if y else 1
    print(f'DEBUG: x, y = {x}, {y}')
    return conf_copy