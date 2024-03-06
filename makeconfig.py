from pathlib import Path
_config = {
    'path_ACE_potential' : '/users/asmith/ACE',
    'path_lammps' : '/users/asmith/programs/lammps/build/lmp',
    'material' : 'SiO2',
    'output_dir_base' : '/users/asmith/grun_out',
    'path_lammps_src': '/users/asmith/programs/lammps/src/',
    'path_venv': '/users/asmith/programs/kelvenv/bin/activate',
}

def config():
    conf_copy = _config.copy()
    # convert strings to Path objects
    for k, v in conf_copy.items():
        if isinstance(v, str) and '/users/' in v:
            conf_copy[k] = Path(v)
    return conf_copy