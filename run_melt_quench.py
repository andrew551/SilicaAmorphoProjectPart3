import sys, os
import datetime
from pathlib import Path
import json

# import local files
import file_conversion
import create_anneal_ctrl


# read the config file
config = {}
with open('config.json', encoding='utf-8') as f:
    config.update(json.load(f))

print(config)
print("\n")

if not config['test_anneal?']:
    raise Exception("config.json says anneal not needed!")


input_struct_path = config['input_struct_path']
nproc = config['ntasks']
lammps_path = config['lammps_path']

output_dir = Path('0_melt_quench').resolve()
config['output_dir'] = output_dir
print(output_dir)
print("\n")


# convert from POSCAR to lammps format
regularised_input_path = f"{output_dir}/POSCAR_regularised.dat"
file_conversion.convert_vasp_to_lammps(input_struct_path, regularised_input_path)


# create anneal ctrl file and save it to a file
command = create_anneal_ctrl.create_anneal_ctrl(regularised_input_path, config)
print(command, flush=True)
print("\n")

lammps_ctrl = f"{output_dir}/lammps_MD_test.ctrl"
with open(lammps_ctrl, 'w') as f:
        f.write(command)

lammps_out = f"{output_dir}/lammps_MD_test.out"


# launch the MD
anneal_command = f"srun -n {nproc} {lammps_path}/lmp -in {lammps_ctrl} > {lammps_out}"
print(anneal_command, flush=True)
os.system(anneal_command)