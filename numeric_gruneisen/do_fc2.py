import os
import sys
sys.path.insert(1, os.path.join(sys.path[0], '..')) # "hack" to add the parent directory to path
import makeconfig # this line needs to be after the sys.path.insert line
config = makeconfig.config()
import numpy as np
import scipy
import datetime
from pathlib import Path
import file_conversion
import json
import datetime
import create_relax_ctrl
from numeric_gruneisen import get_filename
import shutil

if __name__ == '__main__':
    print(config)
    input_file = str(Path('relaxed_model/POSCAR').resolve())
    os.makedirs("6_numeric_gruneisen/affine", exist_ok=True)
    os.makedirs("6_numeric_gruneisen/nonaffine", exist_ok=True)
    print(input_file)
    atoms = file_conversion.read_reg(input_file)

    ## applied volume strains
    volume_strains = np.linspace(-3e-3, 3e-3, 7)
    original_cell = atoms.get_cell()
    # need to copy config.json to deeper folder
    shutil.copyfile('config.json', '6_numeric_gruneisen/affine/config.json')
    shutil.copyfile('config.json', '6_numeric_gruneisen/nonaffine/config.json')
    for vs in volume_strains:
        name_file = get_filename(vs)
        print(f'processing {name_file} ...')
        affine_path = '6_numeric_gruneisen/affine/'+name_file
        nonaffine_path = '6_numeric_gruneisen/nonaffine/'+name_file
        affine_folder = f'6_numeric_gruneisen/affine/{name_file}FC2'
        nonaffine_folder = f'6_numeric_gruneisen/nonaffine/{name_file}FC2'
        os.makedirs(affine_folder, exist_ok=True)
        os.makedirs(nonaffine_folder, exist_ok=True)
        shutil.copyfile(affine_path, affine_folder+'/POSCAR')
        shutil.copyfile(nonaffine_path, nonaffine_folder+'/POSCAR')
        command_fc2 = f'cd {affine_folder}\n\
pwd\n\
echo "start neighbour calculations"\n\
python3 "$AMORPHO_PATH/rap2sparse/calc_NL.py" > neighbour_work.out\n\
echo "start forces calculation"\n\
mpirun -np $SLURM_NTASKS python3 "$AMORPHO_PATH/rap2sparse/RAP2_sparse.py" > rap2_work.out\n\
echo "start collection"\n\
python3 "$AMORPHO_PATH/rap2sparse/collect_RAP2_sparse.py" > collect_work2.out\n\
echo "start conversion"\n\
python3 "$AMORPHO_PATH/rap2sparse/sparse_spRAP_convert.py"\n\
echo "done affine"\n\
now=$(date)\n\
echo "end time : $now"\n\
cd ..\n\
cd ..'
        os.system(command_fc2)
        command_fc2_nonaffine = f'cd {nonaffine_folder}\n\
echo "start neighbour calculations"\n\
python3 "$AMORPHO_PATH/rap2sparse/calc_NL.py" > neighbour_work2.out\n\
echo "start forces calculation"\n\
mpirun -np $SLURM_NTASKS python3 "$AMORPHO_PATH/rap2sparse/RAP2_sparse.py" > rap2_work2.out\n\
echo "start collection"\n\
python3 "$AMORPHO_PATH/rap2sparse/collect_RAP2_sparse.py" > collect_work22.out\n\
echo "start conversion"\n\
python3 "$AMORPHO_PATH/rap2sparse/sparse_spRAP_convert.py"\n\
echo "done nonaffine"\n\
now=$(date)\n\
echo "end time : $now"\n\
cd ..\n\
cd ..'
        os.system(command_fc2_nonaffine)
    