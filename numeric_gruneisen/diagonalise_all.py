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
from diagonalise_python import py_diag, convert_freq

if __name__ == '__main__':
    volume_strains = np.linspace(-3e-3, 3e-3, 7)
    for vs in volume_strains:
        name_file = get_filename(vs)
        print(f'processing {name_file} ...')
        affine_folder = f'6_numeric_gruneisen/affine/{name_file}FC2'
        nonaffine_folder = f'6_numeric_gruneisen/nonaffine/{name_file}FC2'
        py_diag.do_work(f'{affine_folder}/dmat.dat', affine_folder)
        convert_freq.do_work(affine_folder)
        py_diag.do_work(f'{nonaffine_folder}/dmat.dat', nonaffine_folder)
        convert_freq.do_work(nonaffine_folder)