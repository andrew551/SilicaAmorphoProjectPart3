import os
import sys
sys.path.insert(1, os.path.join(sys.path[0], '..')) # "hack" to add the parent directory to path
import makeconfig # this line needs to be after the sys.path.insert line
config = makeconfig.config()
import numpy as np
import scipy
import datetime
from pathlib import Path

def prepare_output_folder(config):
    config['starttime'] = datetime.datetime.now().strftime("%Y%m%d%H%M%S")
    config['output_dir'] = config['output_dir_base'] / 'pydiag_output' / config['starttime'] # note: we are using pathlib.Path object for the paths
    print('output_dir:', config['output_dir'])
    os.makedirs(config['output_dir'])


def do_work(input_file, output_dir):
    output_dir = Path(output_dir)
    dmat_ind = np.loadtxt(input_file)
    print(dmat_ind.shape)
    maxes = np.max(dmat_ind, axis = 0)
    minimums = np.min(dmat_ind, axis = 0)
    print(minimums, maxes, flush=True)
    n = int(maxes[0])
    m = int(minimums[0])
    if not n == int(maxes[1]) or not m == int(minimums[0]):
        raise Exception("not good input!")
    dmat = np.zeros((n, n))
    # allow counting from either 0 or 1 (or any other number)
    dmat_ind[:, :2] -= m
    dmat[dmat_ind[:, 0].astype(int), dmat_ind[:, 1].astype(int)] = dmat_ind[:, 2]
    #prepare_output_folder(config)
    config['output_dir'] = output_dir
    #eigenvalues, eigenvectors = np.linalg.eigh(dmat)
    eigenvalues, eigenvectors = scipy.linalg.eigh(dmat) # eigenvalues in ascending sorted order

    #sort_inds = np.argsort(eigenvalues)
    #eigenvalues = eigenvalues[sort_inds]
    #eigenvectors = eigenvectors[:, sort_inds]
    print(f"saving to: {config['output_dir'] / 'eigenvalues.dat'}")
    np.savetxt(config['output_dir'] / 'eigenvalues.dat', eigenvalues)
    np.savetxt(config['output_dir'] / 'eigenvectors.dat', eigenvectors)

if __name__ == '__main__':
    #input_file = '/mnt/scratch2/q13camb_scratch/adps2/sourecode/SilicaAmorphoProjectPart3/rap2sparse/dmat.dat'
    input_file = str(Path('2_fc2/dmat.dat').resolve())
    output_dir = Path('3_diagonalise_python').resolve()
    do_work(input_file, output_dir)