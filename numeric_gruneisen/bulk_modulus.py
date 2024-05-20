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
import glob
import datetime
import create_relax_ctrl
from numeric_gruneisen import get_filename
import shutil
from numpy.polynomial import polynomial
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.cm as cm
from matplotlib import gridspec
from matplotlib import rc
from matplotlib import rcParams

hplanck = 6.6261e-34
kb = 1.380649e-23

eJ=1.602176487e-19
a_to_m=1.0e-10
amu=1.660599e-27
freq_to_wthz=np.sqrt(eJ/amu/a_to_m**2)/1.0e12
hbar = 1.054571e-34

c_cyan_nature=np.array([68/255, 128/255, 244/255,1.])
c_red=np.array([206./255., 30./255., 65./255.,1.0])
dark_red='#9e0b00'
c_green_nature='#63cf8c'#np.array([96./255., 172./255., 63./255.,1.0])
c_blue_nature=np.array([54./255., 79./255., 156./255.,1.0])
ec=[c_blue_nature,c_green_nature,c_red, 'yellow', 'purple']
plt.rc('axes.formatter', useoffset=False)

'''
B := d^2U/dV^2
'''
def bulk_modulus_helper(volume_strains, volume, energies, natoms=1):
    fig, ax = plt.subplots()
    u0 = np.min(energies)
    energies -= u0
    ax.scatter(volume_strains * volume, np.array(energies), s=20, color=c_green_nature, label='numeric data')
    fitquad = polynomial.polyfit(volume_strains * volume, energies, 2)
    xx = np.linspace(np.min(volume_strains * volume), np.max(volume_strains * volume))
    poly = polynomial.Polynomial(fitquad)
    print(fitquad, poly)
    plt.plot(xx, poly(xx), label='quadratic fit')
    plt.legend()
    B = 2 * volume * poly.coef[2] * eJ * 1e30 /1e9
    ax.set_xlabel('$\Delta V (Å^3)$')
    ax.set_ylabel('$\Delta E (eV)$')
    #plt.grid()
    plt.tick_params(direction='in')
    plt.savefig('6_numeric_gruneisen/bulk_modulus.png', dpi=600, bbox_inches="tight")
    ax.set_title(f'Bulk modulus={B:.3f} GPa,  V={volume} $A^3$')
    plt.savefig('6_numeric_gruneisen/bulk_modulus2.png', dpi=600, bbox_inches="tight")
    return B

def compute_bulk_modulus():
    volume_strains = np.array(config["[6_numeric_gruneisen]volume_strains"])
    atoms_ref = file_conversion.read_reg('relaxed_model/POSCAR')
    atoms_nonaffine_energies = []
    volume_ref = atoms_ref.get_volume()
    
    #n_relax = 18 ## TODO: watch out for this hardcoded constant
    for i, vs in enumerate(volume_strains):
        name_file = get_filename(vs)
        files = glob.glob(f'6_numeric_gruneisen/affine/{name_file}RELAX/steps/struct_{i}_relax_*.out')
        n_relax = [int(x.split('_')[-1][:-4]) for x in files]
        relax_out_file = files[np.argmax(n_relax)]
        print('relax out file: ', relax_out_file)
        #relax_out_file = f'6_numeric_gruneisen/affine/{name_file}RELAX/steps/struct_{i}_relax_{n_relax}.out'
        print("reading energy from output file:" + relax_out_file, flush = True)
        with open(relax_out_file) as f:
            lines = f.readlines()
            ind = -1
            for i, l in enumerate(lines):
                if 'Energy initial, next-to-last, final' in l:
                    ind = i+1
                    break
            atoms_nonaffine_energies.append(list(map(float, lines[ind].split()))[2]) # read energy from out file
    return bulk_modulus_helper(volume_strains, volume_ref, atoms_nonaffine_energies, len(atoms_ref.get_atomic_numbers()))

if __name__ == '__main__':
    compute_bulk_modulus()