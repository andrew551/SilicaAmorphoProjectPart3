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
import shutil
from numpy.polynomial import polynomial
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.cm as cm
from matplotlib import gridspec
from matplotlib import rc
from matplotlib import rcParams
from numpy.polynomial import polynomial

from numeric_gruneisen.compute_gruneisen import match_eigenmodes, compute_material_gruneisen
from solve_SCP import compute_correction, compute_frequency_shifts_smoothed

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

def get_filename(x):
    return f'quartic{(x+1):.3f}'


CUTOFF_FREQUENCY = 5e13
USE_SMOOTH = True

if __name__ == '__main__':
    volume_strains = np.array([-0.001, 0.000, 0.001])
    print('volume strains: ', volume_strains)
    eigenvalues_anh = []
    eigenvalues0 = []
    frequencies0 = []
    frequencies0cm = []
    phi4s = []
    matchings = []
    reference_eigenmodes = np.loadtxt(f'{get_filename(0.0)}/3_diagonalise/eigenvectors.dat')
    reference_eigenvalues = np.loadtxt(f'{get_filename(0.0)}/3_diagonalise/eigenvalues.dat')
    print('matrix shape:', reference_eigenmodes.shape)

    for vs in volume_strains:
        eigenvalues_vs = np.loadtxt(f'{get_filename(vs)}/3_diagonalise/eigenvalues.dat')
        eigenvectors_vs = np.loadtxt(f'{get_filename(vs)}/3_diagonalise/eigenvectors.dat')
        phi4_vs = np.loadtxt(f'{get_filename(vs)}/10_quartic/pet.dat')
        ### delete translational modes ###
        eigenvalues_vs = eigenvalues_vs[3:]
        phi4_vs = phi4_vs[3:, 3:]
        eigenvectors_vs = eigenvectors_vs[3:, 3:]
        matching_i = match_eigenmodes(reference_eigenmodes[3:, 3:], eigenvectors_vs)
        del eigenvectors_vs
        eigenvalues0.append(eigenvalues_vs)
        frequencies0.append(np.sqrt(eigenvalues_vs) * freq_to_wthz * 1e12)
        frequencies0cm.append(np.sqrt(np.abs(eigenvalues_vs)*9.648e27)/18.8e10) # frequencies in cm^-1)
        phi4s.append(phi4_vs)
        matchings.append(matching_i)
    print('f0:', frequencies0[:10])
    Ts = np.arange(1, 1000)
    gamma_mat_T = []
    cv_mat_T = []
    mid = len(volume_strains) // 2

    index_cutoff = np.argmax(frequencies0[mid] > CUTOFF_FREQUENCY)
    print("index cutoff:", index_cutoff, flush=True)
    plt.clf()
    for T in Ts:
        evs_anh_T = []
        for i, vs in enumerate(volume_strains):
            if USE_SMOOTH:
                frequency_shifts = compute_frequency_shifts_smoothed(phi4s[i], frequencies0[i], T)
                frequency_shifts[np.isnan(frequency_shifts)] = 0
                frequency_shifts[:index_cutoff] = frequencies0[mid][:index_cutoff] / frequencies0[mid][index_cutoff + 1] * frequency_shifts[index_cutoff + 1]
                #new_freqs = frequencies0[i] + frequency_shifts
                eigenvalues_vs_anh = ((frequencies0[i] + frequency_shifts) / 1e12 / freq_to_wthz)**2
                if T %10 == 0:
                    print(T, i, frequency_shifts)
                    print(frequencies0[i])
            else:
                new_freqs, correction, eigenvalues_vs_anh = compute_correction(phi4s[i], frequencies0[i], T)
            evs_anh_T.append(eigenvalues_vs_anh[matchings[i]])
            #evs_anh_T.append(eigenvalues0[i][matchings[i]])
        # compute mode gruneisen
        evs_anh_T = np.array(evs_anh_T)
        #print('evsanhshape', evs_anh_T.shape)
        lin_coeff = polynomial.polyfit(volume_strains[1:], evs_anh_T[1:], 1)

        gamma_mode = -1/2 * 1/reference_eigenvalues[3:] * lin_coeff[1, :]
        if T == 100 or T == 300 or T == 1:
            plt.scatter(frequencies0cm[mid], gamma_mode, label=f'T={T}K')
        gamma_mat, cv_mat = compute_material_gruneisen([T], frequencies0cm[i], gamma_mode, evs_anh_T[mid, :])
        gamma_mat_T.append(gamma_mat[0])
        cv_mat_T.append(cv_mat[0])
    plt.legend()
    plt.savefig("10_quartic/modgrun.png", dpi=400)
    plt.clf()
    plt.plot(Ts, gamma_mat_T)
    plt.savefig('10_quartic/gruneisen_quartic.png', dpi=400)

    np.savetxt('10_quartic/gamma_ang.dat', gamma_mat_T)
    np.savetxt('10_quartic/cv_mode_mat.dat', cv_mat_T)