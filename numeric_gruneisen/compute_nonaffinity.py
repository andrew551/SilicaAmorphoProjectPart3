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
import itertools
from numeric_gruneisen import get_filename
import matplotlib.pyplot as plt
from numpy.polynomial import polynomial
import math
from matplotlib.offsetbox import AnchoredText

def wrap_vector(v, cell):
    vec_eqs = [v + cell.T @ np.array([i, j, k]) for (i, j, k) in itertools.product([-1, 0, 1], repeat=3)]

    norms = np.linalg.norm(vec_eqs, axis = 1)
    return vec_eqs[np.argmin(norms)]

if __name__ == '__main__':
    print(config)

    volume_strains = config["[6_numeric_gruneisen]volume_strains"]
    deltaSi_list = []
    deltaO_list = []
    for vs in volume_strains:
        scale_factor = (1+vs)**(1/3)
        atoms_affine = file_conversion.read_reg('6_numeric_gruneisen/affine/'+get_filename(vs))
        atoms_nonaffine = file_conversion.read_reg('6_numeric_gruneisen/nonaffine/'+get_filename(vs))



        deltas_Si = atoms_nonaffine[atoms_affine.get_atomic_numbers() == 14].get_positions() - atoms_affine[atoms_affine.get_atomic_numbers() == 14].get_positions()
        deltas_O = atoms_nonaffine[atoms_affine.get_atomic_numbers() == 8].get_positions() - atoms_affine[atoms_affine.get_atomic_numbers() == 8].get_positions()
        curr_cell = atoms_affine.get_cell()
        deltas_Si = np.array([wrap_vector(v, curr_cell) for v in deltas_Si])
        deltas_O = np.array([wrap_vector(v, curr_cell) for v in deltas_O])

        deltaO_list.append(deltas_O)
        deltaSi_list.append(deltas_Si)

    deltaO_list = np.array(deltaO_list)
    deltaSi_list = np.array(deltaSi_list)

    print(deltaO_list.shape)

    model_unstrained = file_conversion.read_reg('relaxed_model/POSCAR')
    n_Si = np.sum(model_unstrained.get_atomic_numbers() == 14)
    n_O = np.sum(model_unstrained.get_atomic_numbers() == 8)
    l0 = (model_unstrained.get_volume() / len(model_unstrained)) ** (1/3)
    print('l0=', l0)

    mag_deltasO = np.linalg.norm(deltaO_list, axis=2) * np.sign(volume_strains).reshape(-1, 1) / l0 if n_O else None
    mag_deltasSi = np.linalg.norm(deltaSi_list, axis=2) * np.sign(volume_strains).reshape(-1, 1) / l0

    coeff_magO = polynomial.polyfit(volume_strains, mag_deltasO, 1) if n_O else None
    coeff_magSi = polynomial.polyfit(volume_strains, mag_deltasSi, 1)
    print(coeff_magO)
    print(coeff_magSi)
    mean_coeffO  = np.mean(coeff_magO[1, :]) if n_O else math.nan
    mean_coeffSi  = np.mean(coeff_magSi[1, :])

    txt_nonaff_Si = '$\\Gamma_{Si} = ' + f'{mean_coeffSi:.3f}' + '$'
    txt_nonaff_O = '$\\Gamma_{O} = ' + f'{mean_coeffO:.3f}' + '$'
    txt_non_all = '$\\Gamma = ' + f'{(mean_coeffO * n_O + mean_coeffSi * n_Si) / (n_O+n_Si):.3f}' + '$'
    txt_comp = txt_nonaff_Si + (('\n' + txt_nonaff_O + '\n' + txt_non_all) if n_O else '')
    with open("6_numeric_gruneisen/nonaffinity_data.txt", 'w') as f:
        f.write(f'O {mean_coeffO}\n')
        f.write(f'Si {mean_coeffSi}\n')
        f.write(f'all {(mean_coeffO * n_O + mean_coeffSi * n_Si) / (n_O+n_Si)}\n')


    '''
    plt.plot(volume_strains, deltaO_list[:, 40, 2])
    plt.plot(volume_strains, deltaO_list[:, 200, 0])
    plt.savefig("deltas.png")

    '''

    rng = np.random.default_rng()

    indecesO = rng.integers(n_O, size=30) if n_O else []
    indecesSi = rng.integers(n_Si, size=30)
    for _ind, i in enumerate(indecesO):
        if not _ind == 0:
            plt.plot(volume_strains, mag_deltasO[:, i], color = 'blue')
        else:
            plt.plot(volume_strains, mag_deltasO[:, i], color = 'blue', label = 'O')
    
    for _ind, i in enumerate(indecesSi):
        if not _ind == 0:
            plt.plot(volume_strains, mag_deltasSi[:, i], color = 'red', linestyle='--')
        else:
            plt.plot(volume_strains, mag_deltasSi[:, i], color = 'red', linestyle='--', label = 'Si')
    plt.plot(volume_strains, volume_strains, color='black', linestyle = ':')
    plt.xlabel('$\\epsilon$')
    plt.ylabel('$\\vert u \\vert \\slash L_0$')
    plt.legend()
    anchored_text = AnchoredText(txt_comp, loc=4)
    plt.gca().add_artist(anchored_text)
    plt.tick_params(direction='in')
    plt.savefig("sdeltas.pdf", dpi=300, bbox_inches="tight", transparent=True)

    plt.savefig("sdeltas.png", dpi=600, bbox_inches="tight", transparent=False)

    