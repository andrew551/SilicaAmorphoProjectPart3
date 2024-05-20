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
try:
    from numeric_gruneisen.numeric_gruneisen import get_filename
    import numeric_gruneisen.bulk_modulus
except Exception:
    from numeric_gruneisen import get_filename
    import bulk_modulus
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

#return which mode in b each mode in a corresponds to
def match_eigenmodes_old(a, b):
    if not a.shape == b.shape:
        raise ValueError("Input shape mismatch")
    thresh = 0.95
    matching = np.zeros(a.shape[0], dtype=int)
    matched = set()
    for i in range(a.shape[0] - 1, -1, -1):
        best = -999
        best_j = -1
        for j in range(b.shape[0] - 1, -1, -1):
            if j in matched:
                continue
            alignment = np.abs(np.dot(a[:, i], b[:, j]))
            if alignment > best:
                best = alignment
                best_j = j
            if alignment > thresh:
                break # choose this matching quickly
        matched.add(best_j)
        matching[i] = best_j
        print(i, best, best_j)
    return matching

def match_eigenmodes_full(a, b):
    if not a.shape == b.shape:
        raise ValueError("Input shape mismatch")
    dot_products = np.abs(np.einsum("ki,kj->ij", a, b))
    n = a.shape[0]
    min_align = 0.2
    matching = np.zeros(a.shape[0], dtype=int)
    rem_a = set(np.arange(n))
    rem_b = set(np.arange(n))
    for i in range(n):
        max_match = np.unravel_index(np.argmax(dot_products), dot_products.shape)
        if dot_products[max_match] > min_align:
            matching[max_match[0]] = max_match[1]
            dot_products[max_match[0], :] = -2
            dot_products[:, max_match[1]] = -2
            rem_a.remove(max_match[0])
            rem_b.remove(max_match[1])
        else:
            break
    for (ai, bi) in zip(rem_a, rem_b):
        matching[ai] = bi
    return matching
    

def match_eigenmodes(a, b, w=10):
    if not a.shape == b.shape:
        raise ValueError("Input shape mismatch")
    
    n = a.shape[0]
    dot_products = np.zeros(a.shape)
    for i in range(n):
        l0 = max(i-w, 0)
        l1 = min(n, i+w+1)
        xx = np.einsum("k,kj->j", a[:, i], b[:, l0:l1])
        dot_products[i, l0:l1] = np.abs(xx)
    #dot_products = np.abs(np.einsum("ki,kj->ij", a, b))
    args_sorted = np.argsort(dot_products, axis=None)
    min_align = 0.2
    matching = np.zeros(a.shape[0], dtype=int)
    rem_a = set(np.arange(n))
    rem_b = set(np.arange(n))
    for arg in reversed(args_sorted):
        max_match = np.unravel_index(arg, dot_products.shape)
        #print(max_match)
        if dot_products[max_match] < min_align:
            break
        if  max_match[0] in rem_a and max_match[1] in rem_b:
            matching[max_match[0]] = max_match[1]
            rem_a.remove(max_match[0])
            rem_b.remove(max_match[1])
    print('# matched easily: ' + str(n-len(rem_a)))
    if len(rem_a) != len(rem_b):
        raise Exception("runtime error")
    for (ai, bi) in zip(rem_a, rem_b):
        matching[ai] = bi
    return matching

def compute_material_gruneisen(T, frequencies, mode_gruneisen, ev):

    energy_joules1 = frequencies * 0.029979 * hplanck * 1e12 # cm^-1 -> THz -> Joules
    energy_joules = np.sqrt(ev)*freq_to_wthz*1e12*hbar
    print('energy10', energy_joules[:10])
    print('energy10x', energy_joules1[:10])
    grun_T = []
    cv_T = []
    for T_i in T:
        arg = energy_joules / (kb*T_i)
        if T_i == 10:
            print(np.min(arg), np.max(arg))
            print(arg[::50])
            print(arg)
        arg = np.maximum(np.minimum(arg, 100), -100) # clip large values
        mode_occupation = 1/(np.exp(arg) - 1) # Boson mode occupation
        mode_heat_capacity = energy_joules ** 2 / (kb*T_i*T_i) * mode_occupation * (1 + mode_occupation)
        material_gruneisen_i = np.dot(mode_heat_capacity, mode_gruneisen) / np.sum(mode_heat_capacity) # weighted sum of mode Gruneisen by mode heat capacity
        grun_T.append(material_gruneisen_i)
        cv_T.append(np.sum(mode_heat_capacity))
    return np.array(grun_T), np.array(cv_T)
    
    '''
    do i=1,1000
          T=i
          nj=0.0d0 ! boson factor at mode freq []
          do mj=1+3,3*natom
            eng=dsqrt(ww(mj))*freq_to_wthz*1d12*hbar ! in [J]
            arg=eng/(T*kb) ! dimensionless
            if(dabs(arg)>100)arg=100*arg/dabs(arg)
            nj(mj)=1.0d0/(dexp(arg)-1.0d0)
         !   write(*,*)nj(mj)
          enddo

          CV=0.0d0 ! capacity at mode freq positions [J/K]
          do mj=4,3*natom
            CV(mj)=(hbar*dsqrt(ww(mj))*freq_to_wthz*1d12)**2/ &
                   (kb*T*T)*nj(mj)*(1.0d0+nj(mj))  ! [J/K]
        !     write(*,*)cv(mj)
          enddo
    '''

def compute_all():
    atoms_ref = file_conversion.read_reg('relaxed_model/POSCAR')
    volume_strains = np.array(config["[6_numeric_gruneisen]volume_strains"])
    print('volume strains: ', volume_strains)
    eigenvalues_affine = []
    #eigenmodes_affine = []
    eigenvalues_nonaffine = []
    #eigenmodes_nonaffine = []


    reference_eigenmodes = np.loadtxt(f'6_numeric_gruneisen/affine/{get_filename(0.0)}FC2/eigenvectors.dat')
    reference_eigenvalues = np.loadtxt(f'6_numeric_gruneisen/affine/{get_filename(0.0)}FC2/eigenvalues.dat')
    print('matrix shape:', reference_eigenmodes.shape)

    
    for i, vs in enumerate(volume_strains):
        name_file = get_filename(vs)
        print(f'processing {name_file} ...', flush=True)
        affine_folder = f'6_numeric_gruneisen/affine/{name_file}FC2'
        nonaffine_folder = f'6_numeric_gruneisen/nonaffine/{name_file}FC2'
        # open eigenvalues, eigenmodes
        # gamma_k = -1/(2*w^2) * d(w^2)/dV
        eigenvalues_affine_vs = np.loadtxt(affine_folder+'/eigenvalues.dat')
        next_affine_eigenmodes = np.loadtxt(affine_folder+'/eigenvectors.dat')
        matching_affine = match_eigenmodes(reference_eigenmodes, next_affine_eigenmodes)
        del next_affine_eigenmodes
        eigenvalues_nonaffine_vs = np.loadtxt(nonaffine_folder+'/eigenvalues.dat')
        next_nonaffine_eigenmodes= np.loadtxt(nonaffine_folder+'/eigenvectors.dat')
        matching_nonaffine = match_eigenmodes(reference_eigenmodes, next_nonaffine_eigenmodes)
        del next_nonaffine_eigenmodes
        eigenvalues_affine.append(eigenvalues_affine_vs[matching_affine])
        eigenvalues_nonaffine.append(eigenvalues_nonaffine_vs[matching_nonaffine])


    bulk_modulus_computed = bulk_modulus.compute_bulk_modulus()
    print("computed bulk modulus:", bulk_modulus_computed)
    data_affine = np.stack(eigenvalues_affine)
    data_nonaffine = np.stack(eigenvalues_nonaffine)
    print(data_affine.shape)
    print('saving data to ', Path('6_numeric_gruneisen/affine/affine_ev_freq.dat').resolve(), flush=True)
    np.savetxt('6_numeric_gruneisen/affine/affine_ev_freq.dat', data_affine)
    np.savetxt('6_numeric_gruneisen/nonaffine/nonaffine_ev_freq.dat', data_nonaffine)

    #   use strains -0.002 ... 0.002
    skip_grun = 1
    data_affine = data_affine[skip_grun:-skip_grun, :]
    data_nonaffine = data_nonaffine[skip_grun:-skip_grun, :]
    volume_strains = volume_strains[skip_grun:-skip_grun]

    lin_coeff_affine = polynomial.polyfit(volume_strains, data_affine, 1)
    lin_coeff_nonaffine = polynomial.polyfit(volume_strains, data_nonaffine, 1)

    gamma_affine = -1/2 * 1/reference_eigenvalues * lin_coeff_affine[1, :]
    gamma_nonaffine = -1/2 * 1/reference_eigenvalues * lin_coeff_nonaffine[1, :]

    ref_freq = np.sqrt(np.abs(reference_eigenvalues)*9.648e27)/18.8e10 # frequencies in cm^-1
    np.savetxt('6_numeric_gruneisen/affine/gamma_affine.dat', np.c_[ref_freq, gamma_affine])
    np.savetxt('6_numeric_gruneisen/nonaffine/gamma_nonaffine.dat', np.c_[ref_freq, gamma_nonaffine])
    T = np.arange(1, 1001)
    grun_affine, heat_capacity = compute_material_gruneisen(T, ref_freq[3:], gamma_affine[3:], reference_eigenvalues[3:])
    heat_capacity_vol = heat_capacity / (atoms_ref.get_volume()*1e-30)
    specific_heat = heat_capacity_vol / (file_conversion.get_density(atoms_ref)*1e6)
    grun_nonaffine, _ = compute_material_gruneisen(T, ref_freq[3:], gamma_nonaffine[3:], reference_eigenvalues[3:])
    print(grun_affine)
    print(grun_nonaffine)
    np.savetxt('6_numeric_gruneisen/affinegrun_T.dat', np.c_[T, grun_affine])
    np.savetxt('6_numeric_gruneisen/nonaffinegrun_T.dat', np.c_[T, grun_nonaffine])
    np.savetxt('6_numeric_gruneisen/CV_T.dat', np.c_[T, specific_heat])
    np.savetxt('6_numeric_gruneisen/C_per_volume_T.dat', np.c_[T, heat_capacity_vol])

if __name__ == '__main__':
    compute_all()
    exit()

    volume_strains = config["[6_numeric_gruneisen]volume_strains"]
    print('volume strains: ', volume_strains)
    eigenvalues_affine = []
    eigenmodes_affine = []
    eigenvalues_nonaffine = []
    eigenmodes_nonaffine = []

    reference_eigenmodes = np.loadtxt(f'6_numeric_gruneisen/affine/{get_filename(0.0)}FC2/eigenvectors.dat')
    reference_eigenvalues = np.loadtxt(f'6_numeric_gruneisen/affine/{get_filename(0.0)}FC2/eigenvalues.dat')
    print('matrix shape:', reference_eigenmodes.shape)

    atoms_ref = file_conversion.read_reg('relaxed_model/POSCAR')
    atoms_nonaffine_energies = []
    volume_ref = atoms_ref.get_volume()
    for i, vs in enumerate(volume_strains):
        name_file = get_filename(vs)
        print(f'processing {name_file} ...')
        affine_folder = f'6_numeric_gruneisen/affine/{name_file}FC2'
        nonaffine_folder = f'6_numeric_gruneisen/nonaffine/{name_file}FC2'
        # open eigenvalues, eigenmodes
        # gamma_k = -1/(2*w^2) * d(w^2)/dV
        eigenvalues_affine_vs = np.loadtxt(affine_folder+'/eigenvalues.dat')
        eigenmodes_affine.append(np.loadtxt(affine_folder+'/eigenvectors.dat'))
        eigenvalues_nonaffine_vs = np.loadtxt(nonaffine_folder+'/eigenvalues.dat')
        eigenmodes_nonaffine.append(np.loadtxt(nonaffine_folder+'/eigenvectors.dat'))
    
        matching_affine = match_eigenmodes(reference_eigenmodes, eigenmodes_affine[-1])
        matching_nonaffine = match_eigenmodes(reference_eigenmodes, eigenmodes_nonaffine[-1])

        eigenvalues_affine.append(eigenvalues_affine_vs[matching_affine])
        eigenvalues_nonaffine.append(eigenvalues_nonaffine_vs[matching_nonaffine])

        n_relax = 18
        #/mnt/scratch2/q13camb_scratch/adps2/work/silica_DFT192/6_numeric_gruneisen/affine/POSCAR0.997RELAX/steps/struct_0_relax_18.out
        relax_out_file = f'6_numeric_gruneisen/affine/{name_file}RELAX/steps/struct_{i}_relax_{n_relax}.out'
        with open(relax_out_file) as f:
            lines = f.readlines()
            ind = -1
            for i, l in enumerate(lines):
                if 'Energy initial, next-to-last, final' in l:
                    ind = i+1
                    break
            atoms_nonaffine_energies.append(list(map(float, lines[ind].split()))[2]) # read energy from out file

    bulk_modulus_computed = compute_bulk_modulus(volume_strains, volume_ref, atoms_nonaffine_energies)
    
    #exit()

    data_affine = np.stack(eigenvalues_affine)
    data_nonaffine = np.stack(eigenvalues_nonaffine)
    print(data_affine.shape)
    print('saving data to ', Path('6_numeric_gruneisen/affine/affine_ev_freq.dat').resolve())
    np.savetxt('6_numeric_gruneisen/affine/affine_ev_freq.dat', data_affine)
    np.savetxt('6_numeric_gruneisen/nonaffine/nonaffine_ev_freq.dat', data_nonaffine)

    #   use strains -0.002 ... 0.002
    data_affine = data_affine[1:-1, :]
    data_nonaffine = data_nonaffine[1:-1, :]
    volume_strains = volume_strains[1:-1]

    lin_coeff_affine = polynomial.polyfit(volume_strains, data_affine, 1)
    lin_coeff_nonaffine = polynomial.polyfit(volume_strains, data_nonaffine, 1)

    gamma_affine = -1/2 * 1/reference_eigenvalues * lin_coeff_affine[1, :]
    gamma_nonaffine = -1/2 * 1/reference_eigenvalues * lin_coeff_nonaffine[1, :]

    ref_freq = np.sqrt(np.abs(reference_eigenvalues)*9.648e27)/18.8e10 # frequencies in cm^-1
    np.savetxt('6_numeric_gruneisen/affine/gamma_affine.dat', np.c_[ref_freq, gamma_affine])
    np.savetxt('6_numeric_gruneisen/nonaffine/gamma_nonaffine.dat', np.c_[ref_freq, gamma_nonaffine])

    #colors
    cmap = cm.hsv
    matplotlib.rcParams.update({'font.size': 12})
    #matplotlib.rcParams.update({'mathtext.fontset' : 'dejavuserif'})
    matplotlib.rcParams['text.latex.preamble'] = [
        r'\usepackage{siunitx}',
        r'\usepackage{amssymb}',   
        r'\sisetup{detect-all}',   # ...this to force siunitx to actually use your fonts
        r'\usepackage{helvet}',    # set the normal font here
        r'\usepackage{sansmath}',  # load up the sansmath so that math -> helvet
        r'\sansmath'               
    ] 

    fig, ax = plt.subplots()
    ax.set_title('Affine Gruneisen')
    ax.set_xlabel('$\hbar\omega$ $(cm^{-1})$')
    ax.set_ylabel('$\\gamma(\\omega)$')
    ax.scatter(ref_freq[3:], gamma_affine[3:], s=2, color=c_cyan_nature)
    ax.set_ylim((-4, 4))
    plt.grid()
    plt.savefig('6_numeric_gruneisen/affine/fc3grun.png')

    fig, ax = plt.subplots()
    ax.set_title('Non-affine Gruneisen')
    ax.set_xlabel('$\hbar\omega$ $(cm^{-1})$')
    ax.set_ylabel('$\\gamma(\\omega)$')
    ax.scatter(ref_freq[3:], gamma_nonaffine[3:], s=2, color=c_green_nature)
    ax.set_ylim((-4, 4))
    plt.grid()
    plt.savefig('6_numeric_gruneisen/affine/ISgrun.png')

    fig, ax = plt.subplots()
    ax.set_title('Non-affine Gruneisen')
    ax.set_xlabel('$\hbar\omega$ $(cm^{-1})$')
    ax.set_ylabel('$\\gamma(\\omega)$')
    ax.scatter(ref_freq[3:], gamma_nonaffine[3:], s=2, color=c_green_nature)
    ax.set_ylim((-150, 150))
    plt.grid()
    plt.savefig('6_numeric_gruneisen/affine/ISgrun_big.png')

    fig, ax = plt.subplots()
    ax.set_title('Affine Gruneisen')
    ax.set_xlabel('$\hbar\omega$ $(cm^{-1})$')
    ax.set_ylabel('$\\gamma(\\omega)$')
    ax.scatter(ref_freq[3:], gamma_affine[3:], s=2, color=c_cyan_nature)
    ax.set_ylim((-150, 150))
    plt.grid()
    plt.savefig('6_numeric_gruneisen/affine/fc3grun_big.png')

    fig, ax = plt.subplots()
    ax.set_title('Mode Gruneisen')
    ax.set_xlabel('$\hbar\omega$ $(cm^{-1})$')
    ax.set_ylabel('$\\gamma(\\omega)$')
    ax.scatter(ref_freq[3:], gamma_affine[3:], s=2, color=c_cyan_nature, label = 'affine')
    ax.scatter(ref_freq[3:], gamma_nonaffine[3:], s=2, color=c_green_nature, label = 'non-affine')
    ax.set_ylim((-3, 2))
    plt.legend()
    plt.grid()
    plt.savefig('6_numeric_gruneisen/affine/mode_grun.png')

    fig, ax = plt.subplots()
    ax.set_title('IS contribution to Gruneisen')
    ax.set_xlabel('$\hbar\omega$ $(cm^{-1})$')
    ax.set_ylabel('$\\Delta \\gamma(\\omega)$')
    ax.scatter(ref_freq[3:], gamma_nonaffine[3:] - gamma_affine[3:], s=2, color=c_green_nature)
    #ax.set_ylim((-3, 2))
    plt.legend()
    plt.grid()
    plt.savefig('6_numeric_gruneisen/affine/delta_mode_grun.png')

    fig, ax = plt.subplots()
    data_paris_fc3 = np.loadtxt('/mnt/scratch2/q13camb_scratch/adps2/silicon_crys/5_gruneisen/gruneisen_MAIN.dat')
    data_paris_IS = np.loadtxt('/mnt/scratch2/q13camb_scratch/adps2/silicon_crys/5_gruneisen/gruneisen_IS.dat')
    #data_paris_fc3 = np.loadtxt('/mnt/scratch2/q13camb_scratch/adps2/dft_648_silica/5_gruneisen/gruneisen_MAIN.dat')
    #data_paris_IS = np.loadtxt('/mnt/scratch2/q13camb_scratch/adps2/dft_648_silica/5_gruneisen/gruneisen_IS.dat')
    ax.set_title('Mode Gruneisen Comparison (c-Si)')
    ax.set_xlabel('$\hbar\omega$ $(cm^{-1})$')
    ax.set_ylabel('$\\gamma(\\omega)$')
    ax.scatter(data_paris_fc3[:, 0], data_paris_fc3[:, 1], s=2, color = 'yellow', label = 'Paris Code (FC3)')
    ax.scatter(data_paris_IS[:, 0], data_paris_fc3[:, 1] - data_paris_IS[:, 1], s=2, color = 'purple', label = 'Paris Code (FC3 + IS)')
    ax.scatter(ref_freq[3:], gamma_affine[3:], s=2, color=c_cyan_nature, label = 'New method affine')
    ax.scatter(ref_freq[3:], gamma_nonaffine[3:], s=2, color=c_green_nature, label = 'New method non-affine')
    ax.set_ylim((-10, 10))
    plt.legend(fontsize=6)
    plt.grid()
    plt.savefig('6_numeric_gruneisen/affine/mode_grun_paris_comp.png', dpi=600, bbox_inches="tight", transparent=False)
    
    

    #lower_cutoff = ref_freq[6] * 0.029979 * hplanck * 1e12 / kb / 2
    #lower_cutoff = 1
    #print(f'lowest reasonable T: {lower_cutoff}K')
    T = np.arange(1, 1001)
    grun_affine, heat_capacity = compute_material_gruneisen(T, ref_freq[3:], gamma_affine[3:], reference_eigenvalues[3:])
    heat_capacity_vol = heat_capacity / (atoms_ref.get_volume()*1e-30)
    specific_heat = heat_capacity_vol / (file_conversion.get_density(atoms_ref)*1e6)
    grun_nonaffine, _ = compute_material_gruneisen(T, ref_freq[3:], gamma_nonaffine[3:], reference_eigenvalues[3:])
    print(grun_affine)
    print(grun_nonaffine)

    
    grun_paris_fc3, _ = compute_material_gruneisen(T, ref_freq[3:], data_paris_fc3[:, 1], reference_eigenvalues[3:])
    grun_paris_IS, _ = compute_material_gruneisen(T, ref_freq[3:], data_paris_fc3[:, 1] - data_paris_IS[:, 1], reference_eigenvalues[3:])
    print(grun_paris_fc3)
    fig, ax = plt.subplots()
    plt.plot(T[5:],grun_paris_fc3[5:],color='yellow',label='Paris Code (FC3)')
    plt.plot(T[5:], grun_paris_IS[5:], color='purple',label='Paris Code (FC3 + IS)')
    plt.plot(T[5:],grun_affine[5:],color=ec[0],label='affine Gruneisen')
    plt.plot(T[5:], grun_nonaffine[5:], color=ec[1],label='non-affine Gruneisen')
    plt.xlabel(r'$T(K)$')
    plt.ylabel(r'$\gamma$')
    plt.legend(fontsize=6)
    plt.ylim([-8, 3])
    plt.grid()
    fig.savefig("6_numeric_gruneisen/numeric_gruneisen_linear_scale_comp_Paris.png", dpi=600, bbox_inches="tight", transparent=False)
    
    exit()

    '''
    cumul_exp = np.cumsum(grun_nonaffine * heat_capacity_vol / 33.5e9 / 3)
    cumul_exp_a = np.cumsum(grun_affine * heat_capacity_vol / 33.5e9 / 3)
    fig = plt.figure()

    plt.axhline(y=0, color='black')
    plt.plot(T, cumul_exp_a, label = 'affine', color=ec[0])
    plt.plot(T, cumul_exp, label = 'non-affine', color=ec[1])
    plt.xlabel(r'$T(K)$')
    plt.ylabel('thermal expansion')
    plt.title('Cumulative thermal expansion')
    plt.ylim((-5e-4, 12e-4))
    plt.legend()
    plt.grid()
    plt.savefig('6_numeric_gruneisen/cumul_exp.png', dpi=600, bbox_inches="tight", transparent=False)

    print(cumul_exp[299], cumul_exp_a[299])
    exit()
    '''

    fig = plt.figure()

    plt.axhline(y=0, color='black')



    plt.semilogx(T,grun_affine,color=ec[0],markersize=5,label='affine Gruneisen',marker="o")
    plt.semilogx(T,grun_nonaffine,color=ec[1],markersize=5,label='non-affine Gruneisen',marker="o")

    
    
    plt.xlabel(r'$T(K)$')
    plt.ylabel(r'$\gamma$')
    plt.legend()
    plt.ylim([-15, 5])
    fig.savefig("6_numeric_gruneisen/numeric_gruneisen_small.pdf", dpi=300, bbox_inches="tight", transparent=True)

    fig.savefig("6_numeric_gruneisen/numeric_gruneisen_small.png", dpi=600, bbox_inches="tight", transparent=False)

    fig = plt.figure()

    i_cross = np.argmax(np.array(grun_nonaffine) > 0) # get crossing coordinates
    T_cross = T[i_cross]
    
    #plt.annotate(f'Crossing point: {round(T_cross)}K', (T_cross, 0.1))
    plt.plot(T[5:],grun_affine[5:],color=ec[0],markersize=5,label='affine Gruneisen',marker="o")
    plt.plot(T[5:], grun_nonaffine[5:], color=ec[1],markersize=5,label='non-affine Gruneisen',marker="o")
    plt.xlabel(r'$T(K)$')
    plt.ylabel(r'$\gamma$')
    plt.legend()
    plt.ylim([-0.5, 1.5])
    plt.grid()
    fig.savefig("6_numeric_gruneisen/numeric_gruneisen_linear_scale.png", dpi=600, bbox_inches="tight", transparent=False)

    fig = plt.figure()
    #bulk_modulus = 33.5e9 # Silica glass
    bulk_modulus = 98e9 # Silicon
    plt.title("Thermal Expansion of Si (crystalline)")

    plt.plot(T, grun_affine * heat_capacity_vol / bulk_modulus / 3 * 1e6, color=ec[0],markersize=1,label='affine theory (250 atom a-Si)',marker="o")

    plt.plot(T, grun_nonaffine * heat_capacity_vol / bulk_modulus / 3 * 1e6, color=ec[1],markersize=1,label='non-affine theory (250 atom c-Si)',marker="o")

    
    data_white = np.loadtxt('vitreosil_White1973')
    plt.plot(data_white[:, 0], data_white[:, 1], label='White (1973)', linestyle='None', marker='v', markersize=5)
    plt.ylim((-2, 4))
    

    '''
    plt.plot([110+273], [3], label='Witvrouw (1993)', linestyle='None', marker='v', markersize=5)
    plt.plot(np.array([45, 55, 65, 75, 85, 95, 105, 115, 125, 135, 145])+273, [3.6, 3.8, 3.7, 3.3, 3.2, 3.4, 2.9, None, None, 3.2, 3.0], label='Takimoto (2002)', linestyle='None', marker='o', markersize=5)
    '''


    '''
    data_middleton = np.loadtxt('middelman.dat')
    data_watanabe = np.loadtxt('watanabe.dat')
    plt.plot(data_middleton[:, 0], data_middleton[:, 1]*1e-3, label='Middelman (2015)', linestyle='None', marker='v', markersize=2)
    plt.plot(data_watanabe[:, 0], data_watanabe[:, 3], label='Watanabe (2003)', linestyle='None', marker='o', markersize=2)
    '''
    plt.xlabel(r'$T(K)$')
    plt.ylabel(r'$\alpha (K^{-1}) \times 10^{-6}$')
    plt.legend()
    #plt.ylim([-1.5, 1.5])
    #plt.ylim((-55, 8))
    plt.grid()
    fig.savefig("6_numeric_gruneisen/alpha.png", dpi=600, bbox_inches="tight", transparent=False)

    fig = plt.figure()
    plt.plot(T, specific_heat, color=ec[1],markersize=5,label='CV',marker="o")
    plt.xlabel(r'$T(K)$')
    plt.ylabel(r'$C_V (JK^{-1}g^{-1})$')
    plt.legend()
    #plt.ylim([-1.5, 1.5])
    plt.grid()
    fig.savefig("6_numeric_gruneisen/CV.png", dpi=600, bbox_inches="tight", transparent=False)

