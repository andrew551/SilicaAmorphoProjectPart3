import numpy as np
import matplotlib.pyplot as plt
import os
import sys
sys.path.insert(1, os.path.join(sys.path[0], '..')) # "hack" to add the parent directory to path
import file_conversion
#filename = 'frequencies_elpa_new.dat'
#filename = '/mnt/scratch2/q13camb_scratch/adps2/output_folder_chik5001/pydiag_output/20240330232833/eigenvalues.dat'

def make_plot(infile, outfile, outfile2):

    '''
    with open(filename) as f:
        data = np.array([list(map(float, line.split())) for line in f])
    v = data[:, 2]   # third column: frequencies in cm^-1
    x = np.linspace(-v[-1]*1.02, v[-1]*1.02, 4000)

    '''

    try:
        atoms = file_conversion.read_reg('../relaxed_model/POSCAR')
        volume = atoms.get_volume() / 1000 # volume in nm^3
        print("got volume of", volume)
    except Exception:
        volume = 1
        print("failed to read atoms file -- assuming volume of 1")

    ev = np.loadtxt(infile)
    ev = sorted(ev)
    print(ev[:20])
    ev = ev[3:] # remove translational modes
    v = np.sqrt(np.abs(ev)*9.648e27)/18.8e10
    x = np.linspace(-np.max(v)*1.02, np.max(v)*1.02, 4000)


    sigma = 5
    conv = 1/np.sqrt(2*np.pi*sigma**2) * np.exp(-x**2/(2*sigma**2))

    histo, _ = np.histogram(v, bins=np.append(x, x[-1]*2 - x[-2])) # bins have one more element than x as they represent edges
    dos = np.convolve(conv, histo, mode = 'same') / volume / 0.0299792458 / 2 / np.pi
    print(dos)
    plt.tick_params(direction='in')
    plt.plot(x[x > 0], dos[x > 0])
    plt.xlabel('$\hbar\omega$ $(cm^{-1})$')
    plt.ylabel('$g(\omega) (THz^{-1} \cdot nm^{-3})$')
    plt.tick_params(direction='in')
    plt.savefig(outfile, dpi=600, bbox_inches="tight", transparent=False)
    plt.show()


    arr = np.array([x[x>0], dos[x > 0]]).T
    np.savetxt("dos_save.txt", arr)

    fig, ax = plt.subplots()

    sub = np.logical_and(10 < x, x < 200)
    plt.plot(x[sub], dos[sub] / x[sub]**2)
    plt.xlabel('$\hbar\omega$ $(cm^{-1})$')
    plt.ylabel('$g(\omega) / \omega^2 (THz^{-1} cm^2 \cdot nm^{-3})$')
    plt.tick_params(direction='in')
    plt.savefig(outfile2, dpi=600, bbox_inches="tight", transparent=False)
    plt.show()

if __name__ == '__main__':
    make_plot('eigenvalues.dat', 'DOS_plot.png', 'DOS_plot_over_omega_sq.png')