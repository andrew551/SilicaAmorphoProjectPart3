import numpy as np

if __name__ == '__main__':
    input_file =  '/mnt/scratch2/q13camb_scratch/adps2/output_folderDering5184/pydiag_output/20240330232833/eigenvalues.dat'
    ev = np.loadtxt(input_file)
    freq = np.sqrt(np.abs(ev)*9.648e27)/18.8e10
    i = np.arange(ev.shape[0]) + 1
    output = np.c_[i, ev, freq]
    np.savetxt('frequencies.dat', output, fmt='%i %20.15f %20.15f')

    em = np.loadtxt('/mnt/scratch2/q13camb_scratch/adps2/output_folderDering5184/pydiag_output/20240330232833/eigenvectors.dat')
    em = em.T # transpose to save in column-order
    with open("eigenmodes.bin", 'wb') as f:
       f.write(em.tobytes(order='F')) # write file in fortran format