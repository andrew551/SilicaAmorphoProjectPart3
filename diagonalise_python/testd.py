import numpy as np

ev = np.loadtxt('/mnt/scratch2/q13camb_scratch/adps2/output_folder_chik5001/pydiag_output/20240330191530/eigenvalues.dat')
print(ev.shape)
eg = np.loadtxt('/mnt/scratch2/q13camb_scratch/adps2/output_folder_chik5001/pydiag_output/20240330191530/eigenvectors.dat')

print(ev[3])
print(eg[:100, 3])