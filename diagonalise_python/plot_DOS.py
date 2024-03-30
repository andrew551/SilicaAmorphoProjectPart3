import numpy as np
import matplotlib.pyplot as plt

#filename = 'frequencies_elpa_new.dat'
filename = '/mnt/scratch2/q13camb_scratch/adps2/output_folder_chik5001/pydiag_output/20240330033332/eigenvalues.dat'

'''
with open(filename) as f:
    data = np.array([list(map(float, line.split())) for line in f])
v = data[:, 2]   # third column: frequencies in cm^-1
x = np.linspace(-v[-1]*1.02, v[-1]*1.02, 4000)

'''

ev = np.loadtxt(filename)
ev = sorted(ev)
print(ev[:20])
v = np.sqrt(np.abs(ev)*9.648e27)/18.8e10
x = np.linspace(-np.max(v)*1.02, np.max(v)*1.02, 4000)


sigma = 5
conv = 1/np.sqrt(2*np.pi*sigma**2) * np.exp(-x**2/(2*sigma**2))

histo, _ = np.histogram(v, bins=np.append(x, x[-1]*2 - x[-2])) # bins have one more element than x as they represent edges
dos = np.convolve(conv, histo, mode = 'same')

plt.plot(x[x > 0], dos[x > 0])
plt.xlabel('$\hbar\omega$ $(cm^{-1})$')
plt.ylabel('DOS')
plt.savefig('DOS_plot.png', dpi=400)
plt.show()
