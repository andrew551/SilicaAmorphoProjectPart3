import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import beta


with open('3_diagonalise/eigenmodes.bin', 'rb') as f:
    arr1 = np.frombuffer(f.read())
    #first_bytes = [f.read(1) for i in range(100)]
    #print(first_bytes)

ev = np.loadtxt('3_diagonalise/frequencies.dat')
freq = ev[:, 2]#sorted(ev)

print(arr1.shape)
n = round(arr1.size**0.5)
arr1 = arr1.reshape((n, n))
print(arr1.shape)
#arr1 = arr1.T
pr = np.linalg.norm(arr1, axis=1, ord=4)
np.savetxt('5_gruneisen/pr.dat', np.c_[freq, pr])

ev3 = arr1.reshape((n, n//3, 3))
ev3_norms = np.linalg.norm(ev3, axis=2)
print(ev3_norms.shape)

pr3 = np.linalg.norm(ev3_norms, axis=1, ord=4)
np.savetxt('5_gruneisen/pr3.dat', np.c_[freq, pr3])

plt.scatter(freq, pr,s=2)
plt.title("PR ratio")
plt.savefig("pr/PR_all.png")
plt.clf()
plt.scatter(freq, pr3,s=2)
plt.title("PR ratio (3-vector-norm)")
plt.savefig("pr/PR_all3.png")

plt.clf()
#inds = [0, 5, 8, 10, 15, 20, 100, 1000, 2000, 7000, 10000, 15200]
inds = [0, 5, 500, 1000, 2000, 10000, 15200, n-50, n-10, n-1]
for ind in inds:
    if ind >= ev3_norms.shape[0]:
        continue
    mags = ev3_norms[ind, :]
    m_sorted = sorted(mags, reverse=True)
    plt.semilogx(np.arange(len(m_sorted)) + 1, m_sorted, label = f'{ind}')
plt.legend()
plt.savefig("pr/PR_distributions_seimlogx.png")

plt.clf()

for ind in inds:
    if ind >= ev3_norms.shape[0]:
        continue
    mags = ev3_norms[ind, :]
    plt.clf()
    plt.hist(mags*np.sqrt(n//3), bins=100, density=True)
    #a, b, xx, yy = beta.fit(mags**2)
    #xxx = np.linspace(beta.ppf(0.01, a, b), beta.ppf(0.99, a, b), 100)
    #plt.plot(xxx, beta.pdf(xxx, a, b))
    plt.savefig(f'pr/pdf{ind}.png')
#with open('3_diagonalise_fort/eigenmodes2.bin', 'rb') as f:
#    arr2 = np.frombuffer(f.read())
    #first_bytes_f = [f.read(1) for i in range(100)]
    #print(first_bytes_f)

#rint([int.from_bytes(a, "big")-int.from_bytes(b, "big") for a, b in zip(first_bytes, first_bytes_f)])
print(np.linalg.norm(ev3_norms, axis=1))
max_val = np.max(ev3_norms, axis=1)
part = np.partition(ev3_norms, -4, axis=1)
top4 = np.sort(part[:, -4:], axis=1)
print(np.linalg.norm(top4, axis=1))

plt.clf()
plt.scatter(freq, top4[:, -1],s=2)
plt.scatter(freq, top4[:, -2],s=2)
plt.scatter(freq, top4[:, -3],s=2)

plt.savefig("pr/maximum_coeff.png")
plt.clf()
plt.scatter(freq, top4[:, -1] / top4[:, -2], s=2)
plt.savefig("pr/first_second_ratio.png")
