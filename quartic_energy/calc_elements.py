import numpy as np

elements = np.load('10_quartic/elements.npy')
elements_single = np.load('10_quartic/elements2.npy')
with open('10_quartic/u0.txt') as f:
    u0 = float(f.readline())
print(u0)
print(elements.shape)
print(elements_single.shape)

n = elements.shape[0]

pet = np.zeros((n, n))

displacement = 0.5

for i in range(n):
    for j in range(n):
        if j > i:
            for k in range(2):
                for l in range(2):
                    elements[i, j, k, l] = elements[j, i, l, k]

for i in range(n):
    for j in range(n):
        pet[i, j] = (elements[i, j, 0, 0] - 2 * elements_single[j, 0] + elements[i, j, 1, 0] - 2 * elements_single[i, 0] + 4 * u0 - 2 * elements_single[i, 1] + elements[i, j, 0, 1] - 2 * elements_single[j, 1] + elements[i, j, 1, 1])/displacement**4

print(pet)
np.savetxt('10_quartic/pet.dat', pet)

## compute diagonal 2nd derivatives (check eigenvalues consistent)
## yes, these are okay (roughly)
evs = np.zeros(n)
for i in range(n):
    evs[i] = (elements_single[i, 0] - 2*u0 + elements_single[i, 1]) / displacement**2

print(evs)
np.savetxt('10_quartic/evs.dat', evs)

## compute diagonal 4th derivates (1, -4, 6, -4, 1)