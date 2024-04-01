import numpy as np
import h5py
import os
import sys
import scipy.sparse as sp
from scipy.sparse._sparsetools import expandptr
from scipy.sparse.linalg import bicgstab
from ase.io import read

folder="Results"
output_file="fc3.hdf5"
structure="POSCAR"

#iterations for sum rule

files=os.listdir("./"+folder)
n_files=len(files)

atoms=read(structure,format='vasp')
natoms=len(atoms)


fc3=[]
inds=[]
for i in range(n_files):
    with h5py.File(folder+"/fc3_%d.hdf5"%(i), "r") as f:
        data=f["fc3"][()]
        ind=f["indices"][()]
        rank=f["rank"][()]
    fc3.append(data)
    inds.append(ind)
    print("File read for rank ",rank)
    sys.stdout.flush()


fc3=np.vstack(tuple(fc3))
inds=np.hstack(tuple(inds))



print("Read indices are:\n",inds)

print('Number of non-zero elements: ',fc3.shape[0])

##--------------------
##     FUNCIONS
##--------------------




data=fc3[:,6].copy()
inds=fc3[:,:6].astype(int).copy()
no_entries=data.shape[0]

natoms=int(np.max(inds[:,0]))+1

ij=inds[:,[0,2]]

ij=np.unique(ij,axis=0)
N=ij.shape[0]

print("Starting making dictionaries")

keys = dict(list(zip(zip(inds[:,0],inds[:,1],inds[:,2],inds[:,3],inds[:,4],inds[:,5]),np.arange(no_entries))))
ijkeys=dict(list(zip(zip(ij[:,0],ij[:,1]),np.arange(N))))

del ij

print("Making dictionaries finished")


def symmetrising_array(inds,keys):
    trans=np.zeros((6,inds.shape[0]),dtype=int)
    trans[5,:]=np.arange(inds.shape[0])
    for x in range(inds.shape[0]):
        trans[0,x]=keys[(inds[x,4],inds[x,5],inds[x,0],inds[x,1],inds[x,2],inds[x,3])]
        trans[1,x]=keys[(inds[x,4],inds[x,5],inds[x,2],inds[x,3],inds[x,0],inds[x,1])]
        trans[2,x]=keys[(inds[x,2],inds[x,3],inds[x,4],inds[x,5],inds[x,0],inds[x,1])]
        trans[3,x]=keys[(inds[x,0],inds[x,1],inds[x,4],inds[x,5],inds[x,2],inds[x,3])]
        trans[4,x]=keys[(inds[x,2],inds[x,3],inds[x,0],inds[x,1],inds[x,4],inds[x,5])]
    return trans


def fill_translation(ijkeys, no_entries, inds):
    Ai=np.zeros(no_entries)
    for x in range(no_entries):
        Ai[x]=ijkeys[(inds[x,0],inds[x,2])]
    Ai=Ai*27+inds[:,1]*9+inds[:,3]*3+inds[:,5]

    return Ai


def symmetrise(data, trans):
    return (data[trans[0]]+data[trans[1]]+data[trans[2]]+data[trans[3]]+data[trans[4]]+data[trans[5]])/6

trans=symmetrising_array(inds,keys)

Ai= fill_translation(ijkeys, no_entries, inds)
Aj=np.min(trans,axis=0)

A=sp.coo_matrix((np.ones(no_entries),(Ai,Aj)),shape=(27*N,no_entries)).tocsr()

AA=A @ A.transpose(copy=True)

data=symmetrise(data,trans)

dev= A @ data

print("Max deviation before sum rule")
print(np.max(dev))

#change atol and tol to get better convergence
lambdas,info=bicgstab(AA,dev,tol=1e-6,atol=1e-15)

correction=A.T @ lambdas

correction[:]=correction[Aj]

data=data-correction


ASR=np.max(A @ data)

fc3=np.hstack((inds,data.reshape(-1,1)))


##----------------------
##      SUM RULE
##----------------------


print('Mean deviation after ASR application:')
print(np.mean(np.abs(ASR)))
print('Max deviation after ASR application:')
print(np.max(ASR))
print('Maximal value in FC3:')
print(np.max(data))



##----------------------
##     FILE WRITE
##----------------------



with h5py.File(output_file, "w") as f:
    dataset = f.create_dataset(
        "fc3",
        data=fc3,
        #chunks=(24,7),  # Define chunk size
        compression="gzip",           # Use gzip compression
        compression_opts=4            # Set compression level
    )
    dataset=f.create_dataset(
        "inds",
        data=inds,
        #chunks=(24, 6),  # Define chunk size
        compression="gzip",           # Use gzip compression
        compression_opts=4            # Set compression level
    )
    dataset=f.create_dataset(
        "data",
        data=data,
    )

