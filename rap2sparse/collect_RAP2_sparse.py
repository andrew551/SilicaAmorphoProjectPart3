import numpy as np
import h5py
import os
import sys
import scipy.sparse as sp
from scipy.sparse.linalg import bicgstab

#######################
#
#  Collecting partial results of FC2 matrix and calculating acoustic sum rule
#
#######################

folder="Results"
output_file="fc2.hdf5"

files=os.listdir("./"+folder)
n_files=len(files)

data=[]
inds=[]

for i in range(n_files):
    with h5py.File(folder+"/fc2_%d.hdf5"%(i), "r") as f:
        data.append(f["fc2"][()])
        inds.append(f["indices"][()])
        rank=f["rank"][()]
    print("File read for rank ",rank, data[-1].shape, inds[-1].shape, inds[-1])
    sys.stdout.flush()
fc2=np.vstack(data)
inds=np.concatenate(inds) # vstack -> concatenate bug fix
del data

print("Read indices are:\n",inds)

print(f"fc2 shape: {fc2.shape}")
print(f"inds shape: {inds.shape}")

if len(fc2.shape) == 4: # dense case
    print("warning: converting dense matrix output to sparse and back ...")
    fc2_sp = np.zeros((fc2.size, 5))
    fc2_sp[:, :4] = np.array([tuple(x) for x in np.ndindex(fc2.shape)])
    fc2_sp[:, 4] = fc2.flatten()
    fc2 = fc2_sp


natoms=int(fc2[:,0].max())+1
print(natoms)
print(fc2.max(axis=0))


##--------------------
##  SYMMETRIZE
##--------------------

def SM_direct_fromsparse(fc2_sparse):

    
    i=(fc2[:,0]*9+fc2[:,2]*3+fc2[:,3]).astype(np.int64)
    j=(fc2[:,1]*9+fc2[:,3]*3+fc2[:,2]).astype(np.int64)

    print(j.max()*9,i.max()*9,natoms*9)

    coo=sp.coo_matrix((fc2_sparse[:,4],(i,j)),shape=(9*natoms,9*natoms))

    fc2_sp=sp.csr_matrix(coo,copy=False)

    print(type(fc2_sp))

    fc2_sp=(fc2_sp+fc2_sp.transpose())/2

    dev=fc2_sp.sum(axis=1)

    print('Initial maximum and mean deviation:')
    print(dev.max())
    print(dev.mean())

    #create matrix of 1/2 at the same entries as fc2_sp
    halves=fc2_sp.copy()
    halves.data[:]=1/2

    #calculate the diagonals
    diagon=np.array(halves.sum(axis=1)).reshape(1,-1)

    #create diagonals
    diagon=sp.spdiags(diagon,m=9*natoms,n=9*natoms, diags=[0], format='csr').tocsr()

    A=halves+diagon

    print(type(A),type(halves))

    del diagon

    lambdas,info=bicgstab(A,-dev)
    print('BicGStab info:')
    print(info)


    lambdai=halves @ sp.spdiags(lambdas,m=9*natoms,n=9*natoms, diags=[0], format='csr').tocsr()
    print('Max correction')
    print(np.max(lambdas))

    lambdai=lambdai+lambdai.transpose()

    print(type(fc2_sp),type(lambdai))

    fc2_sp=fc2_sp+lambdai
    
    print(type(fc2_sp),type(lambdai))

    print('Max and mean deviation from ASR after correction')
    print(fc2_sp.sum(axis=1).max())
    print(fc2_sp.sum(axis=1).mean())


    fc2_sp=fc2_sp.tocoo()
    i=fc2_sp.row
    j=fc2_sp.col
    data=fc2_sp.data

    id=i//9
    jd=j//9
    ad=(i%9)//3
    bd=i%3

    return np.stack((id,jd,ad,bd,data),axis=1)



#Current code calculates all forces, not only the upper triangular
#fc2=(fc2+np.transpose(fc2,axes=(1,0,3,2)))/2 

def writeout_dense_info(fc2):

    print('Number of non-zero elements: ',np.sum(fc2!=0))

    ##----------------------
    ##      SUM RULE
    ##----------------------

    ASR=np.sum(fc2,axis=1)
    print('Mean and max deviation from ASR:')
    print(np.mean(np.abs(ASR)))
    print(np.max(np.abs(ASR)))
    print("ASR mean ratio to mean maximum of FC2:")
    print(np.mean(np.abs(ASR))/np.mean(np.max(fc2,axis=1)))
    print("Max of FC2")
    print(np.max(np.abs(fc2)))
    print("Mean sum over rows:")
    print(np.mean((ASR)))


##----------------------
##     FILE WRITE
##----------------------

fc2=SM_direct_fromsparse(fc2)



with h5py.File(output_file, "w") as f:
    dataset = f.create_dataset(
        "fc2",
        data=fc2,
        #chunks=(16, 24, 3, 3),  # Define chunk size
        compression="gzip",           # Use gzip compression
        compression_opts=4            # Set compression level
    )

