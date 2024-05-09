import numpy as np
import h5py

def convertRAP3_new(fin_path, fout_ind, fout_val):
    fin=h5py.File(fin_path, 'r')
    force_constants=fin['fc3'][()]
    fin.close()

    print('Reading finished, Writing starting.')

    force_constants[:,:6]+=1 # add one to zero-based indeces
    
    with open(fout_ind, 'wb') as f:
        for i in range(force_constants.shape[0]):
            f.write(force_constants[i, :6].astype(np.int32).tobytes(order='F')) # write file in fortran format, row by row
    with open(fout_val, 'wb') as f:
        f.write(force_constants[:, 6].tobytes(order='F')) # write file in fortran format