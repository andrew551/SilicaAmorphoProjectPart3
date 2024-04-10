import numpy as np
#with open('dmataaa.dat') as f
def convert():
       x = np.loadtxt('d3mataaa.dat')
       print(x.shape)
       print(x[:10, :])
       with open("d3mataaa_ind.bin", 'wb') as f:
              for i in range(x.shape[0]):
                     f.write(x[i, :6].astype(np.int32).tobytes(order='F')) # write file in fortran format
       with open("d3mataaa_val.bin", 'wb') as f:
              f.write(x[:, 6].tobytes(order='F')) # write file in fortran format

if __name__ == '__main__':
       convert()