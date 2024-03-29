import scipy.sparse
import h5py
import numpy as np
from ase.io.vasp import read_vasp

# Define the dimensions of the matrix
atoms = read_vasp('POSCAR')
natoms=len(atoms)
dim_a = natoms
dim_b = natoms
dim_x = 3
dim_y = 3

# Load the sparse.npz matrix
loaded_data = scipy.sparse.load_npz('sparse.npz')

matrix = np.zeros((dim_a,dim_b,dim_x,dim_y),dtype=np.float64)
# Open the .dmat file for writing
for row in range(loaded_data.shape[0]):
    # Get the indices and data for the current row
    row_indices = loaded_data.indices[loaded_data.indptr[row]:loaded_data.indptr[row + 1]]
    row_data = loaded_data.data[loaded_data.indptr[row]:loaded_data.indptr[row + 1]]
    i=int(row/3)
    x=row%3
    # Write "row col value" per row
    for col, value in zip(row_indices, row_data):
        j=int(col/3)
        y=col%3
        matrix[i,j,x,y]=value
        matrix[j,i,y,x]=value

# Check sum rule
# Sum over the first dimension
summed_matrix = np.sum(matrix, axis=0)

# Iterate over the remaining indices and save to a DMAT file
test_file = "sum_rule.dmat"

with open(test_file, "w") as f:
    for j in range(natoms):
        for x in range(3):
            for y in range(3):
                f.write(str(summed_matrix[j, x, y]) + "\n")

output_file = 'original.hdf5'

# Create an HDF5 file in write mode
with h5py.File(output_file, "w") as f:
    # Create a dataset with chunking, compression, and filters
    dataset = f.create_dataset(
        "fc2",
        data=matrix,
        chunks=(24, 24, 3, 3),  # Define chunk size
        compression="gzip",           # Use gzip compression
        compression_opts=4            # Set compression level
    )

print(f"Matrix saved in {output_file}")
