import numpy as np
file_path1 = 'dmat.dat'
output_file_path = 'dmat2.dat'  # Choose a suitable file name
# Compare the ratios of corresponding elements in the arrays
i=0
natoms=1536
dyn_mat_storage=np.zeros((natoms*3*3),dtype=np.float64)

# Read data from the first file and second file
with open(file_path1, 'r') as file1:
    with open(output_file_path, 'w') as output_file:
        for i in range (natoms):
            sum=0
            count=0
            for j in range (natoms*3*3):
                line1 = file1.readline()
                data1 = float(line1.strip().split()[2])
                dyn_mat_storage[j]=data1
            
                if (not dyn_mat_storage[j] == 0):
                    sum+=dyn_mat_storage[j]
                    count+=1
            
            diff=sum/count
            for k in range(3):
                for j in range (natoms*3):
                    if (not dyn_mat_storage[k*natoms*3+j] == 0):
                        dyn_mat_storage[k*natoms*3+j]-=diff
                    formatted_dyn= "{:.15f}".format(dyn_mat_storage[k*natoms*3+j])
                    row =  str(i*3+k+1)+ ' '+str(j+1)+' ' + formatted_dyn+ '\n'
                    # Save the result into the output file
                    output_file.write(row)