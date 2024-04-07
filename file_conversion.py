import numpy as np
from ase.io import read, write
from collections import Counter
import numpy as np
from ase.build.tools import sort
import makeconfig
import traceback
import os
config = makeconfig.config()

def get_cell_parameters(ASE_structure):
    cell = ASE_structure.get_cell()
    a = np.linalg.norm(cell[0])
    b = np.linalg.norm(cell[1])
    c = np.linalg.norm(cell[2])
    alpha = np.arccos(np.dot(cell[1], cell[2]) / (c * b))
    gamma = np.arccos(np.dot(cell[1], cell[0]) / (a * b))
    beta = np.arccos(np.dot(cell[2], cell[0]) / (a * c))
    return a, b, c, alpha, beta, gamma

def write_LAMMPS_structure(structure, filename_lammps, style=None, supercell=(1, 1, 1), by_element=True):
    types = structure.get_atomic_numbers()

    atom_type_unique = np.unique(types, return_index=True)

    # To use unique without sorting
    sort_index = np.argsort(atom_type_unique[1])
    elements_atomic_numbers = np.array(atom_type_unique[0])[sort_index]
    elements_count = np.diff(np.append(np.array(atom_type_unique[1])[sort_index], [len(types)]))
    elements_symbol = ''
    chem_symbol = {}
    chem_symbol[1] = 'H'
    chem_symbol[8] = 'O'
    chem_symbol[6] = 'C'
    chem_symbol[14] = 'Si'
    for atomic_number in elements_atomic_numbers:
        elements_symbol += ' %s ' % (chem_symbol[atomic_number])

    if by_element:
        type_index_unique = np.unique(types, return_index=True)[1]
        # To use unique without sorting
        sort_index = np.argsort(type_index_unique)
        type_index_unique = np.array(type_index_unique)[sort_index]
        count_index_unique = np.diff(np.append(type_index_unique, [len(types)]))
        atom_index = []
        for i, index in enumerate(count_index_unique):
            atom_index += [i for j in range(index)]
    else:
        atom_index = structure.get_atom_type_index()

    atom_index_unique = np.unique(atom_index, return_index=True)[1]

    masses = structure.get_masses()
    # charges = structure.get_charges()

    positions = structure.get_positions()
    number_of_atoms = len(positions)

    lammps_data_file = 'Generated using dynaphopy\n\n'
    lammps_data_file += '{0} atoms\n\n'.format(number_of_atoms)
    lammps_data_file += '{0} atom types\n\n'.format(len(atom_index_unique))

    cell = structure.get_cell()

    # generate lammps oriented lattice vectors
    ax = np.linalg.norm(cell[0])
    bx = np.dot(cell[1], cell[0] / np.linalg.norm(cell[0]))
    by = np.linalg.norm(np.cross(cell[0] / np.linalg.norm(cell[0]), cell[1]))
    cx = np.dot(cell[2], cell[0] / np.linalg.norm(cell[0]))
    cy = (np.dot(cell[1], cell[2]) - bx * cx) / by
    cz = np.sqrt(np.dot(cell[2], cell[2]) - pow(cx, 2) - pow(cy, 2))

    # get rotation matrix (poscar->lammps) from lattice vectors matrix
    lammps_cell = np.array([[ax, bx, cx], [0, by, cy], [0, 0, cz]])
    trans_cell = np.dot(np.linalg.inv(cell), lammps_cell.T)

    # rotate positions to lammps orientation
    positions = np.dot(positions, trans_cell)

    # build lammps lattice vectors
    a, b, c, alpha, beta, gamma = get_cell_parameters(structure)
    xhi = a
    xy = b * np.cos(gamma)
    xz = c * np.cos(beta)
    yhi = np.sqrt(pow(b, 2) - pow(xy, 2))
    yz = (b * c * np.cos(alpha) - xy * xz) / yhi
    zhi = np.sqrt(pow(c, 2) - pow(xz, 2) - pow(yz, 2))

    xhi = xhi + max(0, 0, xy, xz, xy + xz)
    yhi = yhi + max(0, 0, yz)

    if xy > 0:
        xhi -= xy
    if xz > 0:
        xhi -= xz
    if yz > 0:
        yhi -= yz

    # write lammpstrj file
    lammps_data_file += '\n{0:20.10f} {1:20.10f} xlo xhi\n'.format(0, xhi)
    lammps_data_file += '{0:20.10f} {1:20.10f} ylo yhi\n'.format(0, yhi)
    lammps_data_file += '{0:20.10f} {1:20.10f} zlo zhi\n'.format(0, zhi)
    lammps_data_file += '{0:20.10f} {1:20.10f} {2:20.10f} xy xz yz\n\n'.format(xy, xz, yz)

    lammps_data_file += 'Masses\n\n'

    for i, index in enumerate(atom_index_unique):
        lammps_data_file += '{0} {1:20.10f} \n'.format(i + 1, masses[index])

    lammps_data_file += '\nAtoms\n\n'
    charges = None
    # if charges is not None:
    #    for i, row in enumerate(positions):
    #        lammps_data_file += '{0} {1} {2} {3:20.10f} {4:20.10f} {5:20.10f}\n'.format(i + 1, atom_index[i] + 1,
    #                                                                                    charges[i], row[0], row[1], row[2])
    # else:
    if style == None:
        for i, row in enumerate(positions):
            lammps_data_file += '{0} {1} {2:20.10f} {3:20.10f} {4:20.10f}\n'.format(i + 1, atom_index[i] + 1, row[0],
                                                                                    row[1], row[2])
    elif style == 'charge-dummy':
        dummy = 0.12345
        for i, row in enumerate(positions):
            lammps_data_file += '{0} {1} {2:20.10f} {3:20.10f} {4:20.10f} {5:20.10f}\n'.format(i + 1, atom_index[i] + 1, dummy, row[0],
                                                                                    row[1], row[2])
    else:
        raise Exception(f"unimplemented style {style}!")

    tmpf = open(filename_lammps, 'w')
    tmpf.write(lammps_data_file)
    tmpf.close()
    return

'''
fix the atomic numbers to be chemically correct for the material (guess using stoichometry)
'''
def fix_atomic_numbers(x):
    atom_types = x.get_atomic_numbers()
    #print(atom_types)
    counts = Counter(atom_types)
    numbers = list(counts.keys())
    if config['material'] == 'SiO2':
        if not len(numbers) == 2:
            raise Exception(f"ERROR: expected structure with 2 types of atom, got {list(numbers)}")
        
        # assign whichever atom type is more frequent to be oxygen, and the other Si
        if counts[numbers[0]] > counts[numbers[1]]:
            maptypes = {numbers[0] : 8, numbers[1] : 14}
        else:
            maptypes = {numbers[1] : 8, numbers[0] : 14}
    elif config['material'] == 'Si':
        if not len(numbers) == 1:
            raise Exception(f"ERROR: expected structure with 1 types of atom, got {list(numbers)}") 
        # assign all atoms as Si
        maptypes = {numbers[0] : 14}
    else:
        raise Exception(f"unsupported matieral: {config['material']}")    
    new_types = [maptypes[_] for _ in atom_types]
    x.set_atomic_numbers(new_types)
    x = sort(x, tags=x.get_atomic_numbers())
    return x
'''
x -> x_standardised
this function reads the files in the format 
(1) lammps format found on https://github.com/WignerTransport/AmorFo/tree/master/SiO2_structures
or (2) standard vasp format
or (3) lammps dump format
in principle, you can chain as many other different format here as you need them...
'''
def _read_any(file_in):
    if not os.path.isfile(file_in):
        raise Exception(f'ERROR file not found:"{file_in}"')
    try:
        x = read(file_in, format = 'lammps-data', style='charge')
        print('read input file as lammps-data (style=charge)')
        return x
    except Exception:
        pass
        #print("Couldn't read file as lammps-data-charge file type")
    try:
        x = read(file_in, format = 'vasp')
        print('read input file as vasp')
        return x
    except Exception:
        pass
        #print("Couldn't read file as lammps-data-charge file type")
    try:
        x = read(file_in, format = 'lammps-data', style='atomic')
        print(x)
        print('read input file as lammps-data (style=default)')
        return x
    except Exception as e:
        pass
        #print("Couldn't read file as lammps-datafile type")
        #traceback.print_exception(e.__traceback__)
    try:
        x = read(file_in, format = 'lammps-dump-txt', style = 'atomic')
        print(x)
        print('read input file as lammps-dump-txt (style=atomic)')
        return x
    except Exception:
        pass
        #print("Couldn't read file as lammps-datafile type")
    try:
        x = read(file_in, format = None)
        print('read input file using unknown type guess')
        return x
    
    except Exception:
        pass
        #print("Couldn't read file as [unknown:guess] file type")
    raise Exception("ERROR: none of methods tried could read input data")

def read_reg(file_in):
    return fix_atomic_numbers(_read_any(file_in))

def convert_and_regularize_file(file_in, file_out, out_type='lammps', style=None):
    x = read_reg(file_in)
    if out_type == 'lammps':
        write_LAMMPS_structure(x, file_out, style=style)
    elif out_type == 'vasp':
        write_VASP_structure(file_out, x)
    elif out_type == 'xyz':
        write(file_out, x, format='xyz')
    else:
        raise 
    Exception(f'unsupported output type {out_type}')



def write_VASP_structure(filename_POSCAR, structure, scaled=False, supercell=(1, 1, 1)):
    cell = structure.get_cell()

    types = structure.get_atomic_numbers()

    atom_type_unique = np.unique(types, return_index=True)

    # To use unique without sorting
    sort_index = np.argsort(atom_type_unique[1])
    elements_atomic_numbers = np.array(atom_type_unique[0])[sort_index]
    elements_count = np.diff(np.append(np.array(atom_type_unique[1])[sort_index], [len(types)]))
    elements_symbol = ''
    chem_symbol = {}
    chem_symbol[1] = 'H'
    chem_symbol[8] = 'O'
    chem_symbol[6] = 'C'
    chem_symbol[14] = 'Si'
    for atomic_number in elements_atomic_numbers:
        elements_symbol += ' %s ' % (chem_symbol[atomic_number])

    vasp_POSCAR = 'Generated using dynaphopy\n'
    vasp_POSCAR += '1.0\n'
    for row in cell:
        vasp_POSCAR += '{0:20.10f} {1:20.10f} {2:20.10f}\n'.format(*row)
    vasp_POSCAR += ' %s ' % (elements_symbol)
    vasp_POSCAR += ' \n'
    vasp_POSCAR += ' '.join([str(i) for i in elements_count])

    if scaled:
        scaled_positions = structure.get_scaled_positions()
        vasp_POSCAR += '\nDirect\n'
        for row in scaled_positions:
            vasp_POSCAR += '{0:15.15f}   {1:15.15f}   {2:15.15f}\n'.format(*row)

    else:
        positions = structure.get_positions(supercell=supercell)
        vasp_POSCAR += '\nCartesian\n'
        for row in positions:
            vasp_POSCAR += '{0:20.10f} {1:20.10f} {2:20.10f}\n'.format(*row)

    tmpf = open(filename_POSCAR, 'w')
    tmpf.write(vasp_POSCAR)
    tmpf.close()
    return

#x: ase atoms object
def get_density(x):
    v = x.get_volume()
    num_ox = len([n for n in x.get_atomic_numbers() if n == 8])
    num_si = len([n for n in x.get_atomic_numbers() if n == 14])
    m = num_ox*16 + num_si*28.0855
    return m/v/6.022e23*10e23 # Avogadros constant
'''
test conversion function
'''
if __name__ == '__main__':
    x = read_reg('/mnt/scratch2/q13camb_scratch/adps2/output_folder_chik5001/relax_output/20240320212733/relaxed_structure_starting_point.POSCAR')
    x = read_reg('/mnt/scratch2/q13camb_scratch/adps2/output_folder_chik5001/relax_output/20240321211443/relaxed_structure_starting_point.POSCAR')
    x = read_reg('/mnt/scratch2/q13camb_scratch/adps2/output_folder_chik5001/relax_output/20240321211443/relaxed_structure_starting_point.POSCAR')
    x = read_reg('/mnt/scratch2/q13camb_scratch/adps2/output_folder_chik5001/relax_output/20240326012914/relaxed_structure_starting_point.POSCAR')
    x = read_reg('/mnt/scratch2/q13camb_scratch/adps2/output_folder_chik5001/relax_output/20240326012914/steps/struct_1_relax_1.lammps_struct')
    x = read_reg('/mnt/scratch2/q13camb_scratch/adps2/output_folder_chik5001/relax_output/20240329130444/relaxed_structure_starting_point.POSCAR')
    x = read_reg('/mnt/scratch2/q13camb_scratch/adps2/input_folder2/deringher5184/POSCAR_5184 1')
    x = read_reg('/mnt/scratch2/q13camb_scratch/adps2/quenched6000K/1_relax/relax_output/20240403181127/steps/struct_1_relax_8.lammps_struct')
    print('density is', get_density(x))
    '''
    input_struct_path = '/users/asmith/grun_in/models24k/Coords.dat'
    input_struct_path = '/mnt/scratch2/q13camb_scratch/adps2/output_folder1/anneal_output/20240308185606/Coords_ACE_cg_min.dat'
    #output_struct_path = '/users/asmith/grun_in/models24k/Coords_fixed2'
    #input_struct_path = '/users/asmith/grun_in/model1536/POSCAR_1536'
    input_struct_path = '/users/asmith/grun_out/relax_output/20240311184354/steps/final_conf.dump'
    output_struct_path = '/users/asmith/grun_in/model1536/Coords_regxx_2'
    print(f'Running file conversion test with in-file/out-file = {input_struct_path} {output_struct_path}')
    convert_and_regularize_file(input_struct_path, output_struct_path)
    '''