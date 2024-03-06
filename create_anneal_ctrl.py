import makeconfig
import sys, os
import datetime
from pathlib import Path
from ase.io import read, write
from collections import Counter
import spglib
import numpy as np
from ase.build.tools import sort

def create_anneal_ctrl(config, input_struct_path):
    output_dir = config['output_dir']
    command = f'\
# Test of ACE potential for SiO2 system\n\
#\n\
units           metal\n\
boundary        p p p\n\
#\n\
box tilt large\n\
#\n\
atom_style      atomic\n\
#\n\
read_data       {input_struct_path}\n\
#\n\
pair_style      hybrid/overlay pace table spline 6000\n\
pair_coeff      * * pace {config["path_ACE_potential"]}/SiO2-4_24-20-16-12.yace O Si\n\
pair_coeff      1 1 table {config["path_ACE_potential"]}/SiO2-4_24-20-16-12_pairpot.table O_O \n\
pair_coeff      1 2 table {config["path_ACE_potential"]}/SiO2-4_24-20-16-12_pairpot.table O_Si\n\
pair_coeff      2 2 table {config["path_ACE_potential"]}/SiO2-4_24-20-16-12_pairpot.table Si_Si\n\
#\n\
#neighbor        2 bin\n\
dump myDump all custom 1000 {output_dir/"final_relaxation.dump"} id type x y z vx vy vz\n\
dump_modify myDump format line "%d %d %10.8e %10.8e %10.8e %10.8e %10.8e %10.8e"\n\
 \n\
timestep        0.001\n\
thermo          10\n\
thermo_style    custom  step density temp press ke pe etotal\n\
#######################################################\n\
 \n\
fix             F1 all npt temp 1.0 1500.0 0.05 iso 0.0 0.0 1000.0\n\
run             10000\n\
unfix           F1\n\
write_data      {output_dir/"1_heated.dat"}\n\
 \n\
fix             F2 all npt temp 1500.0 1500.0 0.05 iso 0.0 0.0 10.0\n\
run             50000\n\
unfix           F2\n\
write_data      {output_dir/"2_equilibrated.dat"}\n\
 \n\
fix             F5 all npt temp 1500.0 1.0 0.05 iso 0.0 0.0 10.0\n\
run             150000\n\
unfix           F5\n\
write_data      {output_dir/"3_results.dat"}\n\
 \n\
min_style       cg\n\
min_modify      line forcezero\n\
minimize        1.0e-18 1.0e-18 30000 30000\n\
write_data      {output_dir/"Coords_ACE_cg_min.dat"}\n'
    return command


def get_cell_parameters(ASE_structure):
    cell = ASE_structure.get_cell()
    a = np.linalg.norm(cell[0])
    b = np.linalg.norm(cell[1])
    c = np.linalg.norm(cell[2])
    alpha = np.arccos(np.dot(cell[1], cell[2]) / (c * b))
    gamma = np.arccos(np.dot(cell[1], cell[0]) / (a * b))
    beta = np.arccos(np.dot(cell[2], cell[0]) / (a * c))
    return a, b, c, alpha, beta, gamma

def write_LAMMPS_structure(structure, filename_lammps, supercell=(1, 1, 1), by_element=True):
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
    for i, row in enumerate(positions):
        lammps_data_file += '{0} {1} {2:20.10f} {3:20.10f} {4:20.10f}\n'.format(i + 1, atom_index[i] + 1, row[0],
                                                                                row[1], row[2])

    tmpf = open(filename_lammps, 'w')
    tmpf.write(lammps_data_file)
    tmpf.close()
    return


if __name__ == '__main__':
    config = makeconfig.config()
    input_struct_path = '/users/asmith/grun_in/models24k/Coords.dat'
    x = read(input_struct_path, format = 'lammps-data', style='charge')
    print(x)
    atom_types = x.get_atomic_numbers()
    print(atom_types)
    
    counts = Counter(atom_types)
    if counts[1] > counts[2]:
        maptypes = {1 : 8, 2 : 14}
    else:
        maptypes = {2 : 8, 1 : 14}
    new_types = [maptypes[_] for _ in atom_types]
    x.set_atomic_numbers(new_types)
    x = sort(x, tags=x.get_atomic_numbers())
    output_struct_path = '/users/asmith/grun_in/models24k/Coords_fixed'
    #write(output_struct_path, x)
    write_LAMMPS_structure(x, output_struct_path)
    #help(lammps)
    #x = lammps.read_data(input_struct_path)
    '''
    config['starttime'] = datetime.datetime.now().strftime("%Y%m%d%H%M%S")
    config['output_dir'] = config['output_dir_base'] / 'anneal_output' / config['starttime'] # the initial slash cannot be skipped!
    print('output_dir:', config['output_dir'])
    os.makedirs(config['output_dir'])
    command = create_anneal_ctrl(config, input_struct_path)
    print(command)
    with open(config['output_dir'] / 'anneal1500K_silica.lammps_ctrl', 'w') as f:
        f.write(command)
    '''