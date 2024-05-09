import os
import sys
sys.path.insert(1, os.path.join(sys.path[0], '..')) # "hack" to add the parent directory to path
import makeconfig # this line needs to be after the sys.path.insert line
config = makeconfig.config()
import numpy as np
import scipy
import datetime
from pathlib import Path
import file_conversion
import json
import datetime
import create_relax_ctrl

def get_filename(vs):
    return f"POSCAR{(1+vs):.3f}"

def relax_all_structures(volume_strains):
    thresholds = [1e-9,1e-10,1e-11,1e-12,1e-12,1e-12,1e-13,1e-13,1e-13,1e-14,1e-14,1e-14,1e-15,1e-15, 8e-16, 5e-16, 3e-16, 1e-16]
    #thresholds = config['[6_numeric_gruneisen]thresholds']
    number_relaxations = len(thresholds)
    print(f'number relaxations relax={number_relaxations}')
    for number_structure, vs in enumerate(volume_strains):
        input_struct_path = '6_numeric_gruneisen/affine/' + get_filename(vs)
        config['output_dir'] = Path(input_struct_path+'RELAX/').resolve()
        print(f"output: {config['output_dir']}")
        os.makedirs(config['output_dir'] / 'steps', exist_ok=True)
        lammps_structure_file = config['output_dir'] / 'steps' / ('struct_%d_relax_%d.lammps_struct' % (number_structure, 1))
        file_conversion.convert_and_regularize_file(input_struct_path, lammps_structure_file)
        initial_struct_atoms = file_conversion.read_reg(input_struct_path)
        num_atoms_in_primitive_cell = len(initial_struct_atoms.get_atomic_numbers())

        
        temp_dump = config['output_dir'] / 'steps' / 'final_conf.dump'
        for n_relax in range(1, number_relaxations+1):
            lammps_dump_file=config['output_dir'] / 'steps' / ('struct_%d_relax_%d.dump'%(number_structure,n_relax))
            lammps_ctrl_file=config['output_dir'] / 'steps' / ('struct_%d_relax_%d.lammps_ctrl'%(number_structure,n_relax))
            lammps_structure_file=config['output_dir'] / 'steps' / ('struct_%d_relax_%d.lammps_struct'%(number_structure, n_relax))
            lammps_out_file=config['output_dir'] / 'steps' / ('struct_%d_relax_%d.out'%(number_structure, n_relax))
            # do only non volume-change relax
            relax_type = 'relax'

            filename_ctrl = create_relax_ctrl.create_relax_ctrl(lammps_structure_file, lammps_dump_file, lammps_ctrl_file, config,
                                                minimize=relax_type, threshold=thresholds[n_relax-1], keep_cuboidal=True)
            
            cmd = f"srun -n {config['NTASKS']} {config['path_lammps']} -in {filename_ctrl} > {lammps_out_file}" 
            print(cmd, flush=True)
            os.system(cmd)

            print('n_relax=', n_relax, ' type=%s' % (relax_type))
            next_lammps_structure_file=config['output_dir'] / 'steps' / ('struct_%d_relax_%d.lammps_struct'%(number_structure, n_relax+1))
            cmd = f'tail -n {num_atoms_in_primitive_cell+9} {lammps_dump_file} > {temp_dump}'
            print(cmd, flush=True)
            os.system(cmd)
            file_conversion.convert_and_regularize_file(temp_dump, next_lammps_structure_file)
        #convert final one to vasp
        file_conversion.convert_and_regularize_file(next_lammps_structure_file, f"6_numeric_gruneisen/nonaffine/POSCAR{(1+vs):.3f}", out_type='vasp')

if __name__ == '__main__':
    print(config)
    input_file = str(Path('relaxed_model/POSCAR').resolve())
    os.makedirs("6_numeric_gruneisen/affine", exist_ok=True)
    os.makedirs("6_numeric_gruneisen/nonaffine", exist_ok=True)
    print(input_file)
    atoms = file_conversion.read_reg(input_file)

    ## apply volume strains
    #volume_strains = np.linspace(-3e-3, 3e-3, 7)
    volume_strains = config["[6_numeric_gruneisen]volume_strains"]
    original_cell = atoms.get_cell()
    for vs in volume_strains:
        scale_factor = (1+vs)**(1/3)
        atoms.set_cell(original_cell * scale_factor, scale_atoms=True) # apply an affine distortion
        file_conversion.write_file(atoms, '6_numeric_gruneisen/affine/' + get_filename(vs), out_type='vasp')
    
    relax_all_structures(volume_strains)