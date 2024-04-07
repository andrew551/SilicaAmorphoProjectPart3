import os, sys
import file_conversion
from gruneisen_prepare import parse3d_RAP, convert_aaa
import json
from pathlib import Path
import shutil

if __name__ == '__main__':
    prepare_type = sys.argv[1].lower()
    with open("config.json") as f:
        config = json.load(f)
    print(config)
    outpath = Path(config['path_output'])
    if prepare_type == 'anneal':
        pass
    elif prepare_type == 'relax':
        pass
    elif prepare_type == 'fc2_rap':
        pass
    elif prepare_type == 'diagonalise':
        pass
    elif prepare_type == 'fc3_rap':
        pass
    elif prepare_type == 'gruneisen':
        pathd3 = Path('5_gruneisen/fc3/d3mat.dat').resolve()
        if 1:
            if os.path.isdir("5_gruneisen"):
                raise Exception("already exists folder gruneseisen -- delete or move the old folder first")
            os.mkdir("5_gruneisen")
            os.mkdir("5_gruneisen/fc3")
            # make xyz file
            #file_conversion.convert_and_regularize_file("relaxed_model/POSCAR", "5_gruneisen/coords.xyz", out_type="xyz")
            # convert RAP3
            parse3d_RAP.parse3d_RAP(pathd3, 'relaxed_model/POSCAR', 'force_constants/fc3.hdf5')
            # convert to "aaa", "aac"
            os.chdir('5_gruneisen/fc3')
            os.system(f'{config["path_amorpho"]}/bin/gruneisen/convert_d3.exe {pathd3}')
            convert_aaa.convert() #convert aaa.dat -> bin(ind, val)
            os.chdir('..')
            os.chdir('..')
        # make "inputs" file
        atoms = file_conversion.read_reg("relaxed_model/POSCAR")
        natoms = len(atoms)
        file_conversion.convert_and_regularize_file("relaxed_model/POSCAR", "5_gruneisen/coords.dat", out_type="lammps", style="charge-dummy")
        with open('5_gruneisen/inputs', 'w') as f:
            f.write("MAIN                    #  possible values IS or MAIN\n")
            f.write(f"'{Path('5_gruneisen/coords.dat').resolve()}'     #  path to file with the structural model\n")  
            f.write(f"'{Path('3_diagonalise/frequencies.dat').resolve()}'     #  path to file with eigenvalues\n")
            f.write(f"'{Path('3_diagonalise/eigenmodes.bin').resolve()}'     #  path to file with eigenvectors\n")  
            f.write(f"'{Path('5_gruneisen/fc3/').resolve()}/'     #  path to files with FC3 elements\n")
        
        # copy frequencies.dat (for postproc)
        shutil.copyfile('3_diagonalise/frequencies.dat', '5_gruneisen/frequencies.dat')
        # TODO move nlines from fc3/ to gruneisen
        # fix nlines to add back x -> 27*natoms x
        with open('5_gruneisen/fc3/nlines') as f:
            n_fc3_lines = int(f.readlines()[0].strip().split()[0])
        with open('5_gruneisen/nlines', 'w') as f:
            f.write(f'{natoms*27} {n_fc3_lines}\n') 

    else:
        raise Exception(f"bad prepare type {prepare_type}")