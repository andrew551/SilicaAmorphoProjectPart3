import makeconfig
import sys, os
import datetime
from pathlib import Path

import file_conversion

def create_anneal_ctrl(input_struct_path, config):
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
{makeconfig.get_potential_command(config)}\
#\n\
#neighbor        2 bin\n\
dump myDump all custom 1000 {output_dir/"final_relaxation.dump"} id type x y z vx vy vz\n\
dump_modify myDump format line "%d %d %10.8e %10.8e %10.8e %10.8e %10.8e %10.8e"\n\
 \n\
timestep        0.0005\n\
thermo          10\n\
thermo_style    custom  step density temp press ke pe etotal\n\
#######################################################\n\
 \n\
fix             F1 all npt temp 1.0 1500.0 0.01 iso 0.0 0.0 10.0\n\
run             20000\n\
unfix           F1\n\
write_data      {output_dir/"1_heated.dat"}\n\
 \n\
fix             F2 all npt temp 1500.0 1500.0 0.01 iso 0.0 0.0 5.0\n\
run             100000\n\
unfix           F2\n\
write_data      {output_dir/"2_equilibrated.dat"}\n\
 \n\
fix             F5 all npt temp 1500.0 1.0 0.01 iso 0.0 0.0 5.0\n\
run             300000\n\
unfix           F5\n\
write_data      {output_dir/"3_results.dat"}\n\
 \n\
min_style       cg\n\
min_modify      line forcezero\n\
minimize        1.0e-18 1.0e-18 30000 30000\n\
write_data      {output_dir/"Coords_ACE_cg_min.dat"}\n'
    return command
