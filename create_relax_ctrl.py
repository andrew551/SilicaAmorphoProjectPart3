def create_relax_ctrl(lammps_structure_file, lammps_dump_file, lammps_ctrl_file, config, minimize='no', threshold=1e-6):
    path_ACE_potential = config['path_ACE_potential']
    text = '# Test of ACE potential for SiO2 system\n'
    text += '#\n'
    text += 'units           metal\n'
    text += 'boundary        p p p\n'
    text += '#\n'
    text += 'box tilt large\n'
    text += '#\n'
    text += 'atom_style      atomic\n'
    text += '#\n'
    text += 'read_data       %s\n' % (lammps_structure_file)
    text += '#\n'
    text += 'pair_style      hybrid/overlay pace table spline 6000\n'
    text += 'pair_coeff      * * pace %s/SiO2-4_24-20-16-12.yace O Si \n' % (
        path_ACE_potential)
    text += 'pair_coeff      1 1 table %s/SiO2-4_24-20-16-12_pairpot.table O_O \n' % (
        path_ACE_potential)
    text += 'pair_coeff      1 2 table %s/SiO2-4_24-20-16-12_pairpot.table O_Si \n' % (
        path_ACE_potential)
    text += 'pair_coeff      2 2 table %s/SiO2-4_24-20-16-12_pairpot.table Si_Si \n' % (
        path_ACE_potential)
    text += '#\n'
    text += 'neighbor        2 bin\n'  # this corresponds to 4 Angstrom
    text += 'dump myDump all custom 1 %s id type x y z fx fy fz\n' % (lammps_dump_file)
    text += 'dump_modify myDump format line \"%d %d %10.8e %10.8e %10.8e %10.8e %10.8e %10.8e\"\n'
    if minimize == 'vc-relax':
        text += 'thermo_style custom temp pe ke pxx pyy pzz pxy pxz pyz\n'
        text += 'compute 1 all pressure NULL virial\n'
        # text+='fix 1 all box/relax iso 0.0 vmax 0.000001\n'
        text += 'fix 1 all box/relax x 0.0 y 0.0 z 0.0 xy 0.0 yz 0.0 xz 0.0 couple none vmax 0.005\n'
        text += f'minimize {threshold} {threshold} 1000 1000'
    if minimize == 'relax':
        text += 'thermo_style custom temp pe ke pxx pyy pzz pxy pxz pyz\n'
        text += 'compute 1 all pressure NULL virial\n'
        text += 'min_style cg\n'
        text += f'minimize {threshold} {threshold} 1000 1000'
    if minimize == 'no':
        text += 'run 0'
    #filename_ctrl = './%s' % (lammps_ctrl_file)
    filename_ctrl = lammps_ctrl_file
    tmpf = open(filename_ctrl, 'w')
    tmpf.write(text)
    tmpf.close()
    return filename_ctrl