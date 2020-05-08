#!/bin/bash
 
lattice_constants="3.50 3.51 3.52 3.53 3.54 3.55 3.56 3.57 3.58 3.59 3.60 3.61 3.62 3.63 3.64 3.65 3.66 3.67 3.68 3.69 3.70 3.71 3.72"
 
basis_file=BASIS_MOLOPT
potential_file=GTH_POTENTIALS
template_file=diamond_LDA.inp
input_file=diamond_lattice_constant.inp
structure_file=diamond_fd3m.xyz
 
 
for ii in $lattice_constants ; do
    work_dir=lattice_constant_${ii}Ang_LDA
    if [ ! -d $work_dir ] ; then
        mkdir $work_dir
    else
        rm -r $work_dir/*
    fi
    sed -e "s/LT_cutoff/100/g; s/RUN_TYPE ENERGY/RUN_TYPE GEO_OPT/g; s/A 3.57370926 0.0000 0.0000/A ${ii} 0.0000 0.0000/g; s/B 0.00000000 3.57370926 0.00000/B 0.00000000 ${ii} 0.00000/g; s/C 0.00000000 0.000000 3.57370926/C 0.00000000 0.000000 ${ii}/g; /MAX_SCF 1/d" \
        $template_file > $work_dir/$input_file
    cp $basis_file       $work_dir
    cp $potential_file   $work_dir
    cp $structure_file   $work_dir
done
