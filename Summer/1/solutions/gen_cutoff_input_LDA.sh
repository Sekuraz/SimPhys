#!/bin/bash
 
cutoffs="50 100 150 200 250 300 350 400 450"
 
basis_file=BASIS_MOLOPT
potential_file=GTH_POTENTIALS
template_file=diamond_LDA.inp
input_file=diamond_cutoff.inp
structure_file=diamond_fd3m.xyz
 
 
for ii in $cutoffs ; do
    work_dir=cutoff_${ii}Ry_LDA
    if [ ! -d $work_dir ] ; then
        mkdir $work_dir
    else
        rm -r $work_dir/*
    fi
    sed -e "s/LT_cutoff/${ii}/g" \
        $template_file > $work_dir/$input_file
    cp $basis_file       $work_dir
    cp $potential_file   $work_dir
    cp $structure_file   $work_dir
done
