#!/bin/bash

lattice_constants="3.53 3.54 3.55 3.56 3.57 3.58 3.59 3.60 3.61 3.62 3.63 3.64 3.65 3.66 3.67 3.68 3.69 3.70 3.71 3.72"

input_file=diamond_lattice_constant.inp
output_file=diamond_lattice_constant.out

counter=1
max_parallel_calcs=$(expr $no_proc_to_use / $no_proc_per_calc)
for ii in $lattice_constants ; do
    work_dir=lattice_constant_${ii}Ang_LDA
    cd $work_dir
    if [ -f $output_file ] ; then
        rm $output_file
    fi
    cp2k.sopt -i $input_file -o $output_file
    cd ..
    mod_test=$(echo "$counter % $max_parallel_calcs" | bc)
    if [ $mod_test -eq 0 ] ; then
        wait
    fi
    counter=$(expr $counter + 1)
done
wait
