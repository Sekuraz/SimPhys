#!/bin/bash

lattice_constants="3.50 3.51 3.52 3.53 3.54 3.55 3.56 3.57 3.57370926 3.58 3.59 3.60 3.61 3.62"

input_file=diamond_cutoff.inp
output_file=diamond_cutoff.out

counter=1
max_parallel_calcs=$(expr $no_proc_to_use / $no_proc_per_calc)
for ii in $lattice_constants ; do
    work_dir=lattice_constant_${ii}Ang
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
