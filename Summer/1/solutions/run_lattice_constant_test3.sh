#!/bin/bash

lattice_constants="3.85 3.86 3.87 3.88 3.89 3.90 3.91 3.92"

input_file=diamond_lattice_constant.inp
output_file=diamond_lattice_constant.out

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
