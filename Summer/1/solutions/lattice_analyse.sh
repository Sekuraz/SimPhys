#!/bin/bash
 
lattice_constants="3.50 3.51 3.52 3.53 3.54 3.55 3.56 3.57 3.57370926 3.58 3.59 3.60 3.61 3.62"
 
input_file=diamond_lattice_constant.inp
output_file=diamond_lattice_constant.out
plot_file=lattice_constant_data.ssv
 
rel_cutoff=60
 
echo "# Grid cutoff vs total energy" > $plot_file
echo "# Date: $(date)" >> $plot_file
echo "# PWD: $PWD" >> $plot_file
echo "# REL_CUTOFF = $rel_cutoff" >> $plot_file
echo -n "# Cutoff (Ry) | Total Energy (Ha)" >> $plot_file
grid_header=true
for ii in $lattice_constants ; do
    work_dir=lattice_constant_${ii}Ang
    total_energy=$(grep -e '^[ \t]*Total energy' $work_dir/$output_file | awk '{print $3}')
    ngrids=$(grep -e '^[ \t]*QS| Number of grid levels:' $work_dir/$output_file | \
             awk '{print $6}')
    if $grid_header ; then
        for ((igrid=1; igrid <= ngrids; igrid++)) ; do
            printf " | NG on grid %d" $igrid >> $plot_file
        done
        printf "\n" >> $plot_file
        grid_header=false
    fi
    printf "%10.2f  %15.10f" $ii $total_energy >> $plot_file
    for ((igrid=1; igrid <= ngrids; igrid++)) ; do
        grid=$(grep -e '^[ \t]*count for grid' $work_dir/$output_file | \
               awk -v igrid=$igrid '(NR == igrid){print $5}')
        printf "  %6d" $grid >> $plot_file
    done
    printf "\n" >> $plot_file
done
