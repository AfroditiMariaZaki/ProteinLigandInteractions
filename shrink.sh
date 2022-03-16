#!/bin/bash

# This script assumes you have already inflated the system with the InflateGRO command
# provided in the tutorial, and that you have further updated topol.top correctly to
# reflect the number of POPC lipids that were deleted.
#
# This script assumes that the gmx binary is in the $PATH and that all necessary input
# files (coordinates, topology, .mdp) are in the current working directory.
#
# Confirm that you have achieved an appropriate area per lipid by inspecting the area_shrink*.dat
# files along the shrinking iterations.

# prevent accumulation of backup files like mdout.mdp
export GMX_MAXBACKUP=-1

# Run energy minimization on the inflated system
gmx grompp -f em_inflate.mdp -c system_inflated.gro -p topol.top -r system_inflated.gro -o em_inflated.tpr -maxwarn 1
gmx mdrun -deffnm em_inflated -v

# make molecules whole
echo 0 | gmx trjconv -s em_inflated.tpr -f em_inflated.gro -o tmp.gro -pbc mol
mv tmp.gro em_inflated.gro

# loop over 25 shrinking iterations
for curr in {1..24} 
do
    echo "########################################"
    echo "#"
    echo "# RUNNING SHRINKING ITERATION ${curr}..."
    echo "#"
    echo "########################################"

    prev=$((curr - 1))
    if [ $curr -eq 1 ]; then
        if [ ! -e em_inflated.gro ]; then
            echo "em_inflated.gro does not exist! Exiting."
            exit;
        fi
        # special file name if doing the first iteration
        perl inflategro.pl em_inflated.gro 0.95 POPC 0 system_shrink${curr}.gro 5 area_shrink${curr}.dat
    else
        if [ ! -e em_shrink${prev}.gro ]; then
            echo "em_shrink${prev}.gro does not exist! Exiting."
            exit;
        fi
        # otherwise use minimized coordinates from previous iteration
    perl inflategro.pl em_shrink${prev}.gro 0.95 POPC 0 system_shrink${curr}.gro 5 area_shrink${curr}.dat
    fi

    # run grompp and mdrun to carry out energy minimization
    gmx grompp -f em_shrink.mdp -c system_shrink${curr}.gro -r system_shrink${curr}.gro -p topol.top -o em_shrink${curr}.tpr -maxwarn 1
    gmx mdrun -deffnm em_shrink${curr} -v

    # make molecules whole
    echo "0" | gmx trjconv -s em_shrink${curr}.tpr -f em_shrink${curr}.gro -o tmp.gro -pbc mol
    mv tmp.gro em_shrink${curr}.gro

done

exit;
