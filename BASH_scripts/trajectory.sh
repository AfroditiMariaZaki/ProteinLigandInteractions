#!/bin/bash

#Processing of the trajectory to fix problems due to PBC and to remove translational and rotational motion

echo "1 0" | gmx trjconv -f traj_comp.xtc -s topol.tpr -pbc mol -ur compact -center -o traj_center.xtc
echo "1 0" | gmx trjconv -f traj_center.xtc -s topol.tpr -fit rot+trans -o traj_fit.xtc
rm traj_center.xtc

