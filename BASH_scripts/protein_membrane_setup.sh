#!/bin/bash

# Author: Afroditi-Maria Zaki
# Affiliation: SBCB, University of Oxford
# Last Updated: 21.04.2022

#################################################################################################
# This script can be used to setup a ligand-protein-lipid system. 
# It requires a .pdb file for the protein and a preformed equilibrated lipid bilayer (.gro file).
# It will concatenate the configuration files, inflate the bilayer, then gradually shrink it,
# solvate, remove unwanted water molecules from the lipid region, add ions, energy minimise and 
# equilibrate, first in a 5-ns NVT step and then in a 10-ns NPT step.
#################################################################################################


# Using the protein.pdb file, generate the protein .gro file along with the .itp and .top files
gmx pdb2gmx -f protein.pdb -o protein.gro -ignh -ter -ff amber99sb-ildn -ss -water tip3p

# Concatenate the protein and the POPC bilayer in the same box, with the bilayer box dimensions
# Read the dimensions of the protein box
xdim=$(cat popc.gro | awk -F ' ' '{print $1}' | awk 'END{print}')
ydim=$(cat popc.gro | awk -F ' ' '{print $2}' | awk 'END{print}')
zdim=$(cat protein.gro | awk -F ' ' '{print $3}' | awk 'END{print}')
echo $xdim $ydim $zdim

# Compute the centre of the box
halfxdim=$(echo "($xdim*0.5)" | bc)
halfydim=$(echo "($ydim*0.5)" | bc)
# halfzdim was determined by visualization of the protein overlayed withthe POPC bilayer 
# and by trial-and-error until the TMD was properly embedded in the membrane 
halfzdim=8.15
echo $halfxdim $halfydim $halfzdim

# Align the protein along its principal axis
echo "Protein" | gmx editconf -f protein.gro -princ -o protein_princ.gro 
# Rotate the protein so that its principal axis is perpendicular to the membrane plane
gmx editconf -f protein_princ.gro -rotate 0 90 0 -o protein_rotated.gro
# Create new box with same dimensions as POPC bilayer box, and place the protein at the center of the box in the x- and y- dimensions and at the required z- coordinate
# so that the protein TMD is embedded in the membrane
gmx editconf -f protein_rotated.gro -o protein_newbox.gro -box ${xdim} ${ydim} ${zdim} -center ${halfxdim} ${halfydim} ${halfzdim}

sed '$ d' protein_newbox.gro > protein_newbox.gro.tmp
sed '1,2d' popc.gro > popc.gro.tmp

cat protein_newbox.gro.tmp popc.gro.tmp > protein_popc.gro

lines=$(wc -l < protein_popc.gro)
lines=$((lines-3))
sed -i "2s/.*/$lines/" protein_popc.gro
sed -i "1 s/^.*/GABA RDL in POPC/" protein_popc.gro
sed -i "$ s/^.*/$xdim $ydim $zdim/" protein_popc.gro

#Renumber system molecules
#gmx editconf -f protein_popc.gro -o system.gro -resnr 41

# Make .ndx file of protein_lipid system
echo "q" | gmx make_ndx -f protein_popc.gro

# Generate strong_restrainsts .itp files - required for the membrane inflation and shrinking process
chains='A B C D E'

for i in $chains
do
	cp posre_Protein_chain_$i.itp strong_posre_chain_$i.itp
	sed -i 's/1000/100000/g' strong_posre_chain_$i.itp
	sed -i '/100000     1/ s/100000/1000/' strong_posre_chain_$i.itp
done

# Create topology file
rm topol.top
echo "; Include forcefield parameters
#include \"amber99sb-ildn.ff/forcefield.itp\"
#include \"./SLipids_2016/forcefield.ff/forcefield.itp\"

" > topol.top

for i in $chains
do
  echo "; Include chain topologies
#include \"topol_Protein_chain_$i.itp\"

; Strong position restraints for InflateGRO
#ifdef STRONG_POSRES
#include \"strong_posre_chain_$i.itp\"
#endif
  
" >> topol.top
done

echo "; Include lipid topology
#include \"POPC.itp\"

#ifdef POSRES_POPC
#include \"posre_POPC.itp\"
#endif

; Include water topology
#include \"amber99sb-ildn.ff/tip3p.itp\"

#ifdef POSRES_WATER
; Position restraint for each water oxygen
[ position_restraints ]
;  i funct       fcx        fcy        fcz
   1    1       1000       1000       1000
#endif

; Include topology for ions
#include \"amber99sb-ildn.ff/ions.itp\"

[ system ]
; Name
Protein in POPC" >> topol.top

echo "[ molecules ]
; Compound        #mols
Protein_chain_A     1
Protein_chain_B     1
Protein_chain_C     1
Protein_chain_D     1
Protein_chain_E     1
POPC              512" >> topol.top

rm \#*

# Generate system box and place protein_lipid system in the center	
gmx editconf -f protein_popc.gro -o system.gro -box $xdim $ydim  14.5 -center $halfxdim $halfydim 5
#cp protein_popc.gro system.gro

# Use script to inflate lipid membrane
perl inflategro.pl system.gro 4 POPC 14 system_inflated.gro 5 area.dat > inflate.log

# Some lipids are deleted; change number of lipids accordingly in .top file
lipids=$(grep cut-off inflate.log | awk '{print $3}')
lipids=$((512-lipids))

sed -i "s/POPC              512/POPC              $lipids/g" topol.top

# Gradually shrink lipid membrane in 24 steps
# Usually works fine, but make sure this is the optimal number of iterations for the particular system
# Area per Lipid (ApL) should be slightly higher than the exp. ApL
bash shrink.sh

# Change the z-dimension of the box to fit the protein and leave approx. 1 nm separation between the protein and the box limits on both sides.
xdim=$(cat em_shrink24.gro | awk -F ' ' '{print $1}' | awk 'END{print}')
ydim=$(cat em_shrink24.gro | awk -F ' ' '{print $2}' | awk 'END{print}')
zdim=$(cat em_shrink24.gro | awk -F ' ' '{print $3}' | awk 'END{print}')
echo $xdim
echo $ydim

xcoord=$(echo "($xdim*0.5)" | bc)
ycoord=$(echo "($ydim*0.5)" | bc)
zcoord=$(echo "($zdim*0.33)" | bc)

gmx editconf -f em_shrink24.gro -o em_shrink24_box.gro -box $xdim $ydim $zdim -center $xcoord $ycoord $zcoord

# Solvate box
gmx solvate -cp em_shrink24_box.gro -cs -p topol.top -o system_solvated.gro

# Remove water molecules that were inserted in the lipid region
perl water_deletor.pl -in system_solvated.gro -out system_solvated_fixed.gro -ref O32 -middle C316 -nwater 3 > water_deletor.log

# Update topology file with the correct number of water molecules
waters=$(head -n -4 water_deletor.log | awk '{print $1}' | awk END{print})
sed -i '$s/^/;/' topol.top
echo "SOL		$waters" >> topol.top

# Add ions to neutralise and also reach a salt concentration of 0.15 M
gmx grompp -f em.mdp -c system_solvated_fixed.gro -r system_solvated_fixed.gro -p topol.top -o ions.tpr -maxwarn 1
echo "SOL" | gmx genion -s ions.tpr -o system_solvated_ions.gro -p topol.top -pname NA -nname CL -neutral -conc 0.15

# The system is now prepared. Energy minimise it
gmx grompp -f em.mdp -c system_solvated_ions.gro -p topol.top -r system_solvated_ions.gro -o em
gmx mdrun -deffnm em -v

# Make an index file to include one new group that will only contain the protein and the lipids
gmx make_ndx -f em.gro << EOF
1 | r POPC
q
EOF

# Perform an equilibration in the NVT ensemble for 5 ns
gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -n index.ndx -o nvt
gmx mdrun -deffnm nvt -v

# Perform an equilibration in the NPT ensemble for 10 ns
gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -n index.ndx -o npt
gmx mdrun -deffnm npt -v
