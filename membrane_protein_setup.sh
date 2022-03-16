#!/bin/bash


##########################################
# This script can be used to setup a protein-lipid system. 
# It requires a .pdb file for the protein and a preformed equilibrated lipid bilayer (.gro file).
# It will concatenate the configuration files, inflate the bilayer, then gradually shrink it, solvate, 
# remove unwanted water molecules from the lipid region, add ions, energy minimise and equilibrate, first in an NVT step and then in an NPT step.

# Convert .pdb to .gro and select the AMBER99SB-ILDN force field and the TIP3P water model
gmx pdb2gmx -f protein.pdb -o protein.gro -ignh -ter -chainsep id -ff amber99sb-ildn -water tip3p

# Generate the unit cell with the optimal dimensions, rotate protein and translate protein so that it lies in parallel to the z-axis and the TMD domain is embedded in the POPC membrane
# The correct translation of the protein so that its TMD is in the POPC membrane was found manually, by tril-and-error. Surely there is a more clever way to do this.
echo "Protein" | gmx editconf -f protein.gro -o protein_princ.gro -princ -box 12.67761  13.30903 10.00000 -c
gmx editconf -f protein_princ.gro -o protein_aligned.gro -rotate 0 90 0
gmx editconf -f protein_aligned.gro  -o protein_aligned_newbox.gro -box 12.67761  13.30903 10.00000 -c
gmx editconf -f protein_aligned_newbox.gro -o protein_aligned_newbox_translated.gro -translate 0 0 3.8

# Modify and concatenate the protein .gro file and the POPC .gro file
sed '1,2d' POPC.gro > POPC.gro.tmp
tac protein_aligned_newbox_translated.gro | sed "1d" | tac > protein_aligned_newbox_translated.gro.tmp
cat protein_aligned_newbox_translated.gro.tmp POPC.gro.tmp > protein_popc.gro
rm POPC.gro.tmp protein_aligned_newbox_translated.gro.tmp

lines=$(wc -l < protein_popc.gro)
lines=$((lines-3))
sed -i "2s/.*/$lines/" protein_popc.gro

sed -i "1 s/^.*/GABA RDL in POPC/" protein_popc.gro

gmx editconf -f protein_popc.gro -o protein_popc_renumbered.gro -resnr 1

# Make .ndx file of protein_lipid system
echo "q" | gmx make_ndx -f protein_popc_renumbered.gro

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
gmx editconf -f protein_popc_renumbered.gro -o system.gro -box 12.67761  13.30903  14.5 -center 6.338805 6.654515 5

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
halfxdim=$(($xdim/2))
halfydim=$(($ydim/2))
gmx editconf -f em_shrink24.gro -o em_shrink24_box.gro -box $xdim $ydim  14.5 -center $halfxdim $halfydim 5

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
echo "1 | r POPC \n q" | gmx make_ndx -f em.gro 

# Perform an equilibration in the NVT ensemble for 5 ns
gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -n index.ndx -o nvt
gmx mdrun -deffnm nvt -v

# Perform an equilibration in the NPT ensemble for 10 ns
gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -n index.ndx -o npt
gmx mdrun -deffnm npt -v
