#!/bin/bash

#########
## Script for MD system setup
## Original version provided by Lukas S. Stelzl
## modified by Lisa M. Pietrek
#########

sys_name=$1
top=$sys_name".top"
pdb2gmx_name="$sys_name""_pdb2gmx.gro"
#ff_dir > force field directory
ff_dir="/usr/path/amber99sb-star-ildn-q.ff"
#md param
mdp="/usr/path/run_para/minim.mdp"
num_ion_script="/usr/path/script/num_ions.py"

#mkdir for outputfiles  of the pdb2gmx step
mkdir -p pdb2gmx
cd pdb2gmx
ln -s $ff_dir . 

gmx pdb2gmx -f ../*.pdb  -o $pdb2gmx_name -p $top -ignh << EOF
1
1
EOF

cd ..

##define box and solvate system
###new directory for sys solv output files
mkdir -p box
cd box

box_name="$sys_name""_box" 
#editconf: define box
gmx editconf -f ../pdb2gmx/$pdb2gmx_name -o $box_name -c -box 3.0 -bt dodecahedron

sol_name="$sys_name""_sol"
#genbox: solvate box
gmx solvate -cs spc216.gro -cp $box_name -o $sol_name -p ../pdb2gmx/$top

cd ..

#dir for genion> system neutralyzation>> addition of ions
mkdir -p genion
cd genion

prelim_tpr="$sys_name""_prelim.tpr"
#grompp > generation of tpr file for system neutralyzation > assembles coordinates, topology and md params (.mdp file)
gmx grompp -f $mdp -po "$sys_name"".mdp" -c ../box/$sol_name  -o $prelim_tpr -p ../pdb2gmx/$top  -maxwarn 1 


num_water=`grep SOL ../pdb2gmx/"$top"`

charge=`grep "qtot" ../pdb2gmx/"$top" | tail -n 1 | awk '{print $11}'`
#declare -i $charge

echo charge of the peptide is $charge

# http://unix.stackexchange.com/questions/53310/splitting-string-by-the-first-occurrence-of-a-delimiter
num_water=${num_water#*" "}

n_ion=`python $num_ion_script $num_water 0.15`
n_ion=${n_ion%%","*}

#extra_Na=0
NA=$n_ion
CL=$n_ion
solv_sys="$sys_name""_sys.pdb"
echo $CL
echo $NA
#genion cmd: adding ions > replace water with ions to neutralyze system 
gmx genion -p ../pdb2gmx/$top -pname NA -np $NA -nname CL -nn $CL -s $prelim_tpr  -o $solv_sys -neutral << EOF
13
EOF

cd ..

#dir for energy minimalization
mkdir -p em
mkdir -p em/steep1
cd em/steep1

em_name="$sys_name""_steep1"
ln -s $ff_dir . 

#grompp: composition of .tpr file for em 
##input for energy minimization
gmx grompp -f $mdp -po "$em_name"".mdp" -c ../../genion/$solv_sys  -o $em_name".tpr" -p ../../pdb2gmx/$top

#run em
gmx mdrun -v -deffnm $em_name
cd ../../

