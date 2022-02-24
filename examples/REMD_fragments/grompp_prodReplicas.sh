#!/bin/bash

# run gmx gromp for REMD production run

sys=$1

mkdir -p prod
cd prod

for dir in ./*/; do
cd $dir
gmx grompp -f prod*.mdp -c '../../equil/npt/'$sys'_npt.gro' -t '../../equil/npt/'$sys'_npt.cpt'  -p '../../pdb2gmx/'$sys'.top' -o $sys'_prod.tpr'
cd ../ ;
done

cd ../
