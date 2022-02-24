#!/bin/bash

# run gmx gromp for npt equilibration

sys=$1
mdp="/usr/path/run_para/npt.mdp"

mkdir -p equil
mkdir -p equil/npt
cd equil/npt

gmx grompp -f $mdp -c '../nvt/'$sys'_nvt.gro' -r '../nvt/'$sys'_nvt.tpr' -p '../../pdb2gmx/'$sys'.top' -o $sys'_npt.tpr'

cd ../../
