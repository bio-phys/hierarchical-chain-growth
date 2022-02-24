#!/bin/bash

# run gmx gromp for nvt equilibration

sys=$1

mdp='/usr/path/run_para/nvt.mdp '

mkdir equil
mkdir equil/nvt
cd equil/nvt

gmx grompp -f $mdp -c "../../em/steep1/"$sys"_steep1.gro" -r "../../em/steep1/"$sys"_steep1.tpr" -p "../../pdb2gmx/"$sys".top" -o $sys"_nvt.tpr"

cd ../../
