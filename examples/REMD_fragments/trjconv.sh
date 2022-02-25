#!/bin/bash

sys=$1
traj=$sys"_prod.part0001.xtc"

tpr=$sys"_prod.tpr"
#echo " the trajectory to be processed: $traj"
#echo " the .tpr for the trajectory: $tpr"

#######
# uncomment in case of the mdrun did not finish in one run you need to concatenate the different part
#traj=$sys"_prod_concat.xtc"
#gmx trjcat -f *.xtc -o $traj
#


fname=$sys"_100ns"

PBCMOL=$fname"_mod.xtc"
echo " $PBCMOL "

### trjconv: processing
### strip out coordination, correct for periodicy, and alter traj
### EOF 1 1 > write only peptide wo    water into converted traj
gmx trjconv -s  $tpr -f $traj -o $PBCMOL -pbc mol -center -boxcenter tric -ur compact  << EOF
1
1
EOF

FITTED=$fname"_fitted.xtc"
gmx trjconv -s  $tpr -f $PBCMOL -o $FITTED -fit progressive << EOF
1
1
EOF

echo "processed: $traj using $tpr"
echo "produced : $FITTED"

#starting structure
gmx trjconv -s  $tpr -f $FITTED -dump 0 -o $sys"_100ns_nw.gro" << EOF
1
EOF

