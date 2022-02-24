#!/bin/bash

######
## This bash script gives and example workflow we suggest to follow to create the MD fragment library
## by running REMD simulations of the individual fragments.
## You can either chose to run all the subsequent steps manually in the terminal
## or execute this file ones to submit everthing to the cluster with job dependencies 
## (once everything is tested and configured to match your setup, to do so uncomment respective line).
## make sure to check the path to .mdp files or other scripts needed before you run the scripts 
######


dire=$1
sys=$2
s=$3
e=$4

## scripts
path_to_scripts=$5


## provide jobid, if you want to sbmit to the cluster with a dependency for first job
JOBID=$6

cd $dire

#######
## create fragment
#######
for i in $(seq $s $e) ; do
 mkdir -p $i ;
done
tleap -f tleap.txt
# need to quit tleap manually


#######
## setup system and energy min
#######
## on terminal 
for i in $(seq $s $e) ; do
 cd $i
 bash $path_to_scripts"setup_sys.sh" $sys 
 cd .. ;
done
#
## or submit job on cluster to run setup sys and emin for each fragment in a loop
#JOBID=$(sbatch -J setup --dependency=afterok:$JOBID $path_to_scripts"setup_sys.job" $path_to_scripts"setup_sys.sh" $s $e $sys  2>&1  | awk '{print $(NF)}')


#######
## nvt
#######
## run grompp on terminal 
for i in $(seq $s $e) ; do
 cd $i
 bash  $path_to_scripts"gromppNVT.sh" $sys
 cd .. ;
done
#
## or submit job on cluster to run gromppNVT for each fragment in a loop
#JOBID=$(sbatch -J gromppNVT --dependency=afterok:$JOBID $path_to_scripts"single_core.job" $path_to_scripts"gromppNVT.sh" $s $e $sys 2>&1  | awk '{print $(NF)}')
#
#
## run equil
for i in $(seq $s $e) ; do
 cd $i/equil/nvt
 if [[ $i < $e ]] ;
  then sbatch -J nvt$i  --dependency=afterok:$JOBID  $path_to_scripts"nvt.job" $sys 
 else ;
  JOBID=$(sbatch -J nvt$i  --dependency=afterok:$JOBID  $path_to_scripts"nvt.job" $sys 2>&1  | awk '{print $(NF)}')
 fi
 cd ../../../ ;
done


#######
## npt
#######
## run grompp on terminal 
for i in $(seq $s $e) ; do
 cd $i
 bash  $path_to_scripts"gromppNPT.sh" $sys
 cd .. ;
done
#
## or   submit job on cluster to run gromppNPT for each fragment in a loop
#JOBID=$(sbatch -J gromppNPT --dependency=afterok:$JOBID $path_to_scripts"single_core.job" $path_to_scripts"gromppNPT.sh" $s $e $sys 2>&1  | awk '{print $(NF)}')
#
#
## run equil
for i in $(seq $s $e) ; do
 cd $i/equil/npt
 if [[ $i < $e ]]  ;
  then  sbatch -J npt$i  --dependency=afterok:$JOBID  $path_to_scripts"npt.job" $sys
 else 
  JOBID=$(sbatch -J npt$i  --dependency=afterok:$JOBID  $path_to_scripts"npt.job"  $sys 2>&1  | awk '{print $(NF)}')
 fi
 cd ../../../ ;
done


#######
## prod
#######
# run make replikas on terminal
for i in $(seq $s $e) ; do
 cd $i
 mkdir -p prod
 cd prod
 python $path_to_scripts"make_replicas.py"
 cd ../../
done
#
## or submit job on cluster to run make replica for each fragment in a loop
#JOBID=$(sbatch -J prepRepl --dependency=afterok:$JOBID $path_to_scripts"single_core.job" $path_to_scripts"call_make_replikas.sh"  $s $e $sys 2>&1  | awk '{print $(NF)}')
#
#
## run grompp on terminal 
for i in $(seq $s $e) ; do
 cd $i
 bash  $path_to_scripts"grompp_prodReplicas.sh" $sys
 cd .. ;
done
#
## or submit job on cluster to run setup sys and emin for each fragment in a loop
#JOBID=$(sbatch -J gromppP --dependency=afterok:$JOBID $path_to_scripts"single_core.job" $path_to_scripts"grompp_prodReplicas.sh" $s $e $sys 2>&1  | awk '{print $(NF)}')
#
#
## run REMD
for i in $(seq $s $e) ; do
 cd $i/prod
 sbatch -J remd$i  --dependency=afterok:$JOBID $path_to_scripts"re_prod.job" $sys
 cd ../../../ ;
done

