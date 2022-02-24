#!/bin/bash -l
# Standard output and error:
#SBATCH -o ./tjob.out.%j
#SBATCH -e ./tjob.err.%j
# Initial working directory:
#SBATCH -D ./
# Job Name:
#SBATCH -J name
#
# Queue (Partition):
#SBATCH --partition=s.phys  
#SBATCH --gres=gpu:rtx6000:3 
#SBATCH --mem=170000  
#
# Number of nodes and MPI tasks per node:
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --cpus-per-task=2
# SBATCH --ntasks-per-core=2 # hyperthreading
#SBATCH --mail-type=none
#SBATCH --mail-user=<userid>@rzg.mpg.de
#
# Wall clock limit:
#SBATCH --time=24:00:00

#module load purge
module load intel/19.1.3
module load impi/2019.9
module load cuda/11.4
module load gromacs/2019.6


export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MPI_NUM_RANKS=$SLURM_NTASKS_PER_NODE
export OMP_PLACES=cores  ## with enabled hyperthreading this line needs to be commented out

#export OMP_PLACES=threads ## with enabled hyperthreading these two lines need to be uncommented
#export SLURM_HINT=multithread


sys=$1


srun gmx_mpi mdrun -v -pin on  -noappend -ntomp $OMP_NUM_THREADS -maxh 23.9 -deffnm $sys"_prod" -cpi $sys"_prod_prev.cpt" -replex 500 -multidir 278.00 283.54 289.17 294.90 300.72 306.63 312.64 318.75 324.95 331.26 337.67 344.18 350.79 357.52 364.35 371.29 378.35 385.51 392.80 400.22 407.75 415.40 423.17 431.08
