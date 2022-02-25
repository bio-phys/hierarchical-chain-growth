Script-based workflow to setup and run replica exchange moleclular dynamics (REMD) simulations for a set of fragments constituting
a disordered protein.

Key script is **run.sh** which call Python and Bash scripts to setup directories for each fragment. For each fragments REMD 
simulations are setup. 

You can use (and adapt) the example .mdp files (= MD parameter) in the folder "run_para" for the MD simulations.


Paths have to be adapted in scripts. SLURM submission scripts have to be adapted to the cluster you use.

Workflow requires a working installation of GROMACS. 

For subsequent chain growth we can remove water and ions from the REMD trajectories (see trjconv.sh). 
