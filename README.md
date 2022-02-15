# chain-growth
Grow ensembles of disordered biomolecules from fragment libraries

## Overview
Algorithm to assemble full-length chains of disordered proteins / regions from short overlapping fragments
We also provide example scripts
- to prepare the fragments to generate the MD fragment libraries using amber tleap
(- to run REMD simulations of the fragments on a SLURM-based cluster)
- to run hcg or rhcg + corresponding weights for a truncated tau K18

## References
Hierarchical Ensembles of Intrinsically Disordered Proteins at Atomic Resolution in Molecular Dynamics Simulations
Lisa M. Pietrek, Lukas S. Stelzl, and Gerhard Hummer
Journal of Chemical Theory and Computation 2020 16 (1), 725-737, https://pubs.acs.org/doi/abs/10.1021/acs.jctc.9b00809

Global Structure of the Intrinsically Disordered Protein Tau Emerges from its Local Structure
Lukas S. Stelzl, Lisa M. Pietrek, Andrea Holla, Javier S. Oroz, Mateusz Sikora, Jürgen Köfinger, Benjamin Schuler, Markus Zweckstetter, Gerhard Hummer
bioarxiv, http://dx.doi.org/10.1101/2021.11.23.469691
