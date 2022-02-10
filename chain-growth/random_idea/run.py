#!/usr/bin/env python3

from hcgList import make_hcl_l
from fragmentList import generate_fragment_list
#import sys

################
## prepare HCG
################

## input file and path / output path
## path to MD fragments
path0 = 'example/'
## path to store assembled models in
path = 'example/test/'
## file with sequence, format: "fasta" or "PDB"
sequence_f = 'example/short_tau.fasta'

## fragment construction
# length of MD fragments (without the end-capping groups if present)
fragment_length = 5
# length of the residue overlap between subsequent fragments
# == number of the same residues in subsequent fragments
overlap = 2

# generate list of fragments, dictionairy of overlaps between fragments
# generating overlaps_d is necessary, since to match the full-ltngh sequennce
# the overlap between e.g., the two last fragments can vary
fragment_l, overlaps_d = generate_fragment_list(sequence_f, fragment_length, overlap)

## lists for the HCG
n_pairs = len(fragment_l)
## hcg_l : list of paired fragments
## promo_l : list to evaluate if last fragment of level m in hcg_l is promoted to level m+1
hcg_l, promo_l = make_hcl_l(n_pairs)

# maximal number of pairs/full-length models to assemble
kmax = 100
# MD fragments are sampled with or without end-capping groups
capping_groups = True

if weighted:
    from hcgGrow_fct import hierarchical_chain_growth_w as hcg
    ###########
    ## run HCG
    ###########
    hierarchical_chain_growth(hcg_l, promo_l, overlaps_d, path0, path, kmax=kmax,
        capping_groups=capping_groups) #, verbose=True) 
    
else:
    from hcgGrow_fct import hierarchical_chain_growth 

###########
## run HCG
###########
hierarchical_chain_growth(hcg_l, promo_l, overlaps_d, path0, path, kmax=kmax,
        capping_groups=capping_groups) #, verbose=True) 
    

