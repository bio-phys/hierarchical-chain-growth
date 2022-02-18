#!/usr/bin/env python3
"""
run_rhcg
--------
run reweighted hierarchical chain growth
"""
from chain_growth.hcg_list import make_hcl_l
from chain_growth.fragment_list import generate_fragment_list
from chain_growth.hcg_fct import reweighted_hierarchical_chain_growth

################
## prepare HCG
################

## input file and path / output path
# path to MD fragments
path0 = '..'
# path to store assembled models in
path = 'truncated_tauK18_weighted/'
#path weights
path2weights = '../weights'
# file with sequence, format: "fasta" or "PDB"
sequence_f = '../truncated_tauK18.fasta'
# theta for fragment reweighting
theta=10.0

## define fragment construction
# length of MD fragments (without the end-capping groups if present)
fragment_length = 5
# length of the residue overlap between subsequent fragments
# == number of the same residues in subsequent fragments
overlap = 2

# generate list of fragments, dictionary of overlaps between fragments
# generating overlaps_d is necessary, since to match the full-ltngh sequennce
# the overlap between e.g., the two last fragments can vary
fragment_l, overlaps_d = generate_fragment_list(sequence_f, fragment_length, overlap)

## lists for the HCG
n_pairs = len(fragment_l)
print(fragment_l, overlaps_d, n_pairs)
# hcg_l : list of paired fragments
# promo_l : list to evaluate if last fragment of level m in hcg_l is promoted to level m+1
hcg_l, promo_l = make_hcl_l(n_pairs)
print(hcg_l)

# maximal number of pairs/full-length models to assemble
kmax = 10
# MD fragments are sampled with or without end-capping groups
capping_groups = True

###########
## run HCG
###########
reweighted_hierarchical_chain_growth(hcg_l, promo_l, overlaps_d, path0, path, kmax=kmax,
        capping_groups=capping_groups,  path2weights=path2weights, theta=theta)#, verbose=True) 
    

