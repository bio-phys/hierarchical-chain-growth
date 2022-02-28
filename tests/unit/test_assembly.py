
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from chain_growth.hcg_list import make_hcl_l, flatten
from chain_growth.fragment_list import get_sequence_from_fasta, get_sequence_from_pdb, generate_fragment_list
from chain_growth.hcg_fct import hierarchical_chain_growth
import os , pytest, pickle 
import numpy as np
import MDAnalysis as mda

test_dir = os.path.dirname(os.path.abspath(__file__))

def test_make_hcl_l_with_promotion():
    '''Check that the fragment assembly schedule is generated as expected. This test
       is based on Fig. S1 
       https://pubs.acs.org/doi/suppl/10.1021/acs.jctc.9b00809/suppl_file/ct9b00809_si_001.pdf
    '''
    n_fragments = 46 # we could
    hcg_l, promo_l = make_hcl_l(n_fragments)
    
    ##  hcg_l, promo_l index+1 equals level of assembly m.
    ##  index 0, level 1: pairs of fragments (level 0 are fragments)
    
    ## one promotion going from level 1->2
    assert (promo_l[0] == True)
    
    ##  one promotion going from level 4-> 5
    assert (promo_l[4] == True)
    
    ## We also expect that the two unequal fragments for the penulitmative level
    assert( flatten(hcg_l[-2][1]).__len__() == 14)
    ## Out[29]: 14
    assert( flatten(hcg_l[-2][0]).__len__() == 32)
    ## Out[30]: 32

    ## We and finally one 46 mer of fragments at the last level
    assert(flatten(hcg_l[-1][0]).__len__() == 46)
    ## Out[31]: 46
    


