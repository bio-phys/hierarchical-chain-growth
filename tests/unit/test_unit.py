
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from chain_growth.hcg_list import make_hcl_l
from chain_growth.fragment_list import get_sequence_from_fasta, get_sequence_from_pdb, generate_fragment_list
from chain_growth.hcg_fct import hierarchical_chain_growth
import os , pytest, pickle 
import numpy as np
import MDAnalysis as mda

test_dir = os.path.dirname(os.path.abspath(__file__))

def test_make_hcl_l():
    '''Check that the fragment assembly schedule is generated as expected
    '''
