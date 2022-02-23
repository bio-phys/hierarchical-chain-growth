#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from chain_growth.fragment_list import get_sequence_from_fasta, get_sequence_from_pdb
import os , pytest 

def check_sequence():#path_hcg_fulllength, path_ref):
    """
    check_sequence
    --------------
    script to check wether full-length IDP grown via HCG has the correct sequence
    == same as the input sequence
    """
    
    # path to topology file of idp grown with (r)hcg
    file_hcg = '../../examples/run_chain_growth/truncated_tauK18/4/0/pair0.pdb'
    sequence_hcg = get_sequence_from_pdb(file_hcg)
    
    # path to reference sequence == inout sequence
    file_ref =  '../../examples/truncated_tauK18.fasta'
    ## get amino acid sequence from PDB or fasta file in three letter code
    file_type = os.path.basename(file_ref).split('.')[-1]
    if file_type.upper() == 'FASTA':
        sequence_ref = get_sequence_from_fasta(file_ref)
    elif file_type.upper() == 'PDB':
        sequence_ref = get_sequence_from_pdb(file_ref)
        
    #print(sequence_hcg, sequence_ref)
    assert (sequence_hcg == sequence_ref)