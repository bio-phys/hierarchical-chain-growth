#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
check_sequence
--------------
script to check wether full-length IDP grown via HCG has the correct sequence
== same as the input sequence
"""

from fragmentList import get_sequence_from_fasta, get_sequence_from_pdb
import os #, pytest 

# path to topology file of idp grown with (r)hcg
file_hcg = 'example/test/4/0/pair0.pdb'
sequence_hcg = get_sequence_from_pdb(file_hcg)

# path to reference sequence == inout sequence
file_ref =  'example/short_tau.fasta'
## get amino acid sequence from PDB or fasta file in three letter code
file_type = os.path.basename(file_ref).split('.')[-1]
if file_type.upper() == 'FASTA':
    sequence_ref = get_sequence_from_fasta(file_ref)
elif file_type.upper() == 'PDB':
    sequence_ref = get_sequence_from_pdb(file_ref)
    
#print(sequence_hcg, sequence_ref)
assert (sequence_hcg == sequence_ref)
