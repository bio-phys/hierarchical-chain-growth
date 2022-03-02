#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
genFragment
-----------
prepare short overlapping fragments to generate a MD fragment library for (R)HCG
-> script divide fulllength sequence of an IDP into fragments of desired length
-> and overlap between subsequent fragments
=> generate tleap.txt file with fragments covering the sequence of the desired IDP
=> tleap.txt is the input file for Amber Tool "tleap", which can build molecules with a given sequence
"""
from chain_growth.fragment_list import generate_fragment_list

#############
## define the architecture of your fragments
#############

### get the sequence of your IDP ###
## either use a model of your IDP saved in a pdb 
## or save the sequence (eg from fasta file) 
sFile = '../truncated_tauK18.fasta'

## length of the fragments you want to build
fragmentLength = 5

## residue overlap in subsequent fragments
## default is 2
## pro overlap 1: less fragments needed to span your full-length protein sequence
## pro overlap 2: only central residues end up in full-length protein
## in both cases the dihedral angles of the residues
## that end up in the full-length protein are conserved comparably good
overlap = 2



## path to the folder you want to store the output tleap.txt
path2output = './' 


############
## generate fragment sequence with given length and overlap
############

fragmentL, overlap_d = generate_fragment_list(sFile, fragmentLength, overlap)
print(fragmentL, len(fragmentL))


#########
## generate "tleap" input file with fragments to build
#########

# open tleap.txt to store sequence of fragments to generate with tleap
tleap_seq = open('{}/tleap_short_tauK18.txt'.format(path2output) , 'w')

# tleap needs to load a forcefild for generating the fragments
tleap_seq.write('source leaprc.protein.ff14SB\n')

for  i, pep in enumerate(fragmentL):
    p = ''
    print(pep)
    for aa in pep:
        print(aa)
        p = p+' '+aa
    #print( p)
    
    tleap_seq.write('fragment = sequence {ACE %s NME}\nsavepdb fragment %i/fragment.pdb\n\n'%(p,i))
tleap_seq.write('quit')

tleap_seq.close()



