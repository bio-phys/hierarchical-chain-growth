#!/usr/bin/env python3

import MDAnalysis as mda
import os

def get_sequence_from_pdb(file):
    """ generates list with amino acid sequence in three letter code

    Parameters
    ----------
    file : PDB file

    Returns
    -------
    three_list : list
        list with amino acid sequence in three letter code
    """
    
    ## make universe with PDB file
    u = mda.Universe(file)
    ## get amino acid sequence in three letter code - as an array
    three_list = u.residues.resnames
    three_list = three_list.tolist()
    
    return three_list

def get_sequence_from_fasta(file):
    """ generates list with amino acid sequence in three letter code
    
    Parameters
    ----------
    file : fasta file with protein sequence

    Returns
    -------
    three_list : list
        list with amino acid sequence in three letter code
    """

    ## dictionairy to translate one letter to three letter amino acid code
    three_letter ={'V':'VAL', 'I':'ILE', 'L':'LEU', 'E':'GLU', 'Q':'GLN', 'D':'ASP', 
            'N':'ASN', 'H':'HIS', 'W':'TRP', 'F':'PHE', 'Y':'TYR', 'R':'ARG', 'K':'LYS', 
            'S':'SER', 'T':'THR', 'M':'MET', 'A':'ALA',    'G':'GLY', 'P':'PRO', 'C':'CYS'}
    three_list = []
    ## open fasta file and read lines
    with open(file) as data:
        read=data.readlines()    
        for i, l in enumerate(read[1:]):  # exclude header in first line
            for sl in l:
                if sl!='\n': # exclude new line characters
                    ## store each amino acid in three letter code in three_list
                    three_list.append(three_letter[sl.upper()])
    
    return three_list

def generate_fragment_list(file, fragment_length, overlap, n_to_c_term = True):
    """ generates list of fragments of desired length, with desired overlap between subsequent fragments

    Parameters
    ----------
    file : file PDB or fasta format
        to get sequence from
    fragment_length : integer
        desired fragment length 
    overlap : integer
        desired fragment overlap 
        
    Returns
    -------
    fragment_list : list
            list with fragments of desired length and overlap
    overlaps : dictionary
        dict of overlaps between subsequent fragments
        key: number fragment, value: overlap
    """
    
    ## get amino acid sequence from PDB or fasta file in three letter code
    file_type = os.path.basename(file).split('.')[-1]
    if file_type.upper() == 'FASTA':
        three_list = get_sequence_from_fasta(file)
    elif file_type.upper() == 'PDB':
        three_list = get_sequence_from_pdb(file)
    else:
        print('sequence file: unknown file type. Allowed: ".fasta" or ".pdb"')
      
    ## generate list of fragments with desired length
    fragment_list = []
    overlaps = {}
    ## counnter fragment 
    count = 0
    ## loop through range of the sequence length with as step size==fragment_length-overlap -> generate desired overlap
    steps = fragment_length-overlap
    for i in range(0, len(three_list) , steps): 
        ## if generated fragment has desired length as defined with fragment_length
        ## append fragment to fragment_list and "typical" overlap to overlaps
        if len(three_list[i:i+fragment_length]) == fragment_length:
            fragment_list.append(three_list[i:i+fragment_length])
            overlaps[count]=overlap
            count += 1
        ## if reached end of the sequence and last generated fragment has not desired length as defined with fragment_length
        ## go back in sequence such that the length of the new fragment == fragment_length and append to fragment_list
        ## re-evaluate new overlap with previous fragment and apped to overlaps
        else: 
            stepsBackInSeq = fragment_length - len(three_list[i:])
            newOverlap = overlap + stepsBackInSeq
            fragment_list.append(three_list[i-stepsBackInSeq:])
            overlaps[count]=newOverlap
            count += 1
    if overlaps[count-1] >= fragment_length:
        fragment_list.pop(-1)
        del overlaps[count-1]
    if n_to_c_term == False:
        fragment_list.reverse()
        
    return fragment_list, overlaps
