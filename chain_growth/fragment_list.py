#!/usr/bin/env python3
"""
fragmentList
-----------
prepare a list of fragments with desired length + overlap from input sequence
== IDP/ IDR to grow using (r)hcg
"""
import MDAnalysis as mda
import os

class get_sequence:
    def __init__(self, input_sequence):
        """
        Parameters
        ----------
        input_sequence : string
            Either path to sequence file (PDB or fasta) or seqeunce as string in one letter
            amino acid code.

        """
        ## dictionairy to translate one letter to three letter amino acid code
        self.three_letter ={'V':'VAL', 'I':'ILE', 'L':'LEU', 'E':'GLU', 'Q':'GLN', 'D':'ASP', 
                'N':'ASN', 'H':'HIS', 'W':'TRP', 'F':'PHE', 'Y':'TYR', 'R':'ARG', 'K':'LYS', 
                'S':'SER', 'T':'THR', 'M':'MET', 'A':'ALA',    'G':'GLY', 'P':'PRO', 'C':'CYS'}
        self.input_sequence = input_sequence
        

    def _from_pdb(self):
        """
        Returns
        -------
        sequence_three_letter_list : list
            created list with amino acid sequence in three letter code from PDB file
        """

        ## make universe with PDB file
        u = mda.Universe(self.input_sequence)
        ## get amino acid sequence in three letter code - as an array
        sequence_three_letter_list = u.residues.resnames
        sequence_three_letter_list = sequence_three_letter_list.tolist()

        return sequence_three_letter_list

    def _from_fasta(self):
        """
        Returns
        -------
        sequence_three_letter_list : list
            created list with amino acid sequence in three letter code from fasta file
        """

        # sequence_three_letter_list = []
        ## open fasta file and read lines
        with open(self.input_sequence) as data:
            read = data.readlines()    
            sequence_three_letter_list = [self.three_letter[sl.upper()] for i, l in enumerate(read[1:])
                                             for sl in l if sl!='\n']        
        return sequence_three_letter_list

    def _from_string(self):
        """         
        Returns
        -------
        sequence_three_letter_list : list
            created list with amino acid sequence in three letter code from string
        """
        return [self.three_letter[si.upper()] for si in self.input_sequence]

    def three_letter_list(self):
        input_basename = os.path.basename(self.input_sequence).split('.')
        try:
            if len(input_basename) == 2:
                if input_basename[-1].upper() == 'FASTA':
                    sequence_three_letter_list = self._from_fasta()
                elif input_basename[-1].upper() == 'PDB':
                    sequence_three_letter_list = self._from_pdb()

            elif len(input_basename) == 1:
                sequence_three_letter_list = self._from_string()
            return sequence_three_letter_list
        except (NameError , TypeError):
            print('Invalid sequence input. Allowed: ".fasta" ,".pdb" or string')


def generate_fragment_list(input_sequence, fragment_length, overlap, n_to_c_term = True):
    """ generates list of fragments of desired length, with desired overlap between subsequent fragments

    Parameters
    ----------
    input_sequence : string
        path to PDB or fasta file to get sequence from
        or amino acid one letter code as string
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
    sequence_three_letter_list = get_sequence(input_sequence).three_letter_list()
    
    ## generate list of fragments with desired length
    fragment_list = []
    overlaps = {}
    ## counnter fragment 
    count = 0
    ## loop through range of the sequence length with as step size==fragment_length-overlap -> generate desired overlap
    steps = fragment_length-overlap
    for i in range(0, len(sequence_three_letter_list) , steps): 
        ## if generated fragment has desired length as defined with fragment_length
        ## append fragment to fragment_list and "typical" overlap to overlaps
        if len(sequence_three_letter_list[i:i+fragment_length]) == fragment_length:
            fragment_list.append(sequence_three_letter_list[i:i+fragment_length])
            overlaps[count]=overlap
            count += 1
        ## if reached end of the sequence and last generated fragment has not desired length as defined with fragment_length
        ## go back in sequence such that the length of the new fragment == fragment_length and append to fragment_list
        ## re-evaluate new overlap with previous fragment and apped to overlaps
        else: 
            stepsBackInSeq = fragment_length - len(sequence_three_letter_list[i:])
            newOverlap = overlap + stepsBackInSeq
            fragment_list.append(sequence_three_letter_list[i-stepsBackInSeq:])
            overlaps[count]=newOverlap
            count += 1
    if overlaps[count-1] >= fragment_length:
        fragment_list.pop(-1)
        del overlaps[count-1]
    if n_to_c_term == False:
        fragment_list.reverse()
        
    return fragment_list, overlaps
