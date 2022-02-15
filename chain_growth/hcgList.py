#!/usr/bin/env python3

"""
hcgList
-------
prepare list of pairs needed as input for hcg
pairs = fragments or pairs of fragments that are assembled into pairs of fragments or pairs of pairs
"""

import numpy as np
import collections

def pair_fragments(input_ar, verbose=False):
    """ create pairs for each level
    original function from Lukas S. Stelzl, slightly modified by Lisa M. Pietrek
    
    Parameter
    ---------
    input_ar : list
        initial list of MD fragments/pairs of paired fragments
    verbose : Boolean
        set verbosity level
        
    Returns
    -------
    new_pairs : list
        new list with paired fragments/pairs of paired fragments
    promotion : boolean
        evaluation if fragment is promoted to the next hierarchy level
        False : two fragments can be paired 
        True : left-over fragment, that cannot be paired, is promoted 
    """
    promotion = False
    ## create list of paired fragments/pairs of paired fragments
    ## as pairs or list of pair for m > 1
    new_pairs = zip(input_ar[::2], input_ar[1::2])
    new_pairs = [list(p) for p in new_pairs]
    
    ## divide len(input_ar) by 2. If there is remainder (len(input_ar) == uneven) -> promote
    if len(input_ar) % 2:
        if verbose:
            print('promoting')
        ## add fragment/pairs of paired fragments to new_pairs without a partner to pair
        ## to promote it to the next hierarchy level
        promotion = True
        new_pairs.append(input_ar[-1])

    return new_pairs , promotion

def next_power_of_2(N): 
    """ get number of hierarchy levels for chain growth
        -> the next power of 2 from the number of fragments
    function created by Dr. Lukas S. Stelzl

    Parameters
    ----------
    N : integer
        total number of fragments
    
    Returns
    -------
    next power of two for N : integer    
    """
    
    ## calculate number of hierarchy levels m for the number of fragments used
    return 1 if N == 0 else 2**(N - 1).bit_length()

def flatten(pair_list):
    """ flattens list with pairs of fragments / pairs of paired fragments
    function created by Dr. Lukas S. Stelzl
    
    Parameters
    ----------
    pair_list : list 
        list of lists with paired fragments / pairs of paired fragments
    
    Returns
    -------
    pair_list : list
        flattened lists of pairs
    """
    
    if isinstance(pair_list, collections.Iterable):
        return [a for i in pair_list for a in flatten(i)]
    else:
        return [pair_list]

def make_hcl_l(N, n_to_c_term = True):
    """ create input list with fragments/ pairs of fragments to assemble in HCG to get the full-length chain
    original function from Dr. Lukas S. Stelzl, slightly modified by Lisa M. Pietrek

    Parameter
    ---------
    N : integer
        number of fagments
    n_to_c_term : boolean
        direction of growth 
        True : from N-terminus to C-terminus
        False : from C-terminus to N-terminus
        
    Returns:
    hcg_a[:,0] : list 
        list for HCG with lists of fragments assigned to be paired
    hcg_a[:,1] : list
        list with boolean, if here is an promotion in level m (last assembly step)
    """
    next_power2 = next_power_of_2(N)
    # from next power of 2 create max number of assemling levels -> M
    m = np.log2(next_power2).astype(int)
    frag_pair_l = list(np.arange(N))
    hcg_l = []
    
    ## loop through number of levels to generate list of paired fragment/ pairs of paired fragments
    ## and get promotion-evaluation
    ## save all pairs + promotion-evaluation per level in hcg_a
    for level_index in range(m):
        ## reverse initial fragment list to grow from C- to N-terminus
        if level_index == 0 and n_to_c_term == False:
            frag_pair_l.reverse()
        frag_pair_l, promo = pair_fragments(frag_pair_l)
        hcg_l.append([frag_pair_l, promo])
        
    hcg_a  = np.asarray(hcg_l, dtype=object)
   
    return hcg_a[:,0], hcg_a[:,1]


