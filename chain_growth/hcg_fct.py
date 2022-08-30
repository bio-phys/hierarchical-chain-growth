#!/usr/bin/env python3
"""
hcg_fct_multiproc
-----------
core functions for hierarchical chain growth

run hcg in parallel - per level -> m_i loop is paralleliiiized
using multipprocessing.Pool
"""
import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import align
import pathlib2, shutil, os
import MDAnalysis.analysis.distances as distances
from chain_growth.hcg_list import flatten
from multiprocessing import Pool
from functools import partial

def translate_concept(u1, u2, proline_2nd_posi, align_begin1, align_end1, 
                      align_begin2, align_end2):
    """ prepare a dictionary with atoms to align before the assembly step
    
    Parameter
    ---------
    u1 : universe
    u2 : universe
    proline_2nd_posi : boolean
        Whether or not in  u2 2nd residue to align is a proline 
    align_begin1 : integer
        residue index for u1, marks residue to begin the alignment for u1
    align_end1 : integer
        residue index for u1, marks residue to end the alignment for u1
    align_begin2 : integer
        residue index for u2, marks residue to begin the alignment for u2
    align_end2 : integer
        residue index for u2, marks residue to end the alignment for u2
        
    Returns
    -------
    select : dictionary 
        dictionary with atoms that are aligned prior to the assembly step
    """
    ## whether or not a amide proton is there
    if proline_2nd_posi:
        atoms_sel_2nd_res = "name  N"
    else:
        atoms_sel_2nd_res= "name  N or name H"
    
    select={}
    ## atom selection of fragment1 to be aligned
    sel1 = "(resid {} and (name C or name O)) or (resid {} and ({}))".format(
                            u1.atoms.residues[align_begin1].resid, 
                            u1.atoms.residues[align_end1].resid, atoms_sel_2nd_res)
    ## atom selection of fragment2 to be aligned
    sel2 = "(resid {}  and (name C or name O)) or (resid {} and ({}))".format(
                           u2.atoms.residues[align_begin2].resid,
                           u2.atoms.residues[align_end2].resid, atoms_sel_2nd_res)

    ## dict is argument for mda alignto function: arg = "select="
    ## mobile is superimposed on ref with mda.alignto in the assembly step
    select={'mobile': sel1, 'reference': sel2}
    return select


def find_clashes(u1,u2, index1b, index2e, 
                 index1e=-1, clash_radius=2.0):
    """ finds clashes: atoms with a distance below clash_radius
    
    Parameters
    ----------
    u1 : universe
        fragment 1 to pair
    u2 : universe
        fragment 2 to pair
    index1b : integer
        residue index for residues to *exclude* from clash detection
        marks residue to begin exclusion for u1
    index2e : integer
        residue index for residues to *exclude* from clash detection
        marks residue to end exclusion for u2
    index1e : integer
        residue index for residues to *exclude* from clash detection
        marks residue to end exclusion for u1. The default is -1
    clash_distance : float
        max. allowed distance between atoms, everything below is counted as clash
        The  default is 2.0

    Returns
    -------
    clsum : integer
        count of clashes between aligned fragments
        
    NOTE: MDAnalysis is inclusive! resid 1:2 -> selects residues 1+2
    """
    l1 = u1.select_atoms("protein and not (type H) and not (resid {} and backbone) and not resid {}:{}".format(
                                                    u1.atoms.residues[index1b].resid,
                                                    u1.atoms.residues[index1b+1].resid,
                                                    u1.atoms.residues[index1e].resid))

    # atom selection of u2 to scan for clashes
    l2 = u2.select_atoms("protein and not (type H) and not resid 1:{} and not (resid {} and backbone)".format(
                                                    u2.atoms.residues[index2e-1].resid,
                                                    u2.atoms.residues[index2e].resid))

    # -> use mda.distances to generate a matrix with distances
    distmat = distances.distance_array(l1.positions,l2.positions)
    cont = np.less(distmat,clash_radius) # see where distmat < cutoff
    cont = np.where(cont,1,0) # True --> 1, False --> 0
    clsum = np.sum(cont)
    
    return clsum


def merge_universe(u1,u2, merge_begin1, merge_end1, 
                   merge_begin2, merge_end2):
    """ assemble aligned universes and renumber residues in a consecutive manner
    
    Parameters
    ----------
    u1 : universe
        fragment 1 to pair
    u2 : universe
        fragment 2 to pair
    merge_begin1 : integer
        residue index to begin merge for u1
    merge_begin1 : integer
        residue index to end merge for u1
    merge_begin2 : integer
            residue index to begin merge for u2
    merge_end2 : integer
            residue index to end merge for u2
            
    Returns
    ------
    u : universe
        combined universes u1 and u2 with renumbered residues
    """
    
    # atom selection of u1 to be merged
    sel1="resid {}:{}".format(u1.atoms.residues[merge_begin1].resid , u1.atoms.residues[merge_end1].resid)
    sela = u1.select_atoms(sel1)
    # atom selection of u2 to be merged
    sel2="resid {}:{}".format(u2.atoms.residues[merge_begin2].resid , u2.atoms.residues[merge_end2].resid)
    selb = u2.select_atoms(sel2)
    u = mda.core.universe.Merge(sela, selb)
    u_atm = u.select_atoms('all')
    # renumber residue IDs after assembly
    u_atm.residues.resids = np.arange(1,len(u_atm.residues.resids)+1) 

    return u


def get_residue_indices_for_assembly(overlap0, current_overlap, capping_groups,
                                     last_level, verbose):
    """ assign residue indices for alignment, clash detection, assembly 
    -> indices depend on the overlap between subsequent fragments
    
    Parameters
    ----------
    current_overlap : integer
        overlap between the fragments to assemble at current step
    overlap0 : initial overlap, optional
        overlap initially chosen for the fragments without taking into account exceptions. The default is 2
    capping_groups : boolean, optional
        MD fragment are sampled with or without end-capping groups. The default is True
    last_level : boolean, optional
        last level of hierarchical chain growth. The default is False
    
    Returns
    -------
    index_aln_l : list
       residue indices for alignment, used in translate_concept
    index_clash_l : list
       residue indices for clash detection, used in find_clashes
    index_merge_l : list
        residue indices for fragment assembly, used in merge_universe
    """
    ## e: additional factor for end-capping groups
    ## having 1 overlapping residue + align peptide bonds does work only with headgroup!
    e = 1
    if overlap0  == 0:
        e = 0
        align_begin1= -2
        align_end1= -1
        align_begin2= 0     
    else:
        if overlap0 > 1 and capping_groups == False:
            e = 0
        # align peptide bond between two last/ two first residues 
        align_begin1= -(overlap0 + e)
        align_end1= align_begin1+1
        align_begin2= 0 + e
        
    if current_overlap != overlap0:
        align_begin2= np.abs(current_overlap - overlap0) + e
    align_end2= align_begin2+1        
    
        
    # indices to exclude residues from clash calculation
    index1_clashB = align_begin1
    index2_clashE = align_end2
        
    # indicies for assembly of aligned pairs    
    merge_begin1 = 0 # always first residue, maybe change this to make it more flexible??
    merge_end1 = align_begin1
    merge_begin2 = align_end2
    merge_end2 = -1 #always last residue, maybe change this to make it more flexible??
    
    # exclude capping groups
    if capping_groups and last_level:
        merge_begin1 = 1
        merge_end2 = -2
    
    index_aln_l = [align_begin1, align_end1, align_begin2, align_end2]
    index_clash_l = [index1_clashB, index2_clashE]
    index_merge_l = [merge_begin1, merge_end1, merge_begin2, merge_end2]
    
    if verbose:
        print('overlap between fragment 1 & 2:' , current_overlap,
               'align begin / end in fragment 1 & 2, respectively: ', index_aln_l,
               'clash search begin / end in fragment 1 & 2, respectively: ', index_clash_l,
               'merge begin / end in fragment 1 & 2, respectively: ', index_merge_l)
       
    return index_aln_l, index_clash_l, index_merge_l
        
def fragment_assembly(u1, u2, dire, select, index_clash_l, index_merge_l,
         rmsd_cut_off, clash_distance, kmax, ri_l=None, draw_indices=True):
    """ assemble the fragments to pairs
    
    Parameters
    ----------
    u1 : universe
        fragment 1 to pair
    u2 : universe
        fragment 1 to pair
    dire : path
        path to store assembled pair in
    select : dictionary 
        dictionary with resids and atoms
        of fragment 1 and 2 for alignment
    index_clash_l : list
        list with residue indices of fragment 1 and 2 for clash detection
        residues that should be excluded from neighbour search
    index_merge_l : list
        list with residue indices of fragment 1 and 2 for the final assemby
    rmsd_cut_off : float
        cut-off for the RMSD of the fragment alignment
    clash_distance : float
        max. allowed distance between atoms
    kmax : integer
        number of pairs that should be assembled in level m_i
    ri_l : array-like
        array with indices for chosing a specific confoormation of a fragment. The default is None
        if None: draw indices randomly
    draw_indices : booolean
        if new random integers == frame indices are drawn or else taken from a input array.

    Returns
    -------
    if draw_indices:
        rs = list of lists
        successful indices drawn for fragment 1 and 2 during fragment assembly
    else:
        None
    """
    
    k = 0
    assembly_atempt = 0
    
    if np.all(ri_l) is None:
        # array to store random frame indices of successfully assembled fragments / pairs       
        rs = np.zeros((kmax, 2))
    writePDB = True
    while k < kmax:
        if draw_indices:
            # random integer to draw random frame
            r1 = np.random.randint(u1.trajectory.n_frames)
            r2 = np.random.randint(u2.trajectory.n_frames)
        else:
            r1 = int(ri_l[k][0])
            r2 = int(ri_l[k][1])
        # load random frame
        u1.trajectory[r1]
        u2.trajectory[r2]
        
        # rmsd suportimposition of the subsequent fragmenmt         
        old,new = align.alignto(u1, u2, select=select, weights="mass",
                                tol_mass=5., match_atoms=False)
        assembly_atempt += 1
        if new < rmsd_cut_off:
            # calculate clashes
            clashes = find_clashes(u1 , u2 , index1b=index_clash_l[0] , index2e=index_clash_l[1],
                                  clash_radius=clash_distance)
            
            if clashes < 1:
                ## assemble the subsequent fragments 
                u = merge_universe(u1, u2, *index_merge_l)
                
                if writePDB :
                   # save atom positions + topology for first pair in a pdb file
                    u.atoms.write("{}/pair{}.pdb".format(dire,k))
                    at = np.zeros((kmax , u.atoms.n_atoms , 3))
                    writePDB = False

                # save atom positions of subsequehnt pairs in array
                at[k] = u.atoms.positions
                if draw_indices:
                    rs[k] = [r1,r2]
                k+=1
        
    # load and save universe with pair0 as topology and atom positions 
    # of subsequent paired fragment of this assembly step as trajectory
    # n frames = desired number of pairs in this assembly step
    allPairs = mda.Universe("{}/pair0.pdb".format(dire) , at)
    allPairs.atoms.write("{}/pair.xtc".format(dire) , frames='all')
    
    # return random integer array if none was given as input
    if draw_indices:
        return rs
    else:
        return None

def hierarchical_chain_growth(hcg_l, promo_l, overlaps_d, path0, path, kmax, 
             online_fragment_library=False, dict_to_fragment_folder=None,
             rmsd_cut_off=0.6, clash_distance=2.0, capping_groups=True,
             ri_l=None, verbose=False):
    """ perform hierarchical chain growth 
    assemble fragments/ pairs of fragments until reaching the full-length chain
    by calling _loop_func -> does the inner loop and calls fragment assembly
    takes same arguments as hcg function

    Parameters
    ----------
    hcg_l : list
        list for HCG with lists of fragments assigned to be paired
    promo_l : list
        list with boolean, if here is an promotion in level m (last assembly step)
    overlaps_d : dictionary
        dict of overlaps between subsequent fragments
    path0 : string
        path to the MD fragments 
    path : string
        path to folders where the assembled pairs are stored in

    path : path
        path to the MD fragments (and to folders where the assembled pairs are stored in)
    kmax : integer
        number of pairs that should be assembled in level m_i
    online_fragment_library : boolean
         wheter to draw from a fragment library you sampled yourself or 
         a pre-sampled fragment library that is available online for the chain assembly.
         The default is False
    dict_to_fragment_folder : dictionary
        dictionary that transaltes between the input sequence - required fragments - 
        and the code. The default is None
    rmsd_cut_off : float, optional
        cut-off for the RMSD of the fragment alignment. The default is 0.6
    clash_distance : float, optional
        max. allowed distance between atoms. The default is 2.0
    capping_groups : boolean, optional
        MD fragment are sampled with or without end-capping groups. The default is True
    ri_l : array-like
        array with indices for chosing a specific confoormation of a fragment. The default is None
    draw_indices : booolean
        if new random integers == frame indices are drawn or else taken from a input array.
   
    Returns
    -------
    None.
    """
    
    last_level = False
    cpu = os.cpu_count()
    for m , fragment_l in enumerate(hcg_l): 
        # folder to save assembled pair in this level (m+1)
        level = m+1
        # folder to get old pairs from previous level (m)
        previous_level = m
        promotion  = promo_l[m]
 
        try:
            ## if no index list is supplied draw new indices (fragment conformations)
            if np.all(ri_l) is None:
                draw_indices = True
                r_l = []
        except:
            r_l = ri_l[m]
            draw_indices = False
            
        # if MD fragments are sampled with end-capping_groups, 
        # they are removed in the last assembly step 
        if m+1 == len(hcg_l):
            last_level = True
    
        if len(fragment_l) > cpu:
            num_threads = cpu
        else:
            num_threads = len(fragment_l)
        pairs = [(m_i, pair_l) for m_i, pair_l in enumerate(fragment_l)]
        d = {"path0": path0, "path": path, 
             "online_fragment_library": online_fragment_library, "dict_to_fragment_folder": dict_to_fragment_folder,
             "level": level, 'previous_level': previous_level,
             "last_level": last_level, "fragment_l": fragment_l, "rmsd_cutoff": rmsd_cut_off,
             "clash_distance": clash_distance, "kmax": kmax, "r_l":r_l,
             "overlap_d": overlaps_d, "promotion": promotion, "capping_groups": capping_groups,
             "draw_indices" : draw_indices, "verbose": verbose }
        # POOL LOOP application
        with Pool(num_threads) as p: 
            func = partial(_loop_func, d)
            results = p.map(func, pairs)
        if draw_indices:
            for r in results:
                if r is not None:
                    r_l += [r]
                    
            np.save("{}/confIndex_level{}.npy".format(path, m+1), r_l) 
    return None


def _loop_func(variables, pairs):
    """ does the inner loop of hierarchical chain growth and calls fragment assembly
        
    idea for possible implementation - inspired from discussion with / and 
    realized with input from Hendrik Jung    

    Parameters
    ----------
    variables : dictionary
        variables needed for fragment assembly as defined in hcg function.
    pairs : tuple of lists
        pair_l of pairs need to be assembled per level
        list of respective pair ids
   
    Returns
    -------
    if draw_indices:
        rs = list of lists
        successful indices drawn for fragment 1 and 2 during fragment assembly
    else:
        None
    """
    
    # print(variables)
    m_i, pair_l = pairs
    path0 = variables["path0"]
    path = variables["path"]
    online_fragment_library = variables["online_fragment_library"]
    dict_to_fragment_folder = variables["dict_to_fragment_folder"]
    level = variables["level"]
    previous_level = variables["previous_level"]
    last_level = variables["last_level"]
    fragment_l = variables["fragment_l"]
    rmsd_cut_off = variables["rmsd_cutoff"]
    clash_distance = variables["clash_distance"]
    k_max =  variables["kmax"]
    r_l = variables["r_l"]
    draw_indices = variables["draw_indices"]
    overlaps_d = variables["overlap_d"]
    promotion = variables["promotion"]
    capping_groups = variables["capping_groups"]
    verbose = variables["verbose"]
    
    overlap = overlaps_d[0]

    # if MD fragments are sampled with end-capping_groups, 
    # they are removed in the last assembly step 

    proline_2nd_posi = False
    # if promotion of MD fragment in first level DO NOT define old_pair2
    # index ERROR because in this case, pair_l is no list
    if promotion and isinstance(pair_l, list) == False:
        # subfolder fragment 1, 2
        old_pair1  = pair_l
        old_pair2  = pair_l
    else:
        # subfolder fragment 1, 2
        old_pair1 = flatten(pair_l[0])[0]
        old_pair2 = flatten(pair_l[1])[0]
        
    if previous_level == 0 and online_fragment_library:
        ## path0 = path to online fragment library
        path2fragment = path0 
        #print(old_pair1, old_pair2)
        old_pair1_fragmentLib = dict_to_fragment_folder[old_pair1]
        old_pair2_fragmentLib = dict_to_fragment_folder[old_pair2]
        #print(old_pair1_fragmentLib, old_pair2_fragmentLib)
        old_dire1 = '{}/{}'.format(path2fragment, old_pair1_fragmentLib)
        old_dire2 = '{}/{}'.format(path2fragment, old_pair2_fragmentLib)
        #print(old_dire1, old_dire2)
        top = 'fragment.pdb'
        xtc = 'fragment.xtc'
                
    elif previous_level == 0 and online_fragment_library == False:
        previous_level = 'MDfragments'             
        path2fragment = path0
        old_dire1 = '{}/{}/{}'.format(path2fragment, previous_level, old_pair1)
        old_dire2 = '{}/{}/{}'.format(path2fragment, previous_level, old_pair2)
        top = 'pair0.pdb'
        xtc = 'pair.xtc'
    else:
        path2fragment = path 
        old_dire1 = '{}/{}/{}'.format(path2fragment, previous_level, old_pair1)
        old_dire2 = '{}/{}/{}'.format(path2fragment, previous_level, old_pair2)
        top = 'pair0.pdb'
        xtc = 'pair.xtc'
        
    # create folder/path to store assembled pairs
    dire = "{}/{}/{}".format(path, level, old_pair1)
    pathlib2.Path(dire).mkdir(parents=True, exist_ok=True)
 
    # promotion of unpaired fragment to next higher hierarchy level
    if promotion and m_i == (len(fragment_l)-1):
        print('promotion in level, pair', level, old_pair1)
        if not os.path.exists('{}/{}'.format(dire, top)):
            # print('already copied promoted fragment ', 
            #       level, previous_level, old_pair1)

            shutil.copyfile('{}/{}'.format(old_dire1, top),
                                '{}/pair0.pdb'.format(dire))
            shutil.copyfile('{}/{}'.format(old_dire1, xtc),
                                '{}/pair.xtc'.format(dire))
    else:
        # load universe -> load conformations of fragment 1 and 2
        u1 = mda.Universe('{}/{}'.format(old_dire1, top),
                          '{}/{}'.format(old_dire1, xtc))
        u2 = mda.Universe('{}/{}'.format(old_dire2, top),
                          '{}/{}'.format(old_dire2, xtc))
 
        # residue ovearlap MD fragment
        if len(flatten(pair_l[1])) == 1:
            o = overlaps_d[old_pair2]
        # overlap grown pairs, always == overlaps[0]
        # -> overlap that differs "corrected" for when growing pairs with MD fragment
        else:
            o = overlap
 
        # get indices for assembly
        index_aln_l, index_clash_l, index_merge_l = get_residue_indices_for_assembly(
                                                overlap0=overlap, current_overlap=o,
                                                capping_groups=capping_groups, last_level=last_level,
                                                verbose=verbose)
 
        ## check if for u2 residue 2 that is aligned is a proline
        ## only makes a difference if residue overlap is < 2
        res2_u2 = u2.select_atoms('resid {}'.format(u2.atoms.residues[index_aln_l[-1]].resid))
        res2_u2_name = res2_u2.residues.resnames
        if res2_u2_name == 'PRO':
            proline_2nd_posi = True
 
        # create dictionary with desired residues 
        # defined as "mobile" and "ref" you want to superimpose
        select = translate_concept(u1, u2, proline_2nd_posi, *index_aln_l)
        if verbose:
            print('level pair to grow, previous level, old_pair1, old_pair2, overlap',
                  level, previous_level, old_pair1, old_pair2, o)
            print('fragment 1 ',
                  u1.select_atoms('{}'.format(select['mobile'])).residues.resnames)
            print('fragment 2 ',
                  u2.select_atoms('{}'.format(select['reference'])).residues.resnames)
 
        # assemble the fragments into pairs
        # (or pairs into pairs of pairs)
        if  draw_indices:
            rs = fragment_assembly(u1, u2, dire, select, index_clash_l, index_merge_l,
                              rmsd_cut_off, clash_distance,  kmax=k_max)
            ###############
            #r_l.append(rs)
            ###############
            return rs
        else:
            rs = r_l[m_i]
            fragment_assembly(u1, u2, dire, select, index_clash_l, index_merge_l,
                              rmsd_cut_off, clash_distance,  kmax=k_max, ri_l=rs,
                              draw_indices=draw_indices)
            return None
 
 
def reweighted_hierarchical_chain_growth(hcg_l, promo_l, overlaps_d, path0, path, kmax, 
             rmsd_cut_off=0.6, clash_distance=2.0, capping_groups=True,
             ri_l=None,  path2weights='weights/', theta=10.0,  verbose=False):
    """place holder - RHCG in development for multiprocessing"""
    return None