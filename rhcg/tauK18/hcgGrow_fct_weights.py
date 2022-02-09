#!/usr/bin/env python3

import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import align
import pathlib2, shutil, os
import MDAnalysis.analysis.distances as distances
from hcgList import flatten


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
    # (problem with finding correct python c++ libraries for KD neighbour search)
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
        print('overlap in indices fct: ' , current_overlap,
              '1 alnB alnE, 2 alnB alnE ', index_aln_l)
        print('1 clashB, 2 clashE ', index_clash_l)
        print('1 mergeB mergeE, 2 mergeB mergeE ', index_merge_l)
        
    return index_aln_l, index_clash_l, index_merge_l
        
def fragment_assembly(u1, u2, dire, select, index_clash_l, index_merge_l,
         rmsd_cut_off, clash_distance, kmax, w_l, ri_l=None, draw_indices=False,
         chain_weights_prev_l=None):
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
    w_l : list or array-like
        list with weight arrays per for the fragments. The weights define how likely it is to observe
        configuration a (frame a) in comparison to the native ensemble. These weights are
        used to perform weighted drawing, i.e., draw weighted random frames.
        Weights per frame come either from ensemble reweighting (in level m == 1)
        or can be uniform (in level m > 1).
    ri_l : array-like
        array with indices for chosing a specific confoormation of a fragment. 
        if None: draw indices randomly
    draw_indices : booolean
        if new random integers == frame indices are drawn or else taken from a input array.
    chain_weights_prev_l : list or array-like
        array of unnormalized weight products from fragments assembled in prevoius level MODIFY!!!
        eq ??? in stelzl at al JACS Au 2022
        DO NOT USE as fragment weight to draw random frame - not yet normalized!
        
        
    Returns
    -------
    if draw_indices:
        rs : list of lists
            successful indices drawn for fragment 1 and 2 during fragment assembly
        assembled_chain_weights : array
        
    else:
        None
    """
    
    k = 0
    assembly_atempt = 0
    assembled_chain_weights = np.zeros(kmax)
    
    if np.all(ri_l) is None:
        
        rs = np.zeros((kmax, 2))
    else:
        print(ri_l.shape)
    writePDB = True
    while k < kmax:
        if draw_indices:
            # random integer to draw random frame
            r1 = np.random.choice(u1.trajectory.n_frames, p=w_l[0])
            r2 = np.random.choice(u2.trajectory.n_frames, p=w_l[1])
        else:
            # r1 = int(ri_l[k][0])
            # r2 = int(ri_l[k][1])
            r1 = int(ri_l[0,k])
            r2 = int(ri_l[1,k])
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
      #      print(r1, r2, clashes)
            
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
                
                ## calculate weight for assembled chain/ fragment at step k
                ## unnormalized! except for mdfragments
                # print(chain_weights_prev_l[0].shape, chain_weights_prev_l[1].shape, len(chain_weights_prev_l))#,
                      # chain_weights_prev_l[0][r1] , chain_weights_prev_l[1][r2])
                # if np.all(chain_weights_prev_l) is not None:
                    
                  #  print(chain_weights_prev_l[0][r1], chain_weights_prev_l[1][r2])
                assembled_chain_weights[k] = chain_weights_prev_l[0][r1] * chain_weights_prev_l[1][r2]
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
        return rs, assembled_chain_weights
    else:
        return assembled_chain_weights


def hierarchical_chain_growth(hcg_l, promo_l, overlaps_d, path0, path, kmax, 
             rmsd_cut_off=0.6, clash_distance=2.0, capping_groups=True,
             ri_l=None,  path2weights='weights/', theta=10.0,  verbose=False):
    """ perform hierarchical chain growth + ímportance sampling
    assemble fragments (reweighted according to experimental data)
                        or pairs of fragments until reaching the full-length chain

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
    k_max = kmax
    
    for m , fragment_l in enumerate(hcg_l): 
        # folder to save assembled pair in this level (m+1)
        level = m+1
        # folder to get old pairs from previous level (m)
        previous_level = m
        promotion  = promo_l[m]
        overlap = overlaps_d[0]
        
        if level > 1 :# and weighted:
            chain_weights_prev = np.load("{}/chain_weight_level{}.npy".format(path, previous_level), allow_pickle=True)
        else:
            chain_weights_prev = None 

        
        if np.all(ri_l) is None:
            draw_indices = True
            r_l = []
            
        else:
            r_l = ri_l[m]
        # draw_indices = False
        # r_l = ri_l[m]
        assembled_chain_weights = []#np.full((len(fragment_l), kmax), fill_value = np.nan)    
        # if MD fragments are sampled with end-capping_groups, 
        # they are removed in the last assembly step 
        if m+1 == len(hcg_l):
            last_level = True
            k_max = kmax

        # in case one wants to grow only very few full length models
        # grow at least 20 fragments for each pair in lower levels
        # to make sure one still succeeds
        if not last_level and kmax < 20:
            k_max = 20  # more?
        c1 = 0
        for m_i , pair_l in enumerate(fragment_l):
            c2=c1+1
            proline_2nd_posi = False
            # if promotion of MD fragment in first level DO NOT define old_pair2
            # index ERROR because in this case, pair_l is no list
            if promotion and isinstance(pair_l, list) == False:
                # subfolder old fragment 1, 2
                old_pair1  = pair_l
                old_pair2  = pair_l
            else:
                # subfolder old fragment 1, 2
                old_pair1 = flatten(pair_l[0])[0]
                old_pair2 = flatten(pair_l[1])[0]

            if m == 0:
                previous_level = 'MDfragments'             
                path2fragment = path0
            else:
                path2fragment = path 
            old_dire1 = '{}/{}/{}'.format(path2fragment, previous_level, old_pair1)
            old_dire2 = '{}/{}/{}'.format(path2fragment, previous_level, old_pair2)
            # create folder/path to store assembled pairs
            dire = "{}/{}/{}".format(path, level, old_pair1)
            pathlib2.Path(dire).mkdir(parents=True, exist_ok=True)
            
            # print('level, fragment, old_pairs, counter ', level, m_i, old_pair1, old_pair2, c1, c2)#, chain_weights_prev)
            
            # promotion of unpaired fragment to next higher hierarchy level
            if promotion and m_i == (len(fragment_l)-1):
                print('promotion in level, pair',  level, old_pair1)
                
                
                # store weight of promoted fragment for next level
                if level == 1:
                    w  = np.genfromtxt('{}/{}/wopt/w_theta{}_dat.txt'.format(path2weights, old_pair1, theta))
                else:
                    # print('level, fragment, old_pairs, counter, shape chain_weights_prev ', level, m_i, old_pair1,
                    #   old_pair2, c1, c2,  chain_weights_prev.shape)
                    w = chain_weights_prev[c1]
                # assembled_chain_weights[m_i] = w
                assembled_chain_weights.append(w)
                if os.path.exists('{}/pair0.pdb'.format(dire)):
                    # print('already copied promoted fragment ', 
                    #       level, previous_level, old_pair1)
                    continue
                else:
                    shutil.copyfile('{}/pair0.pdb'.format(old_dire1),
                                        '{}/pair0.pdb'.format(dire))
                    shutil.copyfile('{}/pair.xtc'.format(old_dire1),
                                        '{}/pair.xtc'.format(dire))
                    continue
                
            # load universe -> load conformations of fragment 1 and 2
            u1 = mda.Universe('{}/pair0.pdb'.format(old_dire1),
                              '{}/pair.xtc'.format(old_dire1))
            u2 = mda.Universe('{}/pair0.pdb'.format(old_dire2),
                              '{}/pair.xtc'.format(old_dire2))                       
            ## weights 
            if level == 1:
                
                w1  = np.genfromtxt('{}/{}/wopt/w_theta{}_dat.txt'.format(path2weights, old_pair1, theta))
                w2  = np.genfromtxt('{}/{}/wopt/w_theta{}_dat.txt'.format(path2weights, old_pair2, theta))
                chain_weights_prev_l = [w1, w2]
            else:
                # print('level, fragment, old_pairs, counter, shape chain_weights_prev ', level, m_i, old_pair1,
                #       old_pair2, c1, c2, chain_weights_prev.shape)
                w1 = np.ones(u1.trajectory.n_frames)/u1.trajectory.n_frames
                w2 = np.ones(u2.trajectory.n_frames)/u2.trajectory.n_frames
                chain_weights_prev_l = [chain_weights_prev[c1], chain_weights_prev[c2]]
            
            # store w1, w2 in a list as input for fragment_assembly
            w_l = [w1, w2]
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
            
            ## check if for u2 residue 2 that is aligned is ja proline
            res2_u2 = u2.select_atoms('resid {}'.format(u2.atoms.residues[index_aln_l[-1]].resid))
            res2_u2_name = res2_u2.residues.resnames
            if res2_u2_name == 'PRO':
                proline_2nd_posi = True
            
            # create dictionary with desired residues 
            # defined as "mobile" and "ref" you want to superimpose
            select = translate_concept(u1, u2, proline_2nd_posi, *index_aln_l)
            if verbose:
                print('level pair to grow, length old pair, old_pair1, old_pair2, overlap',
                      level, previous_level, old_pair1, old_pair2, o)
                print('fragment 1 ',
                      u1.select_atoms('{}'.format(select['mobile'])).residues.resnames)
                print('fragment 2 ',
                      u2.select_atoms('{}'.format(select['reference'])).residues.resnames)
            
            # assemble the fragments into pairs
            # (or pairs into pairs of pairs)
            if  draw_indices:
                
                rs, acw_mi = fragment_assembly(u1, u2, dire, select, index_clash_l, index_merge_l,
                                  rmsd_cut_off, clash_distance,  kmax=k_max, w_l=w_l, ri_l=ri_l,
                                  draw_indices=draw_indices, chain_weights_prev_l=chain_weights_prev_l)
                r_l.append(rs)
                # assembled_chain_weights[m_i] = acw_mi
                assembled_chain_weights.append(acw_mi)
            else:
                rs = r_l[m_i]
                # rs = r_l[:,m_i,:,]
                acw_mi = fragment_assembly(u1, u2, dire, select, index_clash_l, index_merge_l,
                                  rmsd_cut_off, clash_distance,  kmax=k_max, w_l=w_l, ri_l=rs,
                                  draw_indices=draw_indices, chain_weights_prev_l=chain_weights_prev_l)
                assembled_chain_weights.append(acw_mi)
            c1 += 2
        if draw_indices:
            # print(draw_indices)
            np.save("{}/confIndex_level{}.npy".format(path, level), r_l)
        # assembled_chain_weights = np.asarray(assembled_chain_weights, dtype=object)
        # print(assembled_chain_weights)
        np.save("{}/chain_weight_level{}.npy".format(path, level), assembled_chain_weights)
    return None
