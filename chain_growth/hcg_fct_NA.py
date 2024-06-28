#!/usr/bin/env python3
"""
hcg_fct
-----------
core functions for hierarchical chain growth

run hcg in parallel - per level -> m_i loop is parallelized
using multipprocessing.Pool
for nucleic acids! principle is the same as for IDPs, minor changes for alignment and assembly
"""
import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import align
import pathlib2, shutil, os, itertools
import MDAnalysis.analysis.distances as distances
from chain_growth.hcg_list import flatten, make_hcl_l
from chain_growth.fragment_list import generate_fragment_list
from chain_growth.hcg_fct import get_residue_indices_for_assembly, merge_universe
from multiprocessing import Pool
from functools import partial


def translate_concept_NA(u1, u2, align_begin1, align_end1, 
                      align_begin2, align_end2):
    """ prepare a dictionary with atoms to align before the assembly step
    
    Parameter
    ---------
    u1 : universe
    u2 : universe
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
    
    select={}
    ## atom selection of fragment1 to be aligned
    ## atom selection of fragment2 to be aligned

 
    sel1 = "nucleic and not type H and ((resid {} and name [   O3'  ]) or (resid {} and not name [  O1P , O2P , O2' , C2' , C3' , O3' , C4' , C5' , O4' ] ))".format(
           u1.atoms.residues[align_begin1].resid,
           u1.atoms.residues[align_end1].resid)
    sel2 = "nucleic and not type H and ((resid {} and name [  O3'  ]) or (resid {} and not name [  O1P , O2P , O2' , C2' , C3' , O3' , C4' , C5' , O4' ] ))".format(
                        u2.atoms.residues[align_begin2].resid,
                        u2.atoms.residues[align_end2].resid)
 
    ## dict is argument for mda alignto function: arg = "select="
    ## mobile is superimposed on ref with mda.alignto in the assembly step
    select={'mobile': sel1, 'reference': sel2}
    
    return select

def find_clashes_NA(u1,u2, index1b, index2e, 
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
    index1e : integer
        residue index for residues to *exclude* from clash detection
        marks residue to end exclusion for u1. The default is -1
    index2e : integer
        residue index for residues to *exclude* from clash detection
        marks residue to end exclusion for u2
    clash_distance : float
        max. allowed distance between atoms, everything below is counted as clash
        The  default is 2.0
   
    Returns
    -------
    clsum : integer
        count of clashes between aligned fragments
        
    NOTE: MDAnalysis is inclusive! resid 1:2 -> selects residues 1+2
    """
   
    l1 = u1.select_atoms("all and not (type H) and not " 
                         "( resid {}:{} )".format(                            
                             u1.atoms.residues[index1b].resid, u1.atoms.residues[index1e].resid)) 
    
    l2 = u2.select_atoms("all and not (type H) and not " 
        "resid 1:{}  ".format(u2.atoms.residues[index2e].resid))  
    
    ## -> use mda.distances to generate a matrix with distances
    distmat = distances.distance_array(l1.positions,l2.positions)
    cont = np.less(distmat,clash_radius) # see where distmat < cutoff
    cont = np.where(cont,1,0) # True --> 1, False --> 0
    #print(distmat, cont)
    clsum = np.sum(cont)
    
    return clsum


def fragment_assembly_NA(u1, u2, dire, select, index_clash_l, index_merge_l, w_aln, kmax, 
         rmsd_cut_off, clash_distance, w, ri_l=None, verbose=False):
    """ assemble the fragments to pairs

    Parameters
    ----------
    u1 : universe
        fragment 1 to pair
    u2 : universe
        fragment 2 to pair
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
    kmax : integer
        number of pairs that should be assembled in level m_i
    rmsd_cut_off : float
        cut-off for the RMSD of the fragment alignment
    clash_distance : float
        minimal allowed distance between atoms
    w : list
        list of weights for fragment 1 and 2 .
        Use uniforrm weights.
    ri_l : array-like
        array with indices for chosing a specific confoormation of a fragment. The default is None
        if None: draw indices randomly

    Returns
    -------
    None.
    """

    find_clashes = find_clashes_NA
    k = 0
    assembly_atempt = 0
    clashReject = 0
    rmsdReject = 0
    writePDB = True
    if np.all(ri_l) is None:
        rs = np.zeros((kmax, 2))
 

    while k < kmax:
        # random integer to draw random frame
        if np.all(ri_l) is None:
            r1 = np.random.choice(u1.trajectory.n_frames, p=w[0])
            r2 = np.random.choice(u2.trajectory.n_frames, p=w[1])
        else:
            r1 = int(ri_l[k][0])
            r2 = int(ri_l[k][1])
        # load random frame
        u1.trajectory[r1]
        u2.trajectory[r2]
    
        # rmsd suportimposition of the subsequent fragmenmt         )
        old,new = align.alignto(u1, u2, select=select, weights=w_aln)

        assembly_atempt += 1
        if new < rmsd_cut_off:
            # calculate clashes
            clashes = find_clashes(u1 , u2 , index1b=index_clash_l[0] , index2e=index_clash_l[1] ,
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
               
            
                if np.all(ri_l) is None:
                    rs[k] = [r1,r2]

                k+=1
            else:
                clashReject +=1
                continue
        else:
            rmsdReject +=1
            continue
        
            
    # load and save universe with pair0 as topology and atom positions 
    # of subsequent paired fragment of this assembly step as trajectory
    # n frames = desired number of pairs in this assembly step
    allPairs = mda.Universe("{}/pair0.pdb".format(dire) , at)
    allPairs.atoms.write("{}/pair.xtc".format(dire) , frames='all')
    if np.all(ri_l) is None:
        return rs     #, rmsdReject, clashReject, assembly_atempt
    else:
        return None


def hierarchical_chain_growth_NA(hcg_l, promo_l, overlaps_d, path0, path, kmax, 
            dict_to_fragment_folder=None, rmsd_cut_off=0.6, clash_distance=2.0, capping_groups=True,
            ri_l=None, verbose=False):
    """ perform hierarchical chain growth 
    assemble fragments/ pairs of fragments until reaching the full-length chain

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
    kmax : integer
        number of pairs that should be assembled in level m_i
    dict_to_fragment_folder  : dictionary
        dictonary that interpretes / translates primary sequence and folder atructure in the fragment lirary
    rmsd_cut_off : float, optional
        cut-off for the RMSD of the fragment alignment. The default is 0.6
    clash_distance : float, optional
        max. allowed distance between atoms. The default is 2.0
    capping_groups : boolean, optional
        MD fragment are sampled with or without end-capping groups. The default is True
    ri_l : array-like
        array with indices for chosing a specific confoormation of a fragment. The default is None
    
    Returns
    -------
    None.
    """
    
    last_level = False
    cpu = os.cpu_count()
    k_max = kmax
    ###
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
            # print(ri_l)
        ## if MD fragments are sampled with end-capping_groups, 
        ## they are removed in the last assembly step 
        if level == len(hcg_l):# and capping_groups:
            last_level = True
            k_max = int(kmax/2)
            # print('last level')

        if len(fragment_l) > cpu:
            num_threads = cpu
        else:
            num_threads = len(fragment_l)
        ## construct pair list for level_i
        pairs = [(m_i, pair_l) for m_i, pair_l in enumerate(fragment_l)]
        ## dictionary with variables for _loop_func that performs HCG in parallel for each pair_i in level_i
        d = {"path0": path0, "path": path, "level": level, 'previous_level': previous_level,
             "last_level": last_level, "fragment_l": fragment_l, "rmsd_cutoff": rmsd_cut_off,
             "dict_to_fragment_folder": dict_to_fragment_folder, "clash_distance": clash_distance, "kmax": k_max, "r_l":r_l,
             "overlap_d": overlaps_d, "promotion": promotion, "capping_groups": capping_groups,
             "draw_indices" : draw_indices,  "verbose": verbose }
        # POOL LOOP application
        with Pool(num_threads) as p:
            func = partial(_loop_func, d)
            results = p.map(func, pairs)
        if draw_indices:
            for r in results:
                if r is not None:
                    r_l += [r]
            ## save indices of successfully paired fragments in level_i
            np.save("{}/confIndex_level{}.npy".format(path, m+1), r_l) 
        ###
    return None


def _loop_func(variables, pairs):
    """ does the inner loop of hierarchical chain growth and calls fragment assembly

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
    level = variables["level"]
    previous_level = variables["previous_level"]
    last_level = variables["last_level"]
    fragment_l = variables["fragment_l"]
    rmsd_cut_off = variables["rmsd_cutoff"]
    clash_distance = variables["clash_distance"]
    k_max =  variables["kmax"]
    dict_to_fragment_folder = variables["dict_to_fragment_folder"]
    r_l = variables["r_l"]
    draw_indices = variables["draw_indices"]
    overlaps_d = variables["overlap_d"]
    promotion = variables["promotion"]
    capping_groups = variables["capping_groups"]
    verbose = variables["verbose"]
    overlap = overlaps_d[0]


    ## folder to get old pairs from == length of the pair assembled in m-1 (previous level)
    ## folder to save assembled pair in == length of the pair assembled in level m
    ## -> total number of fragments assembled
    ## without taking into account pairs assembled from previously promoted fragments
    ## except from last level
    if promotion and isinstance(pair_l, list) == False:
        ## subfolder fragment 1 is stored in is called 
        old_pair1  = pair_l
        old_pair2  = pair_l
    else:
        ## subfolder fragment 1 is stored in is called 
        old_pair1 = flatten(pair_l[0])[0]
        old_pair2 = flatten(pair_l[1])[0]           
    
    if previous_level == 0 and dict_to_fragment_folder is not None:
        # print(dict_to_fragment_folder)
        path2fragment = path0 
        old_pair1_fragmentLib = dict_to_fragment_folder[old_pair1]
        old_pair2_fragmentLib = dict_to_fragment_folder[old_pair2]
        if m_i == 0:
            old_dire1 = '{}/5p_terminus/{}'.format(path2fragment, old_pair1_fragmentLib)
        else:
            old_dire1 = '{}/{}'.format(path2fragment, old_pair1_fragmentLib)
        old_dire2 = '{}/{}'.format(path2fragment, old_pair2_fragmentLib)
        top = 'fragment.pdb'
        xtc = 'fragment.xtc'    
    elif previous_level == 0 and dict_to_fragment_folder is None:
        previous_level = 1             
        path2fragment = path0
        old_dire1 = '{}/{}'.format(path2fragment, old_pair1)
        old_dire2 = '{}/{}'.format(path2fragment, old_pair2)
        top = 'pair0.pdb'
        xtc = 'pair.xtc'
    ## levels > 1
    else:
        path2fragment = path
        old_dire1 = '{}/{}/{}'.format(path2fragment, previous_level, old_pair1)
        old_dire2 = '{}/{}/{}'.format(path2fragment, previous_level, old_pair2)
        top = 'pair0.pdb'
        xtc = 'pair.xtc'

    ## create folder/path to store assembled pairs
    dire = "{}/{}/{}".format(path, level, old_pair1)
    pathlib2.Path(dire).mkdir(parents=True, exist_ok=True)
         
   
    ## promotion of unpaired fragment to next higher hierarchy level
    if promotion and m_i == (len(fragment_l)-1):
        if not os.path.exists('{}/{}'.format(dire, top)):

            # print('need to copy promoted fragment', 
            #       length_old_pair, old_pair1, '\n')
            shutil.copyfile('{}/{}'.format(old_dire1, top),
                                '{}/pair0.pdb'.format(dire))
            shutil.copyfile('{}/{}'.format(old_dire1, xtc),
                                '{}/pair.xtc'.format(dire))
        
    else:
        ## load universe -> load conformations of fragment 1 and 2
        u1 = mda.Universe('{}/{}'.format(old_dire1, top),
                          '{}/{}'.format(old_dire1, xtc))
        u2 = mda.Universe('{}/{}'.format(old_dire2, top),
                          '{}/{}'.format(old_dire2, xtc))
    
        ## by default load uniform weights
        w1 = np.ones(u1.trajectory.n_frames)/u1.trajectory.n_frames
        w2 = np.ones(u2.trajectory.n_frames)/u2.trajectory.n_frames
        
        ## ovearlap MD fragment
        if len(flatten(pair_l[1])) == 1:
            o = overlaps_d[old_pair2]
               
            
        ## overlap grown pairs, always == overlaps[0]
        ## -> overlap that differs "corrected" for when growing pairs with MD fragment
        else:
            o = overlap
        
        index_aln_l, index_clash_l, index_merge_l = get_residue_indices_for_assembly(
                                                current_overlap=o, last_level=last_level, overlap0=overlap,
                                                capping_groups=capping_groups, verbose=verbose)
        ## create dictionary with desired residues defined as "mobile" and "ref" you want to superimpose
        select = translate_concept_NA(u1, u2, *index_aln_l)
        ## assign weigths to atoms to be aligned
        l = u1.select_atoms('{}'.format(select['mobile'])).__len__()
        y1 = np.array([2.5, 2.5, 2.5])
        y2 = np.repeat(0.8, l-3)
        y = np.concatenate([y1, y2])
        w_aln = y / y.sum()
        if verbose:
            print('level pair to grow, previous_level, old_pair1, old_pair2, overlap',
                  level, previous_level, old_pair1, old_pair2, o)
            print('fragment 1 ',
                  u1.select_atoms('{}'.format(select['mobile'])).residues.resnames)
            print('fragment 2 ',
                  u2.select_atoms('{}'.format(select['reference'])).residues.resnames)
    
    
        if draw_indices:
            rs  = fragment_assembly_NA(u1, u2, dire, select, index_clash_l, index_merge_l, w_aln=w_aln,
                            rmsd_cut_off=rmsd_cut_off, clash_distance=clash_distance,
                            kmax=k_max, w=[w1,w2], verbose=verbose)
    
            return rs
           
        else:
            rs = r_l[m_i]
            fragment_assembly_NA(u1, u2, dire, select, index_clash_l, index_merge_l, w_aln=w_aln,
                              rmsd_cut_off=rmsd_cut_off, clash_distance=clash_distance, 
                              kmax=k_max, w=[w1,w2], 
                              ri_l=rs, verbose=verbose)
            return None

def run_hcg_NA(sequence, kmax, path0='fragmentLib_RNA/' , path='out/',
                         fragment_length=3, overlap=1, capping_groups=False,
                         clash_distance = 2.0, online_fragment_library=True, rna=True,
                         rmsd_cut_off=0.64, ri_l=None,
                         verbose=False):
    """ Here we process your input sequence for the run with HCG. 
        Other default arguments are defined here that do not to be user defined. 
        We give some explanation of the default arguments. 
        In order to just run HCG and build models it is not mandatory to read this.
   
    Parameters
    ----------
    sequence : string 
        sequence of the IDP to grorw. Can be a string with the nucleic acid sequence,
        the path to a fasta file or a PDB file. Note: Function is NOT case sensitive.
        Only the 4 natural nucleic acid bases are recognized.
    kmax : integer
        Max. number of full-length polymers that should be grown.
    path0 : string
        path to input fragment library. Default is presampled "fragmentLib_RNA/"
        The library can be found on Zenodo: https://zenodo.org/records/8369324
    path : string
        path to output folder the assembled fragments are stored in - folder is created here.
        Default is "out/"
    fragment_length : integer
        length of the fragment. Default is 3 for growing ssRNA from presampled fragmentLib_RNA.
    overlap : integer
        length of the residue overlap between subsequent fragments. Default is 1. 
        Please check for more detailed vdescription of the fragment assembly procedure
    capping_groups : boolean
        chemical capping groups attached to MD fragments. Default is False
    clash_distance : float
        max. allowed distance between atoms, everything below is counted as clash
        The  default is 2.0
    rmsd_cut_off : float, optional
        cut-off for the RMSD of the fragment alignment. The default is 0.64
    online_fragment_library : boolean
        whether polymer is grown from the predsampled MD dimer library or your own library.
        Note: if you  decide to use your own library you need to modify the dict_to_fragment_folder to match your folder structure.
        Default is True
    rna : boolean
        whether single stranded nucleic acid polymer to grow is RNA or DNA. To run with DNA use own fragment library - alignment and clash functions work for DNA, too
        Default is True
        
    Returns
    -------
    None.
    """

    ## generate list of fragments, dictionary of overlaps between fragments
    ## generating overlaps_d is necessary to match the full-length sequence
    fragment_l, overlaps_d = generate_fragment_list(sequence, fragment_length, overlap, NA=True)
    ## number of fragments
    n_pairs = fragment_l.__len__()

    ## hcg_l : list of paired fragments
    ## promo_l : list to evaluate if last fragment of level m in hcg_l is promoted to level m+1
    hcg_l, promo_l = make_hcl_l(n_pairs)

    if online_fragment_library is True: 
        ## what kind of nucleic acid to grow
        ## to run with DNA use own fragment library
        ## alignment and clash functions work for DNA, too
        if rna:  
            aa_l = ['A', 'C', 'G', 'U']
        else:
            aa_l = ['A', 'C', 'G', 'T']

        # create dictionary to translate between the sequence information and folder structure in fragment library
        d = { p : i for i, p in enumerate(itertools.product(aa_l, repeat=fragment_length))}
        ## using this dict we can directly communicate between the order of fragments (== fragment id)
        ## as defined by the sequence and the respective folder_id
        dict_to_fragment_folder = {i : d[tuple(aa_pair)] for i , aa_pair in enumerate(fragment_l)}

    hierarchical_chain_growth_NA(hcg_l, promo_l, overlaps_d, path0, path, kmax,
                  dict_to_fragment_folder=dict_to_fragment_folder,
                  rmsd_cut_off=rmsd_cut_off, clash_distance=clash_distance,
                  capping_groups=capping_groups, ri_l=ri_l,
                  verbose=verbose)

    return None 
