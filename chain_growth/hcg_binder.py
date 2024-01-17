#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
hcg_binder
-----------
functions for calling HCG on a binder instancez

@author: Lisa
"""
import numpy as np
import concurrent.futures
import functools , time , itertools , warnings
warnings.filterwarnings('ignore')
from chain_growth.fragment_list import generate_fragment_list
from chain_growth.hcg_list import make_hcl_l
from chain_growth.hcg_fct import hierarchical_chain_growth

from tqdm.auto import tqdm


def progress_bar(expected_time, increments=10):

    #Here we create a progresss bar for HCG as decorator. Please run the collapsed cells.

    def _progress_bar(func):

        def timed_progress_bar(future, expected_time, increments=10):
            """
            https://stackoverflow.com/questions/59013308/python-progress-bar-for-non-loop-function
            
            Display progress bar for expected_time seconds.
            Complete early if future completes.
            Wait for future if it doesn't complete in expected_time.
            """
            interval = expected_time / increments

            with tqdm(total=increments) as pbar:
                for i in range(increments - 1):
                    if future.done():
                        # finish the progress bar
                        # not sure if there's a cleaner way to do this?
                        pbar.update(increments - i)
                        return
                    else:
                        time.sleep(interval)
                        pbar.update()
                # if the future still hasn't completed, wait for it.
                future.result()
                pbar.update()

        @functools.wraps(func)
        def _func(*args, **kwargs):
            with concurrent.futures.ThreadPoolExecutor(max_workers=1) as pool:
                future = pool.submit(func, *args, **kwargs)
                timed_progress_bar(future, expected_time, increments)

            return future.result()

        return _func

    return _progress_bar

def estimate_run_time(data):
    """ Estimate run time for HCG - depends on kmax = number of full-length chains to grow 
                                     and chain length = number amino acids per chain

        Run time is ~ proportional to kmax * chain length.
        Therefore we fit a linear function for both arguments and multiply them.
        Parameters for functions from a least square fit.
        NOTE: The estimated time may overestimate the true run time for several seconds to minutes - difference more pronounced for large kmax or chain lengths.
        
        Parameters
        ----------
        data : numpy array
            x = kmax
            y = chain length
            
        Returns
        -------
        estimated time as float
        
        Thanks to Johannes Betz.
        
        """
    x, y = data
    def _fit_linear(x, a, b, c):
        return a  + b*x + c*x *np.log(x)
    fit_x = _fit_linear(x, 1.719e-01,  1.131e-02, -1.689e-05)
    fit_y = _fit_linear(y, 1.236e+01, -7.738e-02,  2.485e-02)
    return fit_x * fit_y



def run_hcg_binder(sequence, kmax, path0='dimerLibrary/' , path='out/',
                         fragment_length=2, overlap=0, capping_groups=True,
                         clash_distance = 2.0, online_fragment_library=True,
                         rmsd_cut_off=0.6, ri_l=None, streamlit_progressbar=None,
                         verbose=False):
    """ Here we process your input sequence for the run with HCG. 
        Other default arguments are defined here that are not user defined. 
        We give some explanation of the default arguments. 
        In order to just run HCG and build models it is not mandatory to read this.
        
    Parameters
    ----------
    sequence : string 
        sequence of the IDP to grow. Can be a string with the sequence in one letter amino acid code,
        the path to a fasta file or a PDB file. Note: Function is NOT case sensitive.
        Only the 20 natural amino acids are recognized.
    kmax : integer
        Max. number of full-length polymers that should be grown.
    path0 : string
        path to input fragment library. Default is presampled "dimerLibrary/""
    path : string
        path to output folder the assembled fragments are stored in - folder is created here.
        Default is "out/"
    fragment_length : integer
        length of the fragment. Default is 2 for growing IDPs from presampled fragment library.
    overlap : integer
        length of the residue overlap between subsequent fragments. Default is 0        
    capping_groups : boolean
        chemical capping groups attached to MD fragments. Default is true
    clash_distance : float
        max. allowed distance between atoms, everything below is counted as clash
        The  default is 2.0
    rmsd_cut_off : float, optional
        cut-off for the RMSD of the fragment alignment. The default is 0.6
    online_fragment_library : boolean\
        whether polymer is grown from the predsampled MD dimer library or your own library.
        Note: if you  decide to use your own library you need to modify the dict_to_fragment_folder to match your folder structure.
        Default is True
    streamlit_progressbar : object
        progress bar mapping the progress of HCG, e.g., as used by streamlit (webapp)
        Default is None. 
            If None the progress is mapped by tqdm package.
            
    Returns
    -------
    None.
    """


    ## generate list of fragments, dictionary of overlaps between fragments
    ## generating overlaps_d is necessary to match the full-length sequence
    fragment_l, overlaps_d = generate_fragment_list(sequence, fragment_length, overlap)
    ## number of fragments
    n_pairs = fragment_l.__len__()

    ## hcg_l : list of paired fragments
    ## promo_l : list to evaluate if last fragment of level m in hcg_l is promoted to level m+1
    hcg_l, promo_l = make_hcl_l(n_pairs)

    if online_fragment_library is True: 
        # create dictionary to translate between the sequence information and folder structure in fragment library
        aa_l = ['GLY','ALA', 'VAL', 'LEU', 'ILE', 'THR', 'SER', 'CYS', 'GLN', 'ASN', 'GLU', 'ASP',
                'LYS', 'TRP', 'ARG', 'TYR', 'PHE', 'HIS', 'PRO', 'MET']
        d = { p : i for i, p in enumerate(itertools.product(aa_l, repeat=fragment_length))}
        ## using this dict we can directly communicate between the order of fragments (== fragment id)
        ## as defined by the sequence and the respective folder_id
        dict_to_fragment_folder = {i : d[tuple(aa_pair)] for i , aa_pair in enumerate(fragment_l)}


    # print(dict_to_fragment_folder)
    ## before we run HCG we determine the number of hierarchical levels. You need this number for the progress bar we will generate for the HCG run. 
    ## later we need this number, which is also the folder name the full-length polymer is stored in.
    number_hcg_levels = hcg_l.__len__()
    ## the expected time depends on the length of the protein to grow and kmax
    ## it would be good to haveTODO:  an estimate of this time from these values
    data = np.vstack([kmax, len(sequence)])
    expected_time = estimate_run_time(data)[0] 
    # print(expected_time)
    if streamlit_progressbar is not None:
        hierarchical_chain_growth(hcg_l, promo_l, overlaps_d, path0, path, kmax,
                  dict_to_fragment_folder=dict_to_fragment_folder,
                  rmsd_cut_off=rmsd_cut_off, clash_distance=clash_distance,
                  capping_groups=capping_groups, ri_l=ri_l,
                  streamlit_progressbar=streamlit_progressbar,
                  verbose=verbose)
    else:
        @progress_bar(expected_time=expected_time, increments=number_hcg_levels)
        def hcg(hcg_l, promo_l, overlaps_d, path0, path, kmax,
                dict_to_fragment_folder,
                rmsd_cut_off, clash_distance, capping_groups, ri_l,
                verbose):
            hierarchical_chain_growth(hcg_l, promo_l, overlaps_d, path0, path, kmax,
                dict_to_fragment_folder,
                rmsd_cut_off, clash_distance, capping_groups, ri_l,
                verbose)
            return None
    
        hcg(hcg_l, promo_l, overlaps_d, path0, path, kmax,
                  dict_to_fragment_folder=dict_to_fragment_folder,
                  rmsd_cut_off=rmsd_cut_off, clash_distance=clash_distance,
                  capping_groups=capping_groups, ri_l=ri_l, verbose=verbose)

    return None #number_hcg_levels, path
