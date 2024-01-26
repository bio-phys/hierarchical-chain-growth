#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from chain_growth.hcg_list import make_hcl_l
from chain_growth.fragment_list import generate_fragment_list, get_sequence #  _from_fasta, get_sequence_from_pdb, 
from chain_growth.hcg_fct import hierarchical_chain_growth
import os , pytest, pickle 
import numpy as np
import MDAnalysis as mda
test_dir = os.path.dirname(os.path.abspath(__file__))

def run_hcg():
    """
    run hcg with ri_l
    """
    ################
    ## prepare HCG
    ################
    
    ## input file and path / output path
    # path to MD fragments
    path0 = os.path.join(test_dir, '../../examples/')
    # path to store assembled models in
    path = os.path.join(test_dir, 'test_run_truncated_tauK18/')
    # file with sequence, format: "fasta" or "PDB"
    sequence_f = os.path.join(test_dir, '../../examples/truncated_tauK18.fasta')
    
    ## fragment construction
    # length of MD fragments (without the end-capping groups if present)
    fragment_length = 5
    # length of the residue overlap between subsequent fragments
    # == number of the same residues in subsequent fragments
    overlap = 2
    
    ## generate list of fragments, dictionary of overlaps between fragments
    # generating overlaps_d is necessary, since to match the full-ltngh sequennce
    # the overlap between e.g., the two last fragments can vary
    fragment_l, overlaps_d = generate_fragment_list(sequence_f, fragment_length, overlap)
    
    ## lists for the HCG
    n_pairs = len(fragment_l)
    # hcg_l : list of paired fragments
    # promo_l : list to evaluate if last fragment of level m in hcg_l is promoted to level m+1
    hcg_l, promo_l = make_hcl_l(n_pairs)
    
    # maximal number of pairs/full-length models to assemble
    kmax = 100
    # MD fragments are sampled with or without end-capping groups
    capping_groups = True
    
    
    ###########
    ## run HCG with ri_l
    ###########
    ri_l = pickle.load(open(
        os.path.join(test_dir, '../../examples/run_chain_growth/truncated_tauK18/confIndex_all_levels.pkl'),
        'rb')
    )
    hierarchical_chain_growth(hcg_l, promo_l, overlaps_d, path0, path, kmax=kmax,
            capping_groups=capping_groups, ri_l=ri_l) #, verbose=True) 
 
def calc_radius_of_gyration(pdb_xtc_l):
    pdb = pdb_xtc_l[0]
    xtc = pdb_xtc_l[1]
    u = mda.Universe(pdb, xtc)
    rg_l = [u.atoms.radius_of_gyration() for ts in u.trajectory]
    return np.round(np.average(rg_l), 1)



def calc_end_to_end_distance(pdb_xtc_l):
    pdb = pdb_xtc_l[0]
    xtc = pdb_xtc_l[1]
    u = mda.Universe(pdb, xtc)
    e2e_l = []
    for ts in u.trajectory:
        N = u.select_atoms("name N")[0]
        C = u.select_atoms("name C")[-1]

        r = N.position - C.position
        d = np.linalg.norm(r)
        e2e_l.append(d)
    return np.round(np.average(e2e_l), 1)



def test_sequence():
    """
    check_sequence
    --------------
    script to check wether full-length IDP grown via HCG has the correct sequence
    == same as the input sequence
    """
    run_hcg()    
    # path to topology file of idp grown with (r)hcg
    
    file_hcg =  os.path.join(test_dir, 'test_run_truncated_tauK18/4/0/pair0.pdb')
    sequence_hcg = get_sequence(file_hcg, NA=False).get_sequence_list()

    
    # path to reference sequence == inout sequence
    file_ref = os.path.join(test_dir, '../../examples/truncated_tauK18.fasta')
    ## get amino acid sequence from PDB or fasta file in three letter code
    sequence_ref = get_sequence(file_ref, NA=False).get_sequence_list()
        
    assert (sequence_hcg == sequence_ref)


@pytest.mark.parametrize("pdb_xtc_l, target",
    [
            ([ os.path.join(test_dir, 'test_run_truncated_tauK18/4/0/pair0.pdb'),
               os.path.join(test_dir, 'test_run_truncated_tauK18/4/0/pair.xtc')], 20.6),
            ])
def test_run_hcg_RG(pdb_xtc_l, target):
    """
    test output ensemble for global properties 
    """
    R_G = calc_radius_of_gyration(pdb_xtc_l)
    assert R_G == target
    
    
@pytest.mark.parametrize("pdb_xtc_l, target",
    [
            ([ os.path.join(test_dir, 'test_run_truncated_tauK18/4/0/pair0.pdb'),
               os.path.join(test_dir, 'test_run_truncated_tauK18/4/0/pair.xtc')], 51.0),
            ])
def test_run_hcg_end2end(pdb_xtc_l, target):
    end_to_end_distance = calc_end_to_end_distance(pdb_xtc_l)
    assert  end_to_end_distance == target

    
