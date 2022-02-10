import MDAnalysis as mda
import numpy as np
from MDAnalysis.analysis.rms import rmsd
import MDAnalysis.analysis.distances as distances


def mk_universe(path_gro , path_trj):
    u = mda.Universe(path_gro, path_trj)
    return u



def translate_concept(u1, u2,index1_alnB=-3, index1_alnE=-2, index2_alnB=1, index2_alnE=2 ):
    select={}
     # atom selection of growing peptide u1 to be aligned
    sel1 = "(resid {} and (name C or name O)) or (resid {} and (name N or name H))".format(u1.atoms.residues[index1_alnB].resid, 
                                             u1.atoms.residues[index1_alnE].resid) 
    # atom selection of new peptide u2 to be aligned
    sel2 = "(resid {}  and (name C or name O)) or (resid {} and (name N or name H))".format(u2.atoms.residues[index2_alnB].resid ,
                                               u2.atoms.residues[index2_alnE].resid)
   
    # dict is argument for mda alignto function: arg = "select="
    select={'mobile': sel1, 'reference': sel2}
    return select


def find_clashes(u1,u2, index1b=-3, index2=3, 
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
                                                    u2.atoms.residues[index2-1].resid,
                                                    u2.atoms.residues[index2].resid))
    # (problem with finding correct python c++ libraries for KD neighbour search)
    # -> use mda.distances to generate a matrix with distances
    ## code snippet from Soeren von Buelow
    distmat = distances.distance_array(l1.positions,l2.positions)
    cont = np.less(distmat,clash_radius) # see where distmat < cutoff
    cont = np.where(cont,1,0) # True --> 1, False --> 0
    clsum = np.sum(cont)
    
    return clsum

## clash search a bit less stringent
#def findClashes(u1,u2, index1b=-2, index1e=-1, index2=3, clash_radius=2.7):
#    # atom selection of first fragment
#    l1 = u1.select_atoms("protein and not (type H) and not resid {}:{}".format(u1.atoms.residues[index1b].resid , u1.atoms.residues[index1e].resid)) 
#    # atom selection of subsequent fragmenmt
#    l2 = u2.select_atoms("protein and not (type H) and not resid 1:{}".format(u2.atoms.residues[index2].resid))
#
#    ## code snippet from Soeren von Buelow
#    distmat = distances.distance_array(l1.positions,l2.positions)
#    cont = np.less(distmat,clash_radius) # see where distmat < cutoff
#    cont = np.where(cont,1,0) # True --> 1, False --> 0 (to count pairs)
#    clsum = np.sum(cont)
#    return clsum



def merge_universe(u1,u2, sel1, sel2):
    # atom selection of growing peptide u1 to be merged
    sela = u1.select_atoms(sel1)
    # atom selection of new peptide u2 to be merged
    selb = u2.select_atoms(sel2)
    u = mda.core.universe.Merge(sela, selb)
    u_atm = u.select_atoms('all')
    u_atm.residues.resids = np.arange(1,len(u_atm.residues.resids)+1) 

    return u

