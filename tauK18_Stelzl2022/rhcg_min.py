"""
rhcg_min
---------
development version of rhcg - run example for short tau K18
"""
from lib import translate_concept, find_clashes, merge_universe
import sys
import pathlib2
import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import align

def reweighted_hierarchical_chain_growth(path2pep, path2weights, pair ,start = 0 ,
                          range_end = 14, steps2subsequent_fragment = 2 , max_fragments = 100  , index = -1 , 
                          rmsd_cutOff = 0.6 ,clash_D = 2.0, theta=10.0):    
    
    """ Perform RHCG to grow ensembles of disordered proteins from fragments 

    Parameters
    ----------
    path2pep : string
        path to MD fragment library.
    path2weights : string
        path to MD fragment library.
    pair : integer
        fragment pair you are growing 
    start : integer, optional
        ID of fragment you want to start with. The default is 0.
    range_end : integer, optional
        ID of fragment you want to end with. The default is 14.
    steps2subsequent_fragment : integer, optional
        steps between fragments == number of md fragments you merged in the previous step. 
        The default is 2.
    max_fragments : integer, optional
        maximum number of models you want to grow. The default is 100.
    index : integer, optional
        last residue to exlude from clash search. The default is -1.
    rmsd_cutOff : float, optional
        cut-off for the RMSD of the fragment alignment. The default is 0.6
    clash_distance : float, optional
        max. allowed distance between atoms. The default is 2.0
    theta : float, optional
        which theta to use for the bias. The default is 10.0.
    
    Returns
    -------
    None
    
    """

    # Loop: Loop through number of growing fragments in steps of paired fragments to grow the chain
    print(pair)
    
    # counter for rejection / success / assembly attemtps per level
    reject_clash = 0
    reject_RMSD = 0
    attempts = 0
    success = 0
    # max number of newly paired fragments
    k_max = max_fragments
    numPairs = len(range(start, range_end ,pair))
    
    # array to store the product of weights cW1 * cW2 of assembled fragments / pairs
    ar = np.zeros((numPairs , k_max))
    # array  to store the frame indices of successfully assembled fragments / pairs
    conf_index = np.zeros((numPairs, k_max, 2))
    # counter to fetch the product of weights cW1 * cW2 from the previous level
    c1 = 0
    
    for m_i, i in enumerate(range(start, range_end ,pair)):
        print(i)
        # counter to fetch the product of weights cW1 * cW2 from the previous level
        c2 = c1+1
        #print(c1, c2, ar.shape)
        i2 = int(i+steps2subsequent_fragment)
        dire = "pair{}_{}".format(pair,i)
        k = 0
        pathlib2.Path(dire).mkdir(parents=True, exist_ok=True)
        
        # indices or residue numbers for the assembly functions
        #translate_concept
        index1_alnB=-3
        index1_alnE=-2
        index2_alnB=1
        index2_alnE=2
        # find_clashes
        index1_clashB = -3
        index1_clashE = -1
        index2_4clash = 2
        # merge_universe
        res2merge_begin1 = 1 #resID
        res2merge_end1 = -3 #index
        res2merge_begin2 = 3
        res2merge_end2 = -1

        # id pairs assembled  in previous level
        old_pair1 = int(pair/2)
        old_pair2 = int(old_pair1)
        
        if pair == 2:
            # simulated fragments for level 1== pair 2
            # adapt to run with minimal example
            u1= mda.Universe("{}{}/pair0.pdb".format(path2pep , i) ,
                 "{}{}/pair.xtc".format(path2pep , i))

            u2 = mda.Universe("{}{}/pair0.pdb".format(path2pep , i2) ,
                             "{}{}/pair.xtc".format(path2pep , i2))

        # adapt to run with minimal example
        if pair == 4 and i == 12: 
            print('exception')
            u1 = mda.Universe('pair{}_{}/pair0.pdb'.format(old_pair1, i) ,
                    'pair{}_{}/pair.xtc'.format(old_pair1, i))
            # adapt to run with minimal example
            u2 = mda.Universe("{}{}/pair0.pdb".format(path2pep , i2) ,
                             "{}{}/pair.xtc".format(path2pep , i2))     
        if pair == 16:
            res2merge_begin1 = 2 #resID
            res2merge_end2 = -2 #index
        if pair != 2 and not (pair == 4 and i == 12):
            u1 = mda.Universe('pair{}_{}/pair0.pdb'.format(old_pair1, i) ,
                    'pair{}_{}/pair.xtc'.format(old_pair1, i))
            u2 = mda.Universe('pair{}_{}/pair0.pdb'.format(old_pair2, i2) ,
                    'pair{}_{}/pair.xtc'.format(old_pair2, i2))
            
        writePDB = True
        
        while k < k_max:
            attempts += 1
            
            # load weights to draw random frames according to weights
            # fetch the product of unnormalized weights cW1 * cW2 from the previous level
            if pair == 2:
                # adapt to run with minimal example
                w1 = np.genfromtxt("{}{}/wopt/w_theta{}_dat.txt".format(path2weights, i, theta))
                w2 = np.genfromtxt("{}{}/wopt/w_theta{}_dat.txt".format(path2weights, i2, theta))
                r1 = np.random.choice(len(u1.trajectory),p=w1)
                r2 = np.random.choice(len(u2.trajectory),p=w2)

            elif pair == 4 and i == 12:
                c1 = -1 
                w1 = np.load("ar_pair{}.npy".format(old_pair1))[c1]
                # adapt to run with minimal example
                w2 = np.genfromtxt("{}{}/wopt/w_theta{}_dat.txt".format(path2weights, i2, theta))                
                r1 = np.random.randint(len(u1.trajectory))
                r2 = np.random.choice(len(u2.trajectory),p=w2)

            else:
                w1 = np.load("ar_pair{}.npy".format(old_pair1))[c1]
                w2 = np.load("ar_pair{}.npy".format(old_pair2))[c2]                
                r1 = np.random.randint(len(u1.trajectory))
                r2 = np.random.randint(len(u2.trajectory))
            
            # load random frame
            u1.trajectory[r1]
            u2.trajectory[r2]

            ### grow - fragment assembly
            select = translate_concept(u1, u2,index1_alnB ,index1_alnE , index2_alnB , index2_alnE )
            old,new = align.alignto(u1, u2, select=select, weights="mass")
            if new < rmsd_cutOff:
                clashes = find_clashes(u1 , u2 , index1b=index1_clashB , index1e=index1_clashE, index2=index2_4clash , clash_radius = clash_D)
                if clashes < 1:
                    u = merge_universe(u1, u2, sel1="resid {}:{}".format(res2merge_begin1 , u1.atoms.residues[res2merge_end1].resid ),
                                sel2="resid {}:{}".format(res2merge_begin2 , u2.atoms.residues[res2merge_end2].resid ))

                    ar[m_i , k] = w1[r1] * w2[r2]
                    conf_index[m_i, k] = [r1, r2]
                    
                    if writePDB :
                       # save atom positions in a pdb file as topology file
                        u.atoms.write("{}/pair{}.pdb".format(dire,k))
                        at = np.zeros((k_max , u.atoms.n_atoms , 3))
                        writePDB = False

                    at[k] = u.atoms.positions 
                    k+=1
                    success += 1
                else:
                    reject_clash += 1
                    continue
            else:
                reject_RMSD += 1
        c1 +=2
        allPairs = mda.Universe("{}/pair0.pdb".format(dire) , at)
        allPairs.atoms.write("{}/pair.xtc".format(dire) , frames='all') 
    return None

rcAll = 0
rrAll = 0
start = 0

# fulllength tau K18
range_end =  14
path = sys.argv[1]
path2pep = sys.argv[2]
path2weights = sys.argv[3]
kmax = 100
theta = 10.0

rmsd_cutOff = 0.6
clash_D = 2.0 
level = 1

for i in [2, 4, 8, 16]:
    start = 0
    subsequent_fragment = int(i/2)
    reweighted_hierarchical_chain_growth(path2pep = path2pep, path2weights = path2weights, pair = i,
                                         steps2subsequent_fragment = subsequent_fragment,
                                         max_fragments = kmax , start=start,  range_end = range_end, 
                                         rmsd_cutOff = rmsd_cutOff , clash_D = clash_D, theta=theta)

    level += 1

