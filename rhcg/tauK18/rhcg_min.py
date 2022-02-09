from lib import mk_universe, translate_concept, findClashes, merge_universe
import sys
import pathlib2
import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import align

def merge_fragment_hierar(path2pep, path2weights, pair ,start = 0 ,
                          range_end = 45, steps2subsequent_fragment = 2 , max_fragments = 100  , index = -1 , 
                          rmsd_cutOff = 0.6 ,clash_D = 2.0 ,
                           theta=10.0, biased=True):
    
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
        ID of fragment you want to end with. The default is 45.
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
    biased : boolean, optional
        run rhcg (use weights for fragment conformation choice) or just hcg if == False. 
        The default is True.

    Returns
    -------
    None.

    """

    # Loop: Loop through number of growing fragments in steps of paired fragments to grow the chain
    print(pair)
    m=[]
    reject_clash = 0
    reject_RMSD = 0
    attempts = 0
    success =0
    # max number of newly paired fragments
    k_max = max_fragments
    numPairs = len(range(start, range_end ,pair))
    ar = np.zeros((numPairs , k_max))
    conf_index = np.zeros((2, numPairs , k_max))
    c1 = 0
    #print numPairs
    for ni, i in enumerate(range(start, range_end ,pair)):
        print(i)
#        steps2subsequent_fragment = pair/2
        build_ar = True
        # index for subsequent fragment to add to fragment i
        c2 = c1+1
        #print(c1, c2, ar.shape)
        i2 = int(i+steps2subsequent_fragment)
        dire = "pair{}_{}".format(pair,i)
        k = 0
        pathlib2.Path(dire).mkdir(parents=True, exist_ok=True)
        # list of indices with successful alignment
        matches = []
        #m = []
       
        index1_alnB=-3
        index1_alnE=-2
        index2_alnB=1
        index2_alnE=2

        index1_clashB = -3
        index1_clashE = -1
        index2_4clash = 2

        res2merge_begin1 = 1 #resID
        res2merge_end1 = -3 #index
        res2merge_begin2 = 3
        res2merge_end2 = -1


        old_pair1 = int(pair/2)
        old_pair2 = int(old_pair1)
        if pair == 2:
            # simulated fragments (pentapep) to start the growing           
            # adapt to run with minimal example
            u1= mk_universe("{}{}/pair0.pdb".format(path2pep , i) ,
                 "{}{}/pair.xtc".format(path2pep , i))

            u2 = mk_universe("{}{}/pair0.pdb".format(path2pep , i2) ,
                             "{}{}/pair.xtc".format(path2pep , i2))

        # adapt to run with minimal example
        if pair == 4 and i == 12: 
            print('exception')
            u1 = mda.Universe('pair{}_{}/pair0.pdb'.format(old_pair1, i) ,
                    'pair{}_{}/pair.xtc'.format(old_pair1, i))
            # adapt to run with minimal example
            u2 = mk_universe("{}{}/pair0.pdb".format(path2pep , i2) ,
                             "{}{}/pair.xtc".format(path2pep , i2))

#            index2_alnB=2
#            index2_alnE=3
#         
#            res2merge_begin2 = 4
#            res2merge_end2 = -1
#            index2_4clash = 3

      
        if pair == 16:
            res2merge_begin1 = 2 #resID
            res2merge_end2 = -2 #index
        if pair != 2 and not (pair == 4 and i == 12):
            u1 = mda.Universe('pair{}_{}/pair0.pdb'.format(old_pair1, i) ,
                    'pair{}_{}/pair.xtc'.format(old_pair1, i))
            u2 = mda.Universe('pair{}_{}/pair0.pdb'.format(old_pair2, i2) ,
                    'pair{}_{}/pair.xtc'.format(old_pair2, i2))
         #   print(i, i2)
        while k < k_max:
            attempts += 1
            if k == 0: # or pair == 42 :
                writePDB = True
            else:
                writePDB = False
            if biased:
                if pair == 2:
                    # adapt to run with minimal example
                    w1 = np.genfromtxt("{}{}/wopt/w_theta{}_dat.txt".format(path2weights, i, theta))
                    w2 = np.genfromtxt("{}{}/wopt/w_theta{}_dat.txt".format(path2weights, i2, theta))

                    r1 = np.random.choice(len(u1.trajectory),p=w1)
                    r2 = np.random.choice(len(u2.trajectory),p=w2)
    
                elif pair == 4 and i == 12:
                    c1= -1 # 6 
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
            else:
                w1 = np.ones(len(u1.trajectory))/ len(u1.trajectory)
                w2 = np.ones(len(u2.trajectory))/ len(u2.trajectory)

                r1 = np.random.choice(len(u1.trajectory),p=w1)
                r2 = np.random.choice(len(u2.trajectory),p=w2)
            u1.trajectory[r1]
            u2.trajectory[r2]


            ### grow
            select = translate_concept(u1, u2,index1_alnB ,index1_alnE , index2_alnB , index2_alnE )
            old,new = align.alignto(u1, u2, select=select, weights="mass")
            if new < rmsd_cutOff:
                clashes = findClashes(u1 , u2 , index1b=index1_clashB , index1e=index1_clashE, index2=index2_4clash , clash_radius = clash_D)
                if clashes < 1:
                    u = merge_universe(u1, u2, sel1="resid {}:{}".format(res2merge_begin1 , u1.atoms.residues[res2merge_end1].resid ),
                                sel2="resid {}:{}".format(res2merge_begin2 , u2.atoms.residues[res2merge_end2].resid ))

                    m.append([pair,i,k,r1,r2])
                    ar[ni , k] = w1[r1] * w2[r2]
                    conf_index[0, ni , k] = r1
                    conf_index[1, ni , k] = r2
                    if writePDB :
                       # save atom positions in a pdb file
                        u.atoms.write("{}/pair{}.pdb".format(dire,k))
                    if build_ar:
                        at = np.zeros((k_max , u.atoms.n_atoms , 3))
                        build_ar = False

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
        allPairs.atoms.write("{}/pair.xtc".format(dire) , frames='all') ## new MDA version
    return m, reject_clash, reject_RMSD , attempts , success , ar , conf_index



rcAll = 0
rrAll = 0
start = 0



# fulllength tau K18
range_end =  14
path2pep = sys.argv[1]
path2weights = sys.argv[2]
max_fragments = int(sys.argv[3])
theta=sys.argv[4]

rmsd_cutOff = 0.6
clash_D = 2.0 
level = 0

for i in [2, 4, 8, 16]:
    start = 0
    subsequent_fragment = int(i/2)
    matches, rc ,rr, att, succ, ar, confi = merge_fragment_hierar(path2pep = path2pep, path2weights = path2weights, pair = i , steps2subsequent_fragment = subsequent_fragment ,  
                                               max_fragments = max_fragments , start=start,  range_end = range_end, 
                                               rmsd_cutOff = rmsd_cutOff , clash_D = clash_D, theta=theta)
    rcAll+=rc
    rrAll+=rr
    np.save('ar_pair{}.npy'.format(i) , ar)
    np.save('confIndex_pair{}.npy'.format(i) , confi)

    out = open( "statistics{}.txt".format(level), "w")
    out.write("\nstatistics on level {}\n"
              "stitching attempts:\t\t\t{}\n"
              "rejection because of RMSD:\t\t{}\n"
              "rejection because of clashes:\t\t{}\n"
              "success:\t\t\t\t{}\n"
              "acceptance rate:\t\t\t{}\n".format(level, att, rr, rc, succ, round(float(succ)/float(att),3)))

    if i == 16:
        out.write("rejection because of clashes in total for {} FL models:\t\t{}\n"
                        "rejection because of RMSD in total for {} FL models:\t\t{}".format(max_fragments,rcAll,max_fragments,rrAll))


    out.close()
    level += 1

