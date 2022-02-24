#!/usr/bin/env python3
"""
tlist_prod
----------
Script to create folder with mdp file for  REMD run with the temperature adapted to a temperature as determined for the replicas
Author Lukas S. Stelzl
"""

import os
import errno
import glob
import shutil

def t_list(start_t, step, num_rep):
    t_list = []
    for i in range(num_rep):
        t  = i*3.5 + nump_rep
        t_list.append(t)
    return t_list

def make_tdirs(t_list):
    '''
    http://stackoverflow.com/questions/600268/mkdir-p-functionality-in-python
    '''
    for t in t_list:
        tf= "%.2f" % t
        path=str(tf)
        try:
            os.makedirs(path)
        except OSError as exc: # Python >2.5
               if exc.errno == errno.EEXIST and os.path.isdir(path):
                  pass
               else: raise
      #  os.makedirs(

def change_mdp(default_mdp, t, fn="prod_", replace_t="300"):
    '''
    http://stackoverflow.com/questions/16622754/how-do-you-replace-a-line-of-text-in-a-text-file-python?lq=1
    '''
    output_fn= os.path.join(t,fn+str(t)+".mdp" )
    with open( default_mdp,"r") as mdp:
         with open( output_fn, "w") as o:
              for line in mdp:
                  if line.startswith("ref_t") or line.startswith("ref-t"):
                     l=line.replace(replace_t, str(t))
                     print(l)
                     o.write(l)
                  elif line.startswith("gen-temp")  or line.startswith("gen_temp"):
                       l=line.replace(replace_t, str(t))
                       print(l)
                       o.write(l)
                  else:
                       o.write(line)

def loop_mdp(t_list, default_mdp, fn="prod_", replace_t="300"):
    for t in t_list:
        tf= "%.2f" % t
        change_mdp(default_mdp, tf, fn=fn, replace_t=replace_t)

def loop_start_struc(t_sched_l,struc_dir, rep=64.0):
    pdb_list = glob.glob(struc_dir)
  #  print t_sched_l
    print(pdb_list)
    for i,t in enumerate(t_sched_l):
        tr = "%.2f" % t
        #if i < rep:
        print(i,tr)
        start_pdb_pth = pdb_list[1+i]
        pth, start_pdb = os.path.split(start_pdb_pth)
        shutil.copy2(start_pdb_pth, os.path.join(tr, tr+"_"+start_pdb))
       # else:
         #    "print dealt with the pdbs for all {:}  replicas".format(rep)

def single_equil_stru(struc, t_sched_l, rep=64.0):
    for i,t in enumerate(t_sched_l):
        tr = "%.2f" % t
        #if i < rep:
        print( i,tr)
        pass

#def main()
#    t_list


#if __name__ == '__main__':
#    main()
