#!/usr/bin/env python3
"""
make_replicas.py
----------------
Script to execute tlist.py: Create folder with mdp file for  REMD run with the temperature adapted to a temperature as determined for the replicas
Inputs: 
- Query temperature list - e.g., generated using a temperature generator for REMD
- Desired .mdp file for REMD production run

Original version provided by Lukas S. Stelzl
Slightly modified by Lisa M. Pietrek
"""
from tlist_prod import *
import numpy as np
import sys


ts_l = np.genfromtxt('/usr/path/script/tmp_dir_list.txt')

make_tdirs(ts_l)

default_mdp ="/usr/path/run_para/md_re.mdp"

loop_mdp(ts_l, default_mdp)


