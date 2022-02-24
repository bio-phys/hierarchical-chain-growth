from __future__ import division
import numpy as np
import argparse

def num_ions(Nw, salt_c):
    return Nw*(salt_c/55.5)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('nw', help='The number of water molecules', type=np.float64)
    parser.add_argument('salt_c', help='The desired salt concentration', type=np.float64)
    args = parser.parse_args()
    
    
    o =  num_ions(args.nw, args.salt_c)
    #print o
    print("{0:.0f}, {1:.0f}".format(o, o) )
    

if __name__ == '__main__':
   main()

