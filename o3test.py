#!/home/bmurray/svn/ShellB3/trunk/ShellB3/bin/python
'''
Script used to test our interpolation code

@author: bmurray
'''

from GribFile import GribFile
import gdas_interp as g
import numpy as np

IN_FILE = 'data/gdas1.PGrbF00.060828.18z'
VAR_NAME = 'O3MR'
LAT = 19.949
LON = -62.717

def main():
    print 'For input file: %s' % IN_FILE
    print 'Measuring: %s' % VAR_NAME

    prof = g.vert_interp_ozone(IN_FILE, LAT, LON, 8)
    
    for p in prof:
        print '%9.4f : %.4f' % (p[0], p[1])


if __name__ == '__main__':
    main()

