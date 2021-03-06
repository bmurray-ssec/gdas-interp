#!/home/bmurray/svn/ShellB3/trunk/ShellB3/bin/python
'''
Script used to test our interpolation code

@author: bmurray
'''

import os.path
import time, sys
import gdas_interp as g
import matplotlib.pyplot as plt
from GribFile import GribFile

IN_FILE = 'data/gdas1.PGrbF00.060828.18z'
VAR_NAME = 'Temperature'
#VAR_NAME = 'Relative humidity'
EPSILON = 0.00000001

def main():
    global IN_FILE
    levIdx = 90

    if len(sys.argv) >= 2:
        IN_FILE = sys.argv[1]
    if len(sys.argv) >= 3:
        levIdx = int(sys.argv[2])

    print 'For input file: %s' % IN_FILE
    print 'Measuring: %s' % VAR_NAME
    
    #grb = GribFile('data/gdas1.PGrbF00.060828.18z')

    #pressure = [x + 10 for x in range(20)]

    print 'getting grid data'
    #interp = g.vert_interp_grid(VAR_NAME, filename='data/gdas1.PGrbF00.130901.06z')
    coordGrid, presLevels, profGrid = g.vert_interp_grid(VAR_NAME, filename=IN_FILE)

    plt.imshow(profGrid[levIdx])
    plt.title('"%s"\n %s at pressure %f' % \
            (os.path.basename(IN_FILE), VAR_NAME, presLevels[levIdx]))
    plt.show()

    '''
    nLats = len(coordGrid)
    nLons = len(coordGrid[0])

    for i in range(nLats * nLons):
        lat, lon = divmod(i, nLons)

        col = g.vert_interp(VAR_NAME, coordGrid[lat, lon, 0], coordGrid[lat, lon, 1], filename=IN_FILE)
        gridCol = [ profGrid[j][lat][lon] for j in range(len(profGrid)) ]

        for j in xrange(len(col)):
            col[j] = col[j][1]

        failed = False
        # make sure profiles are of equal lengths
        if len(col) != len(gridCol):
            failed = True
            print 'Single column length does not match grid prof length: [%d -> %d]' % \
                (len(col), len(gridCol[i]))
            exit(1)

        # make sure each profile column at a given (lat,lon) matches the 
        # corresponding grid column at that (lat,lon)
        for j in range(len(col)):
            if abs(col[j] - gridCol[j]) > EPSILON:
                failed = True
                print 'Profiles at (%f, %f) do not match at index %3d [col: %9.6f -> grid: %9.6f]' % \
                    (lat, lon, j, col[j], gridCol[j])
                exit(1)

        if not failed:
            print '%d succeeded' % i        

    # Right now, profGrid is just a list of profiles (nLats*nLons, 101)
    #coordGrid, profGrid = grb.readDataAllLatLon('VAR_NAME')

    #incProfGrid = [[ [] for lat in range(len(coordGrid))] for lon in range(len(coordGrid[0]))]

    #print len(incProfGrid), len(coordGrid)
    #print len(incProfGrid[0]), len(coordGrid[0])
    '''


if __name__ == '__main__':
    main()

