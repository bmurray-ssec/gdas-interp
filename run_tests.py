#!/home/bmurray/svn/ShellB3/trunk/ShellB3/bin/python
'''
Script used to test our interpolation code

@author: bmurray
'''

from GribFile import GribFile
import time
import gdas_interp as g

IN_FILE = 'data/gdas1.PGrbF00.060828.18z'
VAR_NAME = 'Temperature'
#VAR_NAME = 'Relative humidity'
EPSILON = 0.00000001

def main():
    print 'For input file: %s' % IN_FILE
    print 'Measuring: %s' % VAR_NAME
    
    #grb = GribFile('data/gdas1.PGrbF00.060828.18z')

    #pressure = [x + 10 for x in range(20)]

    print 'getting grid data'
    #interp = g.vert_interp_grid(VAR_NAME, filename='data/gdas1.PGrbF00.130901.06z')
    coordGrid, presCol, profGrid = g.vert_interp_grid(VAR_NAME, filename=IN_FILE)

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
                    (coordGrid[lat,lon,0], coordGrid[lat,lon,1], j, col[j], gridCol[j])
                exit(1)

        if not failed:
            print 'lat: %f, lon: %f succeeded' % (coordGrid[lat,lon,0], coordGrid[lat,lon,1])

    # Right now, profGrid is just a list of profiles (nLats*nLons, 101)
    #coordGrid, profGrid = grb.readDataAllLatLon('VAR_NAME')

    #incProfGrid = [[ [] for lat in range(len(coordGrid))] for lon in range(len(coordGrid[0]))]

    #print len(incProfGrid), len(coordGrid)
    #print len(incProfGrid[0]), len(coordGrid[0])


if __name__ == '__main__':
    main()

