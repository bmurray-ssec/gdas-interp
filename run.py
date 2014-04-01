#!/home/bmurray/svn/ShellB3/trunk/ShellB3/bin/python
'''
Script used to test our interpolation code

@author: bmurray
'''

from GribFile import GribFile
import gdas_interp as g
import numpy as np
import matplotlib.pyplot as plt
import pylab
import sys

IN_FILE = 'data/gdas1.PGrbF00.060828.18z'
TEST_FILE = 'test_profile.in'
VAR_NAME = 'Temperature'
#VAR_NAME = 'Relative humidity'

_pressure_levels = \
                   [0.0050,    0.0161,    0.0384,    0.0769,    0.1370, \
         0.2244,    0.3454,    0.5064,    0.7140,    0.9753,    1.2972, \
         1.6872,    2.1526,    2.7009,    3.3398,    4.0770,    4.9204, \
         5.8776,    6.9567,    8.1655,    9.5119,   11.0038,   12.6492, \
        14.4559,   16.4318,   18.5847,   20.9224,   23.4526,   26.1829, \
        29.1210,   32.2744,   35.6505,   39.2566,   43.1001,   47.1882, \
        51.5278,   56.1260,   60.9895,   66.1253,   71.5398,   77.2396, \
        83.2310,   89.5204,   96.1138,  103.0172,  110.2366,  117.7775, \
       125.6456,  133.8462,  142.3848,  151.2664,  160.4959,  170.0784, \
       180.0183,  190.3203,  200.9887,  212.0277,  223.4415,  235.2338, \
       247.4085,  259.9691,  272.9191,  286.2617,  300.0000,  314.1369, \
       328.6753,  343.6176,  358.9665,  374.7241,  390.8926,  407.4738, \
       424.4698,  441.8819,  459.7118,  477.9607,  496.6298,  515.7200, \
       535.2322,  555.1669,  575.5248,  596.3062,  617.5112,  639.1398, \
       661.1920,  683.6673,  706.5654,  729.8857,  753.6275,  777.7897, \
       802.3714,  827.3713,  852.7880,  878.6201,  904.8659,  931.5236, \
       958.5911,  986.0666, 1013.9476, 1042.2319, 1070.9170, 1100.0000]

def get_temp_vector():
    retList = []
    fd = open(TEST_FILE, 'r')
    lines = fd.readlines()
    for line in lines:
        line = line.strip()
        vals = line.split(' ')
        retList.append((float(vals[0]), float(vals[1])))

    retList.sort() 

    presVec = []
    dataVec = []
    for r in retList:
        presVec.append(r[0])
        dataVec.append(r[1])

    return (presVec, dataVec)

def get_rel_hum_vector():
    retList = []
    fd = open(TEST_FILE, 'r')
    lines = fd.readlines()
    for line in lines:
        line = line.strip()
        vals = line.split(' ')
        if float(vals[2]) >= 0:
            retList.append((float(vals[0]), float(vals[2])))

    retList.sort() 

    presVec = []
    dataVec = []
    for r in retList:
        presVec.append(r[0])
        dataVec.append(r[1])

    return (presVec, dataVec)

def main():
    global IN_FILE

    inFile2 = None
    pressure = None
    if len(sys.argv) == 2:
        IN_FILE = sys.argv[1]    
    elif len(sys.argv) == 3:
        pressure = float(sys.argv[2])


    tempProf = get_temp_vector()
    rhProf = get_rel_hum_vector()

    #lat = 19.949
    #lon = -62.717
    lat = -89.0
    lon = 10.0

    #output = g.vert_interp(VAR_NAME, lat, lon, pressure, tempProf=tempProf, rhProf=rhProf)

    output = g.vert_interp(VAR_NAME, lat, lon, pressure, IN_FILE)
    x = []
    y = []

    for o in output:
        print '%10.4f : %7.4f' % (o[0], o[1])
        x.append(o[0])
        y.append(o[1])

    if len(x) > 1:
        plt.plot(x, y, 'b')     
        plt.show()


if __name__ == '__main__':
    main()

