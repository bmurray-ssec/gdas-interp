'''
Created on May 4, 2011

@author: nickb
'''

import numpy as np
import pygrib
from File import File


class GribFile(File):
    """
    Common interface for (reading) Grib files.
    """

    def __init__(self, filename):
        File.__init__(self, filename)
        self.inFile = None
        self.pygrib_ver = [int(x) for x in pygrib.__version__.split('.')]


    def readData(self, varName, **kwargs):
        # only reopen the file if it is not already open
        if self.inFile == None:
            self.inFile = pygrib.open(self.filename)

        iop = None # Indicator of Parameter

        # TODO: GENERALIZE THIS IF POSSIBLE
        # Can't seem to get the name out properly, but 154 definitely looks
        # to be the channel for ozone mixing ratio, see:
        # http://www.cpc.ncep.noaa.gov/products/wesley/opn_gribtable.html
        if varName == "O3MR":
            iop = 154

        if iop is not None:
            records = self.inFile.select(indicatorOfParameter=iop, **kwargs)
        else:
            records = self.inFile.select(name=varName, **kwargs)

        return records


    def readDataLatLon(self, varName, lat, lon, **kwargs):
        """Given a measurement name and (lat,lon) coords, return an array 
        containing the measurement values at that (lat,lon) for each 
        pressure level, sorted by pressure level.
        """
        # TODO: interpolate between lat and lon if necessary
        
        retList = []
        
        records = self._readRecords(varName, **kwargs)

        # get the index of the (lat,lon) pair
        lats, lons = records[0].latlons()   
        latIdx, lonIdx = self._getIndex(lats, lons, lat, lon)

        for record in records:
            # We only care about the pressure levels
            if record.typeOfLevel == 'isobaricInhPa':
                # Append this value as a 2D 1x1 grid, NOTE: right now, just a value
                retList.append((record.level, record.values[latIdx][lonIdx]))

        retList.sort()

        presVec = np.ndarray((len(retList)))
        dataVec = np.ndarray((len(retList)))
        for i in range(len(retList)):
            presVec[i] = retList[i][0]
            dataVec[i] = retList[i][1]

        return presVec, dataVec


    def readDataAllLatLon(self, varName, **kwargs):
        """Read a measurement at all (lat,lon) coords in this file. Return one
        grid with the coordinate mappings, and one grid with the profiles which
        correspond to the coordinate mappings.

        in:
        varName - The name of the measurement for which to get the profiles

        out:
        coordGrid - A 2D grid with a coordinate at coordGrid[latIdx][lonIdx] 
                    corresponding to a profile at profGrid[k][latIdx][lonIdx]
        profGrid  - A 3D grid which maps the coordinate grid to the profiles 
                    themselves. profGrid is indexed by:
                    (pressure level, latIdx, lonIdx)
        """

        records = self._readRecords(varName, **kwargs)

        # get all the valid (lat,lon) pairs and store them into coordGrid
        lats, lons = records[0].latlons()
        coordGrid = self._getLatLonGrid(lats, lons)

        # get profiles and sort by pressure level
        tmpGrid = []
        for record in records:
            if record.typeOfLevel == 'isobaricInhPa':
                # store our 3D profile grid as a numpy.ndarray
                tmpGrid.append((record.level, record.values))

        tmpGrid.sort(key=lambda prof: prof[0])

        # copy profiles into our output format
        presVec = np.ndarray((len(tmpGrid)))
        dataGrid = np.ndarray((len(tmpGrid), len(coordGrid), len(coordGrid[0])))
        for i in range(len(tmpGrid)):
            presVec[i] = tmpGrid[i][0]
            dataGrid[i] = tmpGrid[i][1]

        return coordGrid, presVec, dataGrid


    def _readRecords(self, varName, **kwargs):
        # only reopen the file if it is not already open
        if self.inFile == None:
            self.inFile = pygrib.open(self.filename)

        # read the records
        iop = None
        if varName == "O3MR":
            iop = 154

        records = None
        if iop is not None:
            records = self.inFile.select(indicatorOfParameter=iop, **kwargs)
        else:
            records = self.inFile.select(name=varName, **kwargs)
        return records


    def _getLatLonGrid(self, lats, lons):
        """Given lats and lons objects from pygrib, return the grid of (lat,lon) 
        pairs corresponding to the indices of the grid.
        """
        numLats, numLons = lats.shape

        #coordGrid = [[0 for lon in range(numLons)] for lat in range(numLats)]
        coordGrid = np.ndarray((numLats, numLons, 2))
        for lat in range(len(lats)):
            for lon in range(len(lons[lat])):

                # HACK: deal with a bug in pygrib 1.9.6 and 1.9.8 which flips
                #   the latitude indices in the coordinate grid.
                if self.pygrib_ver >= [1, 9, 6]:
                    coordGrid[lat, lon, 0] = lats[len(lats) - 1 - lat][lon]
                else:
                    coordGrid[lat, lon, 0] = lats[lat][lon]

                coordGrid[lat, lon, 1] = lons[lat][lon]


        return coordGrid

            
    def _getIndex(self, lats, lons, lat, lon):
        numLats, numLons = lats.shape

        # add 1 because the min and max are both inclusive
        latRange = lats.max() - lats.min() + 1
        lonRange = lons.max() - lons.min() + 1

        latIncr = latRange / float(numLats)
        lonIncr = lonRange / float(numLons)
            
        # make sure the index is an integer
        latIndex = int(round((lat - lats.min()) / latIncr))
        lonIndex = int(round((lon - lons.min()) / lonIncr))

        return (latIndex, lonIndex)

