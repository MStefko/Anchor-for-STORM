import os
import copy
import logging, sys
import math
import numpy as np


class beadData:
    """ Container class for the trajectory of a single
    bead tracked by ThunderSTORM.

    Uses custom parsing functions, not pandas. Should probably
    be modified to use pandas, since it is more reliable.
    However, it works for ThunderSTORM files quite well.
    """
    def __init__(self):
        self.filename = None
        self.data = None
        self.drift = None
        self.averaged = False
        self.stabilized = False
        return

    def __init__(self, filename):
        self.filename = filename
        self.averaged = False
        self.stabilized = False
        self._loadData(filename)
        return

    def _loadData(self, filename):
    # Imports data and transforms it to a more readable format.
        self.data = self._loadBeadData(filename)
        self.averaged = False
        self.drift = self._extractXYdata( self.data )
        self.stabilized = False
        return

    def _loadBeadData(self, filename):
    # Import data from a .csv file
        with open(filename, 'r') as f:
            first_line = f.readline().strip()
        #parse variable names
        varNames = self._getVarNames(first_line)
        
        #load raw data from file
        rawdata = np.genfromtxt(filename, delimiter=',')
        
        data = {}
        for i in range(len(varNames)):
            if (varNames[i] == 'id' or varNames[i] == 'frame'):
                data[varNames[i]] = np.array([int(d[i]) for d in rawdata[1:]])
            else:
                data[varNames[i]] = np.array([d[i] for d in rawdata[1:]])
        return data
        
    def _extractXYdata(self, data):
    # Transforms data from a dictionary structure to a structure that
    # can be directly compared with the Xcorrelation data.
        oData = []
        for i in range(len(data['frame'])):
            oData.append( (data['frame'][i], data['x'][i] - data['x'][0], data['y'][i] - data['y'][0]) )
        return oData


    def _getVarNames(self, string):
        # parses variable names from csv file
        varNames = []
        varName = ''
        i = 0
        while i<len(string):
            # scan string for first "
            if (string[i]=='"'):
                i+=1
                # get letters until non-alpha character
                while (string[i].isalpha()):
                    varName += string[i]
                    i+=1
                # add variable name to list
                if (varName):
                    varNames += [varName]
                varName = ''
                # scan for another "
                while (string[i+1] != '"'):
                    i+=1
                # check for end of string and repeat
                if (i+2 == len(string)):
                    break
            i+=1
        return varNames

    def averageData(self, grouping):
        if self.averaged:
            logging.warning("Warning: Bead data already averaged, skipping.")
            return
        i = 0; groupstart = 1
        x = 0; y = 0; n = 0;
        averagedData = []
        while (i<len(self.drift)):
            if (self.drift[i][0]-groupstart < grouping):
                x += self.drift[i][1]; y += self.drift[i][2]; n+=1
                i+=1
            else:
                if n:
                    averagedData.append( (math.floor(groupstart+(grouping/2)), x/n, y/n) )
                x = 0; y = 0; n = 0;
                groupstart += grouping
        if n:
            averagedData.append( (math.floor(groupstart+(grouping/2)), x/n, y/n) )

        self.averaged = True
        self.drift = averagedData
        return

    def stabilizeData(self, frames = 10):
        """ Takes several frames at the beginning, computes the mean value of drift
        for this set, and then translates the trajectory so that this value is 0.
        This is useful, because the first frame, which is usually taken as the 
        0.0 value, can be off the actual center of the bead due to noise.
        By taking more frames (but on a timescale where we don't anticipate
        much drift), we can mitigate the influence of noise.
        """
        if not self.drift:
            logging.error("Error in stabilizing: Load data properly first.")
            return
        if not self.averaged:
            logging.error("Error in stabilizing: Average the data first.")
            return
        if len(self.drift)<frames:
            logging.error("Error: Trajectory shorter than stabilization frame.")
            return
        if frames>10:
            logging.warning("Warning: Long stabilization timeframe. ("+str(frames)+\
                "frames) Are you sure you know what you are doing?")
        # Determine the stabilized zero-values
        start = self.drift[0:frames]
        xs = np.array([a[1] for a in start])
        ys = np.array([a[2] for a in start])
        X = np.mean(xs); Y = np.mean(ys)

        # Subtract them from the whole trajectory
        traj = []
        for a in self.drift:
            traj += [(a[0], a[1]-X, a[2]-Y)]
        self.drift = traj
        self.stabilized = True
        return
