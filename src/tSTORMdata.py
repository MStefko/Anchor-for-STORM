from beadData import *
from ImageStack import *

import os
import copy
from operator import add
import logging, sys
import math
import numpy as np
from scipy import misc
from scipy.ndimage.filters import gaussian_filter
from scipy.fftpack import fft2, ifft2
from scipy.interpolate import RectBivariateSpline, griddata, interp2d
from scipy.optimize import minimize
if not 'DISPLAY' in os.environ.keys():
    import matplotlib
    matplotlib.use('Agg')
from matplotlib import pyplot as plt
import matplotlib.cm as cm



class tSTORMdata:
    """ Container class for holding data from ThunderSTORM bead tracking algorithm. Multiple beads
    can (and should) be tracked.
    Input data must satisfy following requirements:
     - Each file contains data for only one bead.
     - The bead is tracked in each frame (a sparse set of frames missing should be OK, but
        can cause unexpected behavior). Two beads in one frame will cause errors.
    Parameters:
     - filenames (list of strings):
        Paths to .csv files containing data about beads. One and only one bead per file.
     - grouping (int):
        After loading the data, frames are grouped together and averaged, so a smoother and
        smaller dataset is then acquired. Helps remove noise.
        """
    def __init__(self, filenames = None, grouping = None, stabilizing = None):
        self.beads = []
        self.grouping = None # Set to a value if grouping is used
        self.stabilized = None # Set to a value if stabilizing is used.
        self.meanDrift = None
        self.driftVariance = None
        if filenames:
            self.loadData(filenames)
        if grouping:
            self.averageData(grouping)
            logging.info("Averaging bead data upon loading into groups of "+str(grouping))
        if stabilizing:
            self.stabilizeData(frames = stabilizing)
            logging.info("Stabilizing bead data by first "+str(stabilizing)+" frames.")
        return

    def loadData(self, filenames):
        if isinstance(filenames, str):
            filenames = [filenames]
        for filename in filenames:
            self.beads.append( beadData(filename) )
        self.meanDrift = None
        self.driftVariance = None
        return

    def clearData(self):
        self.beads = []
        self.meanDrift = None
        self.driftVariance = None
        return

    def averageData(self, grouping = 20):
        self.grouping = grouping
        for bead in self.beads:
            bead.averageData(grouping)
        self.meanDrift = None
        self.driftVariance = None
        return

    def stabilizeData(self, frames = 10):
        self.stabilized = frames
        for bead in self.beads:
            bead.stabilizeData(frames = frames)
        self.meanDrift = None
        self.driftVariance = None
        return

    def getBeadDrift(self, no):
        return self.beads[no].drift

    def getMeanDrift(self):
        """ Computes the mean drift of all loaded beads, frame-by-frame. """
        if self.meanDrift:
            return self.meanDrift
        if not self.grouping:
            logging.error("Error: Drifts are not grouped - unstable behavior, cancelling!")
            return None
        if not self.beads:
            logging.error("Error: No beads loaded.")
            return None
        if not all(len(i.drift) == len(self.beads[0].drift) for i in self.beads):
            logging.error("Error: Drift data not of same length. (Sparsely loaded or overloaded bead?)")
            return None
        meanDrift = []
        # Cycle over each frame
        for i in range(len(self.beads[0].drift)):
            x=0; y=0; n=None
            # Cycle over each bead
            for bead in self.beads:
                x += bead.drift[i][1]
                y += bead.drift[i][2]
                n = bead.drift[i][0]
            # Divide data by number of beads
            meanDrift.append( (n, x/len(self.beads), y/len(self.beads)) )
        self.meanDrift = meanDrift
        return self.meanDrift

    def getDriftVariance(self):
        """ Computes the variance of drift data between different beads
        frame-by-frame.
        """
        if self.driftVariance:
            return self.driftVariance
        if not self.getMeanDrift():
            return None
        if len(self.beads) == 1:
            logging.error("Error: 1 bead loaded, variance unusable.")
            return None
        driftVariance = []
        for i in range(len(self.beads[0].drift)):
            Vx=0; Vy=0; n=None
            for bead in self.beads:
                Vx += (bead.drift[i][1] - self.meanDrift[i][1])**2
                Vy += (bead.drift[i][2] - self.meanDrift[i][2])**2
                n = bead.drift[i][0]
            driftVariance.append( (n, Vx/len(self.beads), Vy/len(self.beads)) )
        self.driftVariance = driftVariance
        return self.driftVariance

    def plotDrifts(self, skip = 1, fig = 1):
        if fig:
            plt.figure(fig)
        plt.subplot(211)
        for bead in self.beads:
            n = [a[0] for a in bead.drift]
            x = [a[1] for a in bead.drift]
            plt.plot(n[::skip], x[::skip], '.')
        plt.ylabel('X Drift [nm]')

        plt.subplot(212)
        for bead in self.beads:
            n = [a[0] for a in bead.drift]
            y = [a[2] for a in bead.drift]
            plt.plot(n[::skip], y[::skip], '.')
        plt.xlabel('Frame')
        plt.ylabel('Y Drift [nm]')
        #plt.show()
        
    def plotMeanDrift(self, skip = 1, fig = 2):
        if fig:
            plt.figure(fig)
        d = self.getMeanDrift()
        plt.plot(d[0][::skip], d[1][::skip],'r.',label='X')
        plt.plot(d[0][::skip],d[2][::skip],'g.',label='Y')
        plt.xlabel('Frame')
        plt.ylabel('Drift [nm]')
        plt.legend()
        #plt.show()

    def plotVariance(self, skip = 1, fig = 3):
        if fig:
            plt.figure(fig)
        d = self.getDriftVariance()

        n = [a[0] for a in d]
        x = [a[1] for a in d]
        y = [a[2] for a in d]
        plt.plot(n[::skip], x[::skip],'r.',label='X')
        plt.plot(n[::skip],y[::skip],'g.',label='Y')
        plt.xlabel('Frame')
        plt.ylabel('Drift Variance [nm]')
        plt.legend()
        #plt.show()

    def plotMeanDriftVsXcorr(self, stack, skip = 1, fig = 4):
        if fig:
            plt.figure(fig)
        d = self.getMeanDrift()
        xcorr = stack.getDrift()
        if (len(d) != len(xcorr)):
            logging.error("Error: Lengths do not match up! Bead: " + str(len(d)) +\
                ", Xcorr: " + str(len(xcorr)))
            return

        n = [a[0] for a in d]
        xB = [a[1] for a in d]
        xX = [a[1] for a in xcorr]

        plt.subplot(211)
        plt.plot(n[::skip], xB[::skip], 'b.',label='tSTORM')
        plt.plot(n[::skip], xX[::skip], 'r.',label = 'Xcorr')
        plt.ylabel('X Drift [nm]')
        plt.legend()

        yB = [a[2] for a in d]
        yX = [a[2] for a in xcorr]

        plt.subplot(212)
        plt.plot(n[::skip], yB[::skip], 'b.')
        plt.plot(n[::skip], yX[::skip], 'r.')
        plt.xlabel('Frame')
        plt.ylabel('Y Drift [nm]')
        plt.legend()
        #plt.show()

    def plotError(self, stack, skip = 1, fig = 5):
        """ Plots the difference between bead and xcorrelation drift frame-by-frame.
        """
        if fig:
            plt.figure(fig)
        d = self.getMeanDrift()
        xcorr = stack.drift
        if (len(d) != len(xcorr)):
            logging.error("Error: Lengths do not match up! Bead: " + str(len(d)) +\
                ", Xcorr: " + str(len(xcorr)))
            return

        n = [a[0] for a in d]
        xB = np.array([a[1] for a in d])
        xX = np.array([a[1] for a in xcorr])
        x = xB - xX

        yB = np.array([a[2] for a in d])
        yX = np.array([a[2] for a in xcorr])
        y = yB - yX
        
        logging.info("Mean error [nm]: \n  x: " + str(np.mean(x)) +\
            "\n  y: " + str(np.mean(y)) + "\nStandard deviation [nm]: \n  x: " +\
            str(np.std(x)) + "\n  y: " + str(np.std(y)))
        plt.plot(n[::skip], x[::skip], 'b.',label='X')
        plt.plot(n[::skip], y[::skip], 'r.',label='Y')
        plt.xlabel('Frame')
        plt.ylabel('Drift Error [nm]')
        plt.legend()
        #plt.show()