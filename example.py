import logging, sys
logging.basicConfig(stream=sys.stdout, level=logging.DEBUG) # can switch to logging.INFO

# Keeps matplotlib working on ssh sessions without a display.
import os
if not 'DISPLAY' in os.environ.keys():
    import matplotlib
    matplotlib.use('Agg')
from matplotlib import pyplot as plt

import Anchor

""" Recap of analysis parameters:
 - ImageStackName: name (path) of .tif file with bright-field image data
 - StormFileName: name (path) of .csv file with STORM fluorescence
    localization data
 - BeadFileNames (optional): string (or array of strings) with path 
    to .csv files containing ThunderSTORM data for each tracked bead 
    (one file per one and only one bead)
 - FrameGrouping: number of images in a single group for averaging (recommended around 20)
 - BeadStabilizing: number of intial frames (after averaging!) for which the stabilization
    is performed (the mean position of these frames is set to 0) (recommended around 10)
 - StackUpscaling: upscaling factor of images before cross-correlation (not recommended
    - use None)
 - StackSmoothing: True/False, decides if gaussian smoothing is used on images before
    determining drift (recommended True)
 - StackCropping: number of edge pixels cropped from mask images before calculating
    crosscorrelation (recommended more than 5 pixels, but no more than 5-10% of image size)
 - StackStabilizing: equivalent of BeadStabilizing for bright-field stack images
 - StackNmPerPixel: property of camera (usually 100)
 - FrameSetup: used for conversion between Stack frame IDs and STORM frame IDs
    - StormFrames: array of at least 2 STORM frame IDs
    - StackFrames: array of same size, with corresponding stack frame IDs
    """ 

AnalysisParameters = \
    {'ImageStackName' : 'stack.tif',\
     'StormFileName' : 'storm.csv',\
     'BeadFileNames' : ['beads/bead1.csv', 'beads/bead2.csv',\
                        'beads/bead3.csv'],\
     'FrameGrouping' : 20,\
     'BeadStabilizing' : 10,\
     'StackUpscaling' : None,\
     'StackSmoothing' : True,\
     'StackCropping' : 32,\
     'StackStabilizing' : 10,\
     'StackNmPerPixel' : 100,\
     'FrameSetup' : {'StormFrames' : [500, 19999],\
                     'StackFrames' : [  1, 10000]},\
     'OutputName' : 'stormdata_corrected.csv',\
     }


print(AnalysisParameters)
A = Anchor.Analyser(AnalysisParameters)
 
A.startAnalysis()

#plt.plot()
plt.savefig('drift.png')

A.tStorm.plotMeanDriftVsXcorr(A.im_stack)
plt.savefig('vs.png')

A.tStorm.plotError(A.im_stack)
#plt.plot()
plt.savefig('error.png')
