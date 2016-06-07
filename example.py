import logging, sys
logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)

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
 - StackUpscaling: upscaling factor of images before cross-correlation (not recommended)
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

A1 = \
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
     'figname' : 'orig.png'}

A2 = \
    {'ImageStackName' : 'stack.tif',\
     'StormFileName' : 'storm.csv',\
     'BeadFileNames' : ['beads/bead1.csv', 'beads/bead2.csv',\
                        'beads/bead3.csv'],\
     'FrameGrouping' : 20,\
     'BeadStabilizing' : None,\
     'StackUpscaling' : None,\
     'StackSmoothing' : True,\
     'StackCropping' : 32,\
     'StackStabilizing' : None,\
     'StackNmPerPixel' : 100,\
     'FrameSetup' : {'StormFrames' : [500, 19999],\
                     'StackFrames' : [  1, 10000]},\
     'OutputName' : 'stormdata_corrected.csv',\
     'figname' : 'nostabilizing.png'}

"""A3 = \
    {'ImageStackName' : 'stack.tif',\
     'StormFileName' : 'storm.csv',\
     'BeadFileNames' : ['beads/bead1.csv', 'beads/bead2.csv',\
                        'beads/bead3.csv'],\
     'FrameGrouping' : None,\
     'BeadStabilizing' : 10,\
     'StackUpscaling' : None,\
     'StackSmoothing' : True,\
     'StackCropping' : 32,\
     'StackStabilizing' : 10,\
     'StackNmPerPixel' : 100,\
     'FrameSetup' : {'StormFrames' : [500, 19999],\
                     'StackFrames' : [  1, 10000]},\
     'OutputName' : 'aaa.csv',\
     'figname':'nogrouping.png'}
"""
A4 = \
    {'ImageStackName' : 'stack.tif',\
     'StormFileName' : 'storm.csv',\
     'BeadFileNames' : ['beads/bead1.csv', 'beads/bead2.csv',\
                        'beads/bead3.csv'],\
     'FrameGrouping' : 20,\
     'BeadStabilizing' : 10,\
     'StackUpscaling' : None,\
     'StackSmoothing' : False,\
     'StackCropping' : 32,\
     'StackStabilizing' : 10,\
     'StackNmPerPixel' : 100,\
     'FrameSetup' : {'StormFrames' : [500, 19999],\
                     'StackFrames' : [  1, 10000]},\
     'OutputName' : 'stormdata_corrected.csv',\
     'figname':'nosmoothing.png'}

A5 = \
    {'ImageStackName' : 'stack.tif',\
     'StormFileName' : 'storm.csv',\
     'BeadFileNames' : ['beads/bead1.csv', 'beads/bead2.csv',\
                        'beads/bead3.csv'],\
     'FrameGrouping' : 20,\
     'BeadStabilizing' : 10,\
     'StackUpscaling' : 3,\
     'StackSmoothing' : True,\
     'StackCropping' : 32,\
     'StackStabilizing' : 10,\
     'StackNmPerPixel' : 100,\
     'FrameSetup' : {'StormFrames' : [500, 19999],\
                     'StackFrames' : [  1, 10000]},\
     'OutputName' : 'stormdata_corrected.csv',\
     'figname':'upscaling.png'}

A6 = \
    {'ImageStackName' : 'stack.tif',\
     'StormFileName' : 'storm.csv',\
     'BeadFileNames' : ['beads/bead1.csv', 'beads/bead2.csv',\
                        'beads/bead3.csv'],\
     'FrameGrouping' : 20,\
     'BeadStabilizing' : 10,\
     'StackUpscaling' : None,\
     'StackSmoothing' : True,\
     'StackCropping' : None,\
     'StackStabilizing' : 10,\
     'StackNmPerPixel' : 100,\
     'FrameSetup' : {'StormFrames' : [500, 19999],\
                     'StackFrames' : [  1, 10000]},\
     'OutputName' : 'stormdata_corrected.csv',\
     'figname':'nocropping.png'}

for a in [A1,A2,A4,A5,A6]:
 print(a)
 for i in range(8):
  plt.figure(i)
  plt.clf()
 A = Anchor.Analyser(a)
 
 A.startAnalysis()
 #plt.plot()
 plt.savefig('drift_'+a['figname'])

 A.tStorm.plotMeanDriftVsXcorr(A.im_stack)
 plt.savefig('vs_'+a['figname'])


 A.tStorm.plotError(A.im_stack)
 #plt.plot()
 plt.savefig('error_'+a['figname'])
 del(A)
