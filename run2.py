import logging, sys
logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)

import os
if not 'DISPLAY' in os.environ.keys():
    import matplotlib
    matplotlib.use('Agg')
from matplotlib import pyplot as plt

from Analyser import *

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

A = Analyser(AnalysisParameters)
A.startAnalysis()
plt.savefig('1.png')

A.tStorm.plotDrifts()
plt.savefig('2.png')

A.tStorm.plotMeanDriftVsXcorr(A.im_stack)
plt.savefig('3.png')

A.tStorm.plotError(A.im_stack)
plt.savefig('4.png')