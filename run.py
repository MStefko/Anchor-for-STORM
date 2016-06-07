import logging, sys
logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)

from ImageStack import *
from tSTORMdata import *
from STORM import *
from Analyser import *
import os

if not 'DISPLAY' in os.environ.keys():
    import matplotlib
    matplotlib.use('Agg')
from matplotlib import pyplot as plt

im_stack = ImageStack()
im_stack.importImagesTif('stack.tif')

im_stack.averaging = 20
im_stack.cropping = 32
im_stack.smoothing = True

im_stack.analyse()
plt.savefig('1.png')

tStorm = tSTORMdata(filenames = ['bead1.csv', 'bead2.csv', 'bead3.csv'],\
                    grouping = 20)
logging.info("Stabilizing data!")
tStorm.stabilizeData(frames = 10)

tStorm.plotDrifts()
plt.savefig('2.png')

tStorm.plotMeanDriftVsXcorr(im_stack)
plt.savefig('3.png')

tStorm.plotError(im_stack)
plt.savefig('4.png')

storm = STORM(filename = 'storm.csv')

A = Analyser(im_stack = im_stack,\
             tStorm = tStorm,\
             STORM = storm)
A.setConvRatios(stackframes = [1, 10000], stackgrouping = 20,\
                STORMframes = [500, 19999])
A.calculateStackTrajectory()

newdata = A.getCorrectedData()
newdata.to_csv('newdata.csv')

olddata = pd.read_csv('storm.csv')
plt.figure(10)
plt.plot(newdata['frame'], newdata['x']-olddata['x'],'r.')
plt.savefig('5.png')
