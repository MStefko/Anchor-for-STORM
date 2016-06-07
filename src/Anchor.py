from ImageStack import *
from tSTORMdata import *
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
from STORM import *

import logging, sys

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


class Analyser:
    """ Encapsulator class, which loads 3 types of data:
         - im_stack: A .tif (or .png) stack of bright-field images for cross-correlation.
         - tStorm: A set of .csv files containing trajectories of beads tracked by ThunderSTORM
         - STORM: Localization data of fluorescent markers.
        Analysis can then be carried out on this data, and corrected ouput saved to a file. """

    def __init__(self, im_stack = None, tStorm = None, STORM = None):
        self.im_stack = im_stack # Class ImageStack()
        self.tStorm = tStorm     # Class tSTORMdata()
        self.STORM = STORM       # Class STORM()
        self.f_x = None   # Interpolation function for X coordinate.
        self.f_y = None   # Interpolation function for Y coordinate.
        self._poly = None # Carries information about conversion from im_stack to STORM frame IDs.
        self.outputName = None # Output of corrected data .csv file
        self.correctedData = None # pandas frame

    def __init__(self, Params):
        keys = Params.keys()
        if not 'ImageStackName' in keys:
            raise KeyError('.tif image stack is required for analysis. (Key: ImageStackName)')
        if not 'StormFileName' in keys:
            raise KeyError('STORM .csv file required. (Key: StormFileName)')

        if 'BeadFileNames' in keys:
            beadGrouping = Params['FrameGrouping'] if ('FrameGrouping' in keys) else None
            beadStabilizing = Params['BeadStabilizing'] if ('BeadStabilizing' in keys) else None
            self.tStorm = tSTORMdata(filenames = Params['BeadFileNames'],\
                grouping = beadGrouping, stabilizing = beadStabilizing)
        else:
            self.tStorm = None


        self.im_stack = ImageStack(path = Params['ImageStackName'])
        self.im_stack.averaging = Params['FrameGrouping'] if ('FrameGrouping' in keys) else None
        self.im_stack.upscaling = Params['StackUpscaling'] if ('StackUpscaling' in keys) else None
        self.im_stack.smoothing = Params['StackSmoothing'] if ('StackSmoothing' in keys) else None
        self.im_stack.cropping = Params['StackCropping'] if ('StackCropping' in keys) else None
        self.im_stack.stabilizing = Params['StackStabilizing'] if ('StackStabilizing' in keys) else None
        self.im_stack.nmPerPixel = Params['StackNmPerPixel'] if ('StackNmPerPixel' in keys) else 100

        self.STORM = STORM(filename = Params['StormFileName'])

        if 'FrameSetup' in keys:
            s = Params['FrameSetup']
            self.setConvRatios(stackframes = s['StackFrames'],
                               stackgrouping = Params['FrameGrouping'],
                               STORMframes = s['StormFrames'])

        self.outputName = Params['OutputName'] if ('OutputName' in keys) else None

        self.f_x = None   # Interpolation function for X coordinate.
        self.f_y = None 
        self.correctedData = None # pandas frame


    def setConvRatios(self, stackframes = None, stackgrouping = None, STORMframes = None):
        """ This function sets up the correct conversion relation between IDs of bright-field
            and STORM frames.
            Input:
            - stackframes <array of at least 2 ints>:
                Sequence numbers of at least 2 images from the .tif stack, starting from 1.
            - STORMframes <array of ints of same dimensions as stackframes>:
                Corresponding frame IDs in the STORM .csv file for each of the stack ID entries.
            - stackgrouping <int>:
                If grouping (i.e. averaging) is used on the .tif stack, this number scales the
                conversion properly. """
        if stackgrouping:
            # Rescale the stackframes properly.
            stackframes = ((np.array(stackframes)-1)/stackgrouping) + 1
        self._poly = np.polyfit(np.array(STORMframes),\
                                np.array(stackframes), 1)
        return

    def _getStackIdFromSTORMId(self, STORM_ids):
        """ Input STORM frame ID, get corresponding im_stack frame ID. """
        if self._poly is None:
            logging.error("Error: Conversion ratios not set.")
            return None
        return np.polyval(self._poly, np.array(STORM_ids))

    def _getSTORMIdFromStackId(self, stack_ids):
        """ Reverse of the above """
        # invert the polynomial: y=ax+b -> x=y/a-b/a
        p = np.array( [1./self._poly[0], -self._poly[1]/self._poly[0]] )
        return np.polyval(p, np.array(stack_ids))

    def _getStackDrift(self):
        return self.im_stack.getDrift()

    def calculateStackTrajectory(self):
        """ Sets up the linear interpolation of the im_stack (crosscorrelation)
            trajectory. This is then used to correct the STORM data. """
        logging.info("Setting up interpolation functions from stack drift data.")
        d = self._getStackDrift()
        ids = [a[0] for a in d]
        logging.debug("Stack real IDs:")
        logging.debug(ids)
        xs = [a[1] for a in d]
        ys = [a[2] for a in d]
        self.f_x = interp1d(ids, xs, kind='linear')
        self.f_y = interp1d(ids, ys, kind='linear')

    def _getDriftForSTORMId(self, ID):
        stackId = self._getStackIdFromSTORMId(ID)
        return (self.f_x(stackId), self.f_y(stackId))

    def getCorrectedData(self):
        """ Heart of the algorithm, calculates the trajectory, and outputs
            corrected STORM data. """
        if self.f_x is None:
            # interpolate stack drift trajectory
            self.calculateStackTrajectory()
        # Load stormIDs and convert them to stackIDs, which are then plugged into f_x and f_y.
        storm_IDs = self.STORM.data['frame'].as_matrix()
        stack_IDs = self._getStackIdFromSTORMId(storm_IDs)
        logging.info("Calculating corresponding stack IDs from STORM IDs using"+\
            " the interpolation functions.")
        logging.debug("Stack new IDs:")
        logging.debug(stack_IDs)
        try:
            # Calculate the trajectory.
            diff_x = self.f_x(stack_IDs)
            diff_y = self.f_y(stack_IDs)
        except ValueError:
            """ If any STORM frames reach out of bounds of known drift trajectory,
                print out a warning that the trajectories have to be extrapolated.
                For edge cases (i.e. only outlying by 1 or 2 frames), this should
                not be a problem. Otherwise, something is wrong. """

            logging.warning("STORM frame(s) out of bounds of cross-correlation!")
            logging.warning("These STORM frames are out of bounds:")
            wrongframes = [a for a in np.unique(storm_IDs) if (self._getStackIdFromSTORMId(a)<min(self.f_x.x)\
                                            or self._getStackIdFromSTORMId(a)>max(self.f_x.x))]
            logging.warning(wrongframes)
            logging.warning("Their stack equivalents with current settings:")
            logging.warning(self._getStackIdFromSTORMId(wrongframes))
            logging.warning("Current bounds:")
            logging.warning([min(self.f_x.x), max(self.f_x.x)])
            logging.warning("Drift values for these frames will now be extrapolated."+\
                "If their values are close to bounds, no action is necessary.\n")
            self.f_x.bounds_error = self.f_y.bounds_error = False
            self.f_x.fill_value = self.f_y.fill_value = 'extrapolate'
            diff_x = self.f_x(stack_IDs)
            diff_y = self.f_y(stack_IDs)

        # Get old X and Y coordinates of all localizations.
        X_old = self.STORM.data['x'].as_matrix()
        Y_old = self.STORM.data['y'].as_matrix()
        assert(len(X_old)==len(diff_x))
        # Subtract the corresponding drift trajectories to get corrected values.
        X_new = pd.Series(np.add(X_old, -1.0*diff_x))
        Y_new = pd.Series(np.add(Y_old, -1.0*diff_y))

        # Copy the STORM dataset and overwrite the X and Y values.
        newdata = copy.copy(self.STORM.data)
        newdata['x'] = X_new
        newdata['y'] = Y_new
        return newdata

    def startAnalysis(self):
        self.im_stack.importImages()
        self.im_stack.analyse()

        self.correctedData = self.getCorrectedData()
        if self.outputName:
            self.correctedData.to_csv(self.outputName)
            logging.info("Corrected data saved to: "+self.outputName)


    def saveCsvFromPd(self, data, filename = 'data.csv'):
        # Most useless function ever.
        data.to_csv(filename)
        return