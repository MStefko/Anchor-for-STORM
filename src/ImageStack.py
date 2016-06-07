import libtiff
import os
import copy
import logging, sys
from operator import add
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

class ImageStack:
    """ A container class for an image stack with all the necessary
    tools to determine drift """
    def __init__(self, path = None):
    # Initialization of the class.
        self.originalImages = []
        self.images = []
        self.upscaling = None
        self.averaging = None
        self.cropping = None
        self.smoothing = None
        self.stabilizing = None
        self.nmPerPixel = 100


        self.path = path
        self.drift = []

    def importImages(self, path = None, skip = None, startsWith = None, saveOriginal = None):
        if not path:
            path = self.path
        if path.lower().endswith(('.tif', '.tiff')):
            self.importImagesTif(path = path, skip = skip, saveOriginal = saveOriginal)
        else:
            self.importImagesPng(path = path, skip = skip, saveOriginal = saveOriginal, startsWith = startsWith)
        return

    def importImagesPng(self, pathToFolder = None, skip = None, startsWith = None, saveOriginal = None):
        """ Import images from a folder containing the .png images.
        After finishing, the image stack is loaded to self.images
        Parameters:
         - pathToFolder (string): 
                Path to desired folder.
                If empty, current working directory is used.
         - skip (int):
                If set to x, only each x-th image will be loaded.
                Useful for working with large datasets.
         - startsWith (string):
                Only loads images which have a filename starting with <startsWith>
         - saveOriginal (bool):
                Also copies the image stack to self.originalImages, so you can reload it easier later. """
        logging.info("Started import.")
        if not pathToFolder:
            pathToFolder = os.getcwd()
        else:
            pathToFolder = os.path.normcase(pathToFolder)
        self.path = pathToFolder
        filenames = [f for f in os.listdir(pathToFolder) if f.endswith(".png")]
        filenames.sort()
        if startsWith:
            filenames = [f for f in filenames if f.startswith(startsWith)]
        if skip:
            filenames = filenames[0: :skip]

        self.images = [misc.imread(os.path.join(pathToFolder,name)) for name in filenames]
        if saveOriginal:
            self.originalImages = copy.deepcopy(self.images)
        logging.info("Finished import: " + str(len(self.images)) + " images.")
        return

    def importImagesTif(self, path = None, skip = None, saveOriginal = None):
        """ Import from a .tif image stack. """
        logging.info("Started import.")
        if not path:
            if not self.path:
                raise ValueError('Path for import not set!')
            path = self.path
        f = libtiff.TiffFile(path)
        self.images = [a.get_image().astype('float') for a in f.get_tiff_array()]
        if saveOriginal:
            logging.info("Started copying.")
            self.originalImages = copy.deepcopy(self.images)
            logging.info("Finished copying.")
        logging.info("Finished import: " + str(len(self.images)) + " images.")
        return

    def _smoothImages(self):
        """ Applies a 2D gaussian smoothing filter to the images.
        Sigma: 0.5 pixels"""
        logging.info("Smoothing images.")
        smoothed = [gaussian_filter(im, 0.5) for im in self.images]
        self.images = smoothed
        return

    def _normalizeImages(self):
        """ Each image is individually centered around the value 0.0, and then normalized so that its
        standard deviation is 1.0. """
        logging.info("Normalizing images.")
        normalized = [((im-np.mean(im))/np.std(im)) for im in self.images]
        self.images = normalized
        return

    def _makeImagesSquare(self):
        """ Each image is cropped so that its dimensions are square. The top-left corner is preserved. """
        logging.info("Cropping images to square shape.")
        dims = self.images[0].shape
        mindim = min(dims)
        squared = []
        for im in self.images:
            squared += [np.copy(im[0:mindim,0:mindim])]
        self.images = squared
        return

    def _upscaleImages(self, interpolationType='linear'):
        """ Upscale images using the interpolation method defined in <interpolationType>
        Parameters:
         - interpolationType (string): ['linear', 'cubic', 'quintic'] """
        logging.info("Upscaling images by factor of: " + str(self.upscaling))
        assert(self.upscaling>1.0)
        scale = self.upscaling
        upsampled = []

        for im in self.images:
            x_old = np.array(range(im.shape[0]))
            y_old = np.array(range(im.shape[1]))
            f_upsamp = interp2d(x_old, y_old, im, kind = interpolationType, copy = False)
            x_new = np.arange(0.0, float(im.shape[0]), 1.0/scale)
            y_new = np.arange(0.0, float(im.shape[1]), 1.0/scale)
            upsampled += [f_upsamp(x_new, y_new)]

        self.images = upsampled
        return

    def _averageImages(self):
        """ Puts the images into groups of size defined in <self.averaging>.
        Images of each group are then averaged together, and form one resulting image.
        Size of the stack thus gets reduced by a factor of <self.averaging>. """
        logging.info("Grouping and averaging images in groups of " + str(self.averaging))
        assert(self.averaging)
        group = self.averaging
        averaged = []

        incr = np.zeros(self.images[0].shape)
        count = 0; i = 0
        while (i<len(self.images)):
            incr = np.add(self.images[i], incr)
            count+=1; i+=1
            if ((i%group)==0):
                averaged += [np.divide(incr,count)]
                incr = np.zeros(self.images[0].shape)
                count = 0
        if count:
            averaged += [np.divide(incr,count)]
        self.images = averaged
        return

    def _getCrossCorrelation(self, ref, mask, fft_ref = False, fft_mask = False):
        """ Computes the cross correlation between reference and mask images.
        For parameter description, refer to <self._getDriftValue()> """
        # Images should be square and of same dimensions at this point.
        assert(ref.shape==mask.shape)

        if not fft_ref:
            four_ref = fft2(ref)
        else:
            four_ref = ref

        if not fft_mask:
            # Crop the mask and replace the edges with 0.0 values. Helps the crosscorrelation.
            if self.cropping:
                size = min(mask.shape)
                crop = self.cropping
                mask_cropped = np.copy(mask[crop:(size-crop), crop:(size-crop)])
                mask_padded = np.pad(mask_cropped, crop, mode='constant')
                four_mask = fft2(mask_padded)
            else:
                four_mask = fft2(mask)           
        else:
            four_mask = mask

        # Conjugate the mask.
        four_mask_conj = np.conjugate(four_mask)
        # Compute pointwise product of reference and mask.
        product = np.multiply(four_mask_conj, four_ref)
        # Compute ifft of this product
        xcorr = ifft2(product)
        # Take the absolute value
        xcorr_abs = np.absolute(xcorr)
        return xcorr_abs

    def _getDriftValue(self, ref, mask, fft_ref = False, fft_mask = False):
        """ 
        Takes a reference and mask image, and returns the computed drift in nm.
        Parameters:
         - ref (square np.array): Either a raw image, or already computed fft of the reference image.
            If the fft has been computed, the flag <fft_ref> must be set to True.
         - mask (square np.array, same dimensions as ref): Same as for <ref>. If the fft of the images
            is passed, <fft_mask> must be set to True.
         - fft_ref (bool): Flag to mark if the fft has already been computed for <ref>.
         - fft_mask (bool): Ditto for <mask>. """
        # Compute the crosscorrelation. The -1.0 is to use the minimize() function, since
        # we are looking for the maximum of the crosscorrelation.
        xcorr = -1.0 * self._getCrossCorrelation(ref, mask, fft_ref = fft_ref, fft_mask = fft_mask)
        # Get the approximate location of maximum by just looking and the largest value.
        approx = np.unravel_index(np.argmin(xcorr), xcorr.shape)
        # Make a 3x3 periodic tile of the xcorrelation array, and work with the center one.
        # This is convenient, so that when the maximum is close to the edge, we don't have
        # to worry about going out of bounds.
        tiled = np.tile(xcorr, (3,3))#                                  X X X
        # Translate the maximum coordinates into the middle tile.  O -> X O X
        approx = np.add(xcorr.shape, approx)#                           X X X
        # Get the coordinate ranges of the array.
        x = np.arange(tiled.shape[0])
        y = np.arange(tiled.shape[1])

        # Set low and high bounds of interpolation.
        # Low bound: 5 pixels under approximate estimate.
        x_l = approx[0]-5; y_l = approx[1]-5
        # High bound: 5 pixels above approximation.
        x_h = approx[0]+5; y_h = approx[1]+5
        # This gives us a 10x10 pixel region to interpolate.
        f_upsamp = RectBivariateSpline(x[x_l:x_h], y[y_l:y_h], tiled[x_l:x_h,y_l:y_h])

        # Find the minimum of the interpolation.
        minimum_loc = minimize(lambda v: f_upsamp(v[0],v[1]), approx, bounds = ((x_l, x_h), (y_l, y_h))).x
        # Subtract the dimensions of the image, so that we move again from the
        # 3x3 tiled array into ordinary coordinates. Note the -1.0 .
        minimum_loc = np.add(-1.0 * minimum_loc, np.array(xcorr.shape))
        # Transform from (upscaled) pixel coordinates to nanometers.
        nm = self.nmPerPixel
        scale = self.upscaling if self.upscaling else 1.0
        drift_nm = [float(a)*nm/scale for a in minimum_loc]
        logging.debug(drift_nm)
        # Reverse the order so that x is first and y second. 
        # For some reason it otherwise returns as (y,x). 
        return drift_nm[::-1] 

    def getDriftValues(self):
        """ Wrapper function for the drift determination algorithm.
        The first image of the stack is picked as the reference. """
        logging.info("Started Xcorrelation.")
        # Take first image of the stack as the reference, and pre-compute the fft2, so that
        # we don't have to compute it at each call of _getDriftValue(). Remember to set the
        # <fft_ref> flag to True.
        ref = fft2(self.images[0])
        # Take the rest of the stack as masks.
        masks = self.images[1:]
        # Drift of the first image is 0 (duh)
        drifts = [(1, 0.0, 0.0)]
        # We use i to save the position in the stack together with the drift values.
        i = 2
        for mask in masks:
            # Call the function for each mask.
            drifts += [ (i,) + tuple(self._getDriftValue(ref, mask, fft_ref = True, fft_mask = False))]
            i += 1
        self.drift = drifts
        logging.info("Finished Xcorrelation.")
        return

    def getDriftValuesStep(self):
        """ Instead of choosing the first image as reference, here we always compare the neighboring
        two images. """
        drifts = [(1, 0.0, 0.0)]
        for i in range(len(self.images)-1):
            # The incremental drift needs to be added to the previous one.
            drifts += [ (i+2,) + tuple(map(add, self._getDriftValue(self.images[i], self.images[i+1]), drifts[-1][1:])) ]

        self.drift = drifts
        return

    def getDrift(self):
        if not self.drift:
            self.getDriftValues()
        return self.drift

    def _stabilizeDriftValues(self):
        if not self.stabilizing:
            raise ValueError('Stabilizing factor not set.')
        self.getDrift()
        frames = self.stabilizing
        x_shift = np.mean( [a[1] for a in self.drift[0:frames] ] )
        y_shift = np.mean( [a[2] for a in self.drift[0:frames] ] )
        correctedX = [(a[1] - x_shift) for a in self.drift]
        correctedY = [(a[2] - y_shift) for a in self.drift]

        newDrift = []
        for i in range(len(self.drift)):
            newDrift += [ (self.drift[i][0], correctedX[i], correctedY[i])]
        self.drift = newDrift
        return

    def displayDrift(self, poly = None):
        """ Convenient way to plot the drift.
        Parameters:
         - poly (int): Optional degree of polynomial interpolation. """
        xs = np.array([d[1] for d in self.drift])
        ys = np.array([d[2] for d in self.drift])
        # Proper scaling of the x axis
        skip = self.averaging if self.averaging else 1.0        
        n = skip * np.arange(len(xs))
        if poly:
            px = np.polyfit(n,xs,poly)
            py = np.polyfit(n,ys,poly)
            lines = plt.plot(n,xs,'ro',n,ys,'bo',\
                n,np.polyval(px,n),'g',n,np.polyval(py,n),'y')
            plt.setp(lines, 'linewidth', 3.0)
        else:
            plt.plot(n,xs,'r.',label='X')
            plt.plot(n,ys,'b.',label='Y')
        plt.xlabel('Frame')
        plt.ylabel('Drift [nm]')
        plt.title('Xcorrelation')
        plt.legend()
        #plt.show()
        #return lines

    def setDefaultSettings(self):
        self.upscaling = 1
        self.averaging = 20
        self.cropping = 32

    def analyse(self):
        """ Default analysis algorithm. """
        if self.averaging:
            self._averageImages()
        self._makeImagesSquare()
        if self.smoothing:
            self._smoothImages()
        self._normalizeImages()
        if self.upscaling > 1:
            self._upscaleImages()
        self.getDriftValues()
        if self.stabilizing:
            self._stabilizeDriftValues()
        self.displayDrift()
