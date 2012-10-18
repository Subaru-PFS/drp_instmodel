#!/usr/bin/env python

import pydebug    

import os
import pyfits
import cPickle as pickle

import numpy as np

import scipy.ndimage
import scipy.signal
import scipy.interpolate as spInterp

import psf

class SplinedPsf(psf.Psf):
    def __init__(self, detector, spotType='jeg'):
        """ Create or read in our persisted form. By default use JEG's models. 

        Parameters
        ----------

        detector : a Detector object
           carries the properties of the detector (and camera).
        spot_type : str, optional
           Whether to load 'jeg' or 'zemax' spot images.
        """

        psf.Psf.__init__(self, detector)
        
        # The locations at which we have PSFs
        self.wave = []
        self.fiber = []
        self.spots = []

        self.coeffs = []
        self.xcCoeffs = []
        self.ycCoeffs = []

        self.loadFromFile(spotType)

    def __str__(self):
        nSpots = len(self.spots)
        if nSpots > 0:
            return ("%return<%d spots; %d wavelengths (%0.2fAA to %0.2fAA); %d fibers>" %
                    (self.__class__.__name__, nSpots,
                     len(np.unique(self.wave)), np.min(self.wave), np.max(self.wave),
                     len(np.unique(self.fiber))))
        else:
            return ("%s<0 spots>" %
                    (self.__class__.__name__))

    @staticmethod
    def psfFileDir(band, spotType):
        dataRoot = os.environ.get('PFS_INSTDATA_DIR', '.')
        return os.path.join(dataRoot, 'data', 'spots', spotType, band)
        
    @staticmethod
    def psfSpotFile(band, spotType):
        return os.path.join(SplinedPsf.psfFileDir(band, spotType),
                            'psfSplines.pck')
        
    @property 
    def imshape(self):
        return self.spots[0].shape

    def traceCenters(self, fibers, waves):
        """ Return the pixel centers for the given fibers and wavelengths """

        fibers = np.array(fibers)
        waves = np.array(waves)
        
        fiberIdx = np.argsort(fibers)
        waveIdx = np.argsort(waves)
             
        x = self.xcCoeffs(fibers[fiberIdx], waves[waveIdx])
        y = self.ycCoeffs(fibers[fiberIdx], waves[waveIdx])

        return (x[fiberIdx][:, waveIdx], 
                y[fiberIdx][:, waveIdx])
    
    def wavesForRows(self, fibers, rows=None, waveRange=None, pixelScale=None):
        """ Return our best estimate for the wavelength at the given row centers.

        Returns:
          rows   - the rows which we evaluated. [NROWS]
          waves  - for each fiber, the wavelengths at rows [NFIBERS, NROWS]
    
        """

        if pixelScale == None:
            pixelScale = self.detector.config['pixelScale']            
        
        if waveRange == None:
            waveRange = self.wave.min(), self.wave.max()

        # Assume that the full spectrum fits on the detector.
        minY, maxY = self.ycCoeffs(fibers, waveRange)[0]
        doReorder = minY > maxY
        if doReorder:
            minY, maxY = maxY, minY
        if minY < 0:
            print("one of wavelengths %s maps below the deector (%0.5f mm)" % (waveRange, minY))
        if maxY >self.detector.config['ccdSize'][1]:
            print("one of wavelengths %s maps above the deector (%0.5f mm)" % (waveRange, maxY))
        
        minRow = int(minY/pixelScale)
        maxRow = int(maxY/pixelScale)

        if rows == None:
            rows = np.arange(minRow, maxRow+1)

        # Invert the spline into a row->wave map. 
        # Just use a linear interpolation based on the evaluation near the pixels.
        allWaves = np.linspace(waveRange[0], waveRange[1], maxRow-minRow)

        waves = []
        for f in fibers:
            allWaveRows = self.ycCoeffs([f], allWaves)[0] / pixelScale

            if doReorder:
                allWaveRows0 = allWaveRows[::-1]
                allWaves0 = allWaves[::-1]
            else:
                allWaveRows0 = allWaveRows
                allWaves0 = allWaves
                
            waveFunc = spInterp.interp1d(allWaveRows0, allWaves0, 'linear', bounds_error=False)
            fiberWaves = waveFunc(rows[1:])
            waves.append(fiberWaves)

        return rows, waves
    
    def psfsAt(self, fibers, waves=None, everyNthPsf=1, usePsfs=None):
        """ Return a stack of PSFs, instantiated on a rectangular grid.

        Paramaters
        ----------
        fibers : array_like, int
           the fiber IDs we want PSFs for
        waves : array_like, AA, optional
           the wavelengths we want PSFs as. If not set, uses the native wavelength grid.
        everyNthPsf : bool, optional
           how many wavelengths get the same PSF -- this is a performance optimization.  The
           default is for each wavelength to get its own PSF.
        usePsfs : array_like , optional
           
        Returns
        -------
        psfImages : array
           a 3D array of images [len(fibers)*len(waves), imshape.x, imshape.y]
        psfIds : array
           the (fiberId, wavelength) positions the images are instantiated for.
        psfPos : array
           the (x,y) positions the PSF images should be centered on.

        Notes
        -----
        The everyNthPsf arg should be waveSpacing (and fiberSpacing) or something.
        
        """

        if waves == None:
            waves = np.unique(self.wave)
            
        waveSign = 1 if waves[-1] > waves[0] else -1
        interpWaves = waves[::waveSign*everyNthPsf]
        
        centers = [(x,y) for x in fibers for y in waves]
        if usePsfs != None:
            newImages = usePsfs
        else:
            newImages = np.zeros(shape=(len(fibers)*len(interpWaves),
                                        self.imshape[0],
                                        self.imshape[1]))
            for ix in range(self.imshape[0]):
                for iy in range(self.imshape[1]):
                    newImages[:, iy, ix] = self.coeffs[iy, ix](fibers, interpWaves).flat

        # OK, optionally replicate the sparsely sampled PSF grid into the full one.
        if everyNthPsf > 1:
            newImageList = []
            for w_i in range(len(interpWaves)):
                for i_i in range(0, everyNthPsf):
                    newImageList.append(newImages[w_i])
            finalImages = newImageList[:len(waves)+1]
        else:
            finalImages = newImages

        print "psfsAt: for %d waves and useNthPsf=%d, returned %d unique psfs" % (len(waves),
                                                                                  everyNthPsf,
                                                                                  len(interpWaves))
        return finalImages, centers, self.traceCenters(fibers, waves)

    def makeComb(self, waves, nth=1, hackScale=1000):
        """ Return a functor which returns True for every nth item in the initialization list. """
        class Comb(object):
            def __init__(self, waves, nth):
                self.waves = np.asarray(waves)

            def __call__(self, x):
                return hackScale * (x in self.waves)

            def getNativeValues(self):
                return np.asarray((self.waves, 
                                   hackScale * np.ones(len(self.waves))),
                                   dtype='f4').T
            
        return Comb(waves, nth)
    
    def fiberImage(self, fiber, spectrum, outImg=None, waveRange=None, everyNthPsf=1, returnUnbinned=False):
        """ Return an interpolated image of a fiber """

        # Evaluate at highest resolution
        pixelScale = self.spotScale
        everyNthPsf *= int(self.detector.config['pixelScale'] / self.spotScale)
        
        if waveRange == None:
            waveRange = self.wave.min(), self.wave.max()
        minX, maxX = self.xcCoeffs([fiber], waveRange)[0]
        minY, maxY = self.ycCoeffs([fiber], waveRange)[0]

        minRow = minY/pixelScale
        maxRow = maxY/pixelScale
        minCol = minX/pixelScale

        # Generalize this... CPL
        if minRow > maxRow:
            minRow, maxRow = maxRow, minRow

        # Get the wavelengths for the fiber pixels.
        allPixelRows, allPixelWaves = self.wavesForRows([fiber], waveRange=waveRange, 
                                                        pixelScale=pixelScale)
        pixelWaves = allPixelWaves[0]
        pixelFlux = spectrum(pixelWaves)
        
        # Get the PSFs and their locations on the oversampled pixels.
        fiberPsfs, psfIds, centers = self.psfsAt([fiber], pixelWaves, everyNthPsf=everyNthPsf)
        xCenters, yCenters = [c[0] for c in centers]
        
        psfToSpotRatio = self.detector.config['pixelScale'] / pixelScale
        outImgSpotPixelScale = round(psfToSpotRatio)

        # pixels
        traceWidth = int((xCenters.max() - xCenters.min())/pixelScale + 0.5)
        traceHeight = int((yCenters.max() - yCenters.min())/pixelScale + 0.5)
        spotWidth = fiberPsfs[0].shape[-1]

        # We want out fiber image to be sized in units of outImgSpotPixelScale.
        fiHeight = (traceHeight + spotWidth)
        fiWidth = (traceWidth + spotWidth)
        fiHeight += outImgSpotPixelScale - fiHeight%outImgSpotPixelScale
        fiWidth += outImgSpotPixelScale - fiWidth%outImgSpotPixelScale
        fiberImage = np.zeros((fiHeight, fiWidth), dtype='f4')

        # mm
        fiberImageOffset = np.asarray((xCenters.min(), yCenters.min()))

        # pixels
        outImgOffset = fiberImageOffset / self.detector.config['pixelScale']

        # Adjust offset by fractional detector pixel
        fiberImageOffset -= (outImgOffset - outImgOffset.astype('i4')) * self.detector.config['pixelScale']
        outImgOffset = outImgOffset.astype('i4')
        
        if outImg == None:
            outExp = self.detector.simBias()
            outImg = outExp.image

        # construct the oversampled fiber image
        lasty = 0
        for i in range(len(pixelWaves)):
            specWave = pixelWaves[i]
            specFlux = pixelFlux[i]
            
            rawPsf = fiberPsfs[i]

            # in mm
            xc = xCenters[i]
            yc = yCenters[i]
            xoffset = xc - fiberImageOffset[0]
            yoffset = yc - fiberImageOffset[1]

            # pix offset
            xPixOffset = xoffset / pixelScale
            yPixOffset = yoffset / pixelScale

            # Keep the shift to the smallest fraction possible, or rather keep the integer steps 
            # exactly +/- 1.
            inty = round(yPixOffset)
            fracy = yPixOffset - inty

            intx = round(xPixOffset)

            if i % 1000 in range(2):
                print("%5d %6.1f xc: %3.4f yc: %3.4f %3.10f %d %d %3.10f %0.1f" % (i, specWave, xc, yc, yPixOffset, lasty, inty, fracy, specFlux))
            lasty = inty
            
            # Assume we are well enough oversampled to ignore fractional pixel shifts.
            spot = specFlux * rawPsf
            
            self.placeSubimage(fiberImage, spot, intx, inty)
                        
            # bin psf to ccd pixels, shift by fractional pixel only.
            #psf = self.scalePsf(rawPsf, -fracx, -fracy, doDetails=doDetails)
            #self.placeSubimage(outImg, spot, intx, inty)

        # transfer flux from oversampled fiber image to final resolution output image
        resampledFiber = self.addOversampledImage(fiberImage, outImg, outImgOffset, outImgSpotPixelScale)

        if returnUnbinned:
            return outImg, minRow, minCol, fiberImage, resampledFiber
        else:
            return outImg, minRow, minCol

    def addOversampledImage(self, inImg, outImg, outOffset, outScale):
        """ Add the outScale-oversampled inImg to outImg at the given offset. """

        resampled = self.rebin(inImg, inImg.shape[0]/outScale, inImg.shape[1]/outScale)
        self.placeSubimage(outImg, resampled, *outOffset)

        return resampled
    
    def rebin(self, a, *args):
        """ rebin(a, *new_axis_sizes) taken from scip cookbook. """
        
        shape = a.shape
        lenShape = len(shape)
        factor = np.asarray(shape)/np.asarray(args)
        evList = ['a.reshape('] + \
                 ['args[%d],factor[%d],'%(i,i) for i in range(lenShape)] + \
                 [')'] + ['.sum(%d)'%(i+1) for i in range(lenShape)]

        return eval(''.join(evList))
    
    def placeSubimage(self, img, subImg, xc, yc):
        parentx, childx = self.trimSpan((0, img.shape[1]-1),
                                        (0, subImg.shape[1]-1),
                                        xc)
        
        parenty, childy = self.trimSpan((0, img.shape[0]-1),
                                        (0, subImg.shape[0]-1),
                                        yc)

        try:
            img[parenty, parentx] += subImg[childy, childx]
        except ValueError, e:
            print "failed to place child at (%d,%d): %s" % (xc,yc, e)
    
    def trimSpan(self, _parent, _child, childOffset=None):
        """ Return the final overlapping span as slices. """ 

        parent = list(_parent)
        child = list(_child)
        
        # Put child on parent's frame
        if childOffset:
            child[0] += childOffset
            child[1] += childOffset

        # Trim child to parent frame
        if child[0] < parent[0]:
            child[0] = parent[0]
        if child[1] > parent[1]:
            child[1] = parent[1]

        # Trim parent to child frame
        if parent[0] < child[0]:
            parent[0] = child[0]
        if parent[1] > child[1]:
            parent[1] = child[1]

        # Optionally put child back in its own frame.
        if childOffset:
            child[0] -= childOffset
            child[1] -= childOffset

        if parent[1]-parent[0] != child[1]-child[0]:
            pydebug.debug_here()
            
        return (slice(parent[0], parent[1]+1),
                slice(child[0], child[1]+1))

    def _frac(self, n):
        """ Return the fractional part of a float. """

        return np.modf(n)[0]
    
    def _sincfunc(self, x, dx):
        dampfac = 3.25
        if dx != 0.0:
            return np.exp( -((x+dx)/dampfac)**2 ) * np.sin( np.pi*(x+dx) ) / (np.pi * (x+dx))
        else:
            xx = np.zeros(len(x))
            xx[len(x)/2] = 1.0
            return xx

    def _sincshift(self, image, dx, dy):
        """ UNUSED & UNTESTED from boss -- for comparison only, etc. etc. """
        
        sincrad = 10
        s = np.arange(-sincrad, sincrad+1)
        sincx = self._sincfunc(s, dx)

        #- If we're shifting just in x, do faster 1D convolution with wraparound                                
	    #- WARNING: can introduce edge effects if PSF isn't nearly 0 at edge                                    
        if abs(dy) < 1e-6:
            newimage = scipy.signal.convolve(image.ravel(), sincx, mode='same')
            return newimage.reshape(image.shape)

        sincy = self._sincfunc(s, dy)
        kernel = np.outer(sincy, sincx)
        newimage = scipy.signal.convolve2d(image, kernel, mode='same')
        return newimage

    def loadFromFile(self, spotType, filename=None):
        if not filename:
            filename = SplinedPsf.psfSpotFile(self.detector.band, spotType)

        with open(filename, 'r') as pfile:
            d = pickle.load(pfile)

        self.spots = d['spots']
        self.wave = d['wave']
        self.fiber = d['fiber']
        self.coeffs = d['coeffs']
        self.spotScale = d['spotScale']

        # Make fiber -> (xc,yc) and wave->(xc,yc) splines. We have to do this here, because
        # the detector geometry matters.
        
        # Shift offsets to origin.
        xc = d['xc'] + self.detector.config['ccdSize'][1] * self.detector.config['pixelScale'] / 2
        yc = d['yc'] + self.detector.config['ccdSize'][0] * self.detector.config['pixelScale'] / 2

        # XXX - This fiber/wave indexing scheme is not safe in general
        xx = np.unique(self.fiber)
        yy = np.unique(self.wave)

        xcCoeffs = spInterp.RectBivariateSpline(xx, yy, xc.reshape(len(xx), len(yy)))
        ycCoeffs = spInterp.RectBivariateSpline(xx, yy, yc.reshape(len(xx), len(yy)))
        
        print "WARNING, swapping x and y centers: spectra disperse along columns, but the optomechanical view is that dispersion is along X"
        self.xc = d['yc']
        self.yc = d['xc']
        self.xcCoeffs = ycCoeffs
        self.ycCoeffs = xcCoeffs

def constructSplinesFromSpots(band, spotType='zemax'):
    """ NASTY function to construct a pickle file holding the per-pixel spot splines. These should probably be 
    in the FITS file, but I haven't worked out how to do that (.knots -> binary table?)
    """
    
    psfFilepath = SplinedPsf.psfSpotFile(band, spotType)

    spotsFilepath = os.path.join(SplinedPsf.psfFileDir(band, spotType),
                                 'spots.fits')
    
    sf = pyfits.open(spotsFilepath, mode='readonly')
    s = sf[1].data
    
    outputDict = {}
    outputDict['wave'] = s['wavelength']
    outputDict['fiber'] = s['fiberIdx']
    outputDict['spotScale'] = sf[0].header['pixscale']
    outputDict['xc'] = s['spot_xc']
    outputDict['yc'] = s['spot_yc']

    spots = s['spot'].astype('float32')
    outputDict['spots'] = spots

    imshape = spots.shape[1:]

    # Make a spline for each pixel. The spots are centered on the output grid,
    # and we track the xc,yc offset separately.

    # XXX - This fiber/wave indexing scheme is not safe in general
    xx = np.unique(outputDict['fiber'])
    yy = np.unique(outputDict['wave'])

    coeffs = np.zeros(imshape, dtype='O')
    for ix in range(imshape[0]):
        print "splining col %d" % (ix)
        for iy in range(imshape[1]):
            coeffs[iy, ix] = spInterp.RectBivariateSpline(xx, yy,
                                                          spots[:, iy, ix].reshape(len(xx), len(yy)))
    outputDict['coeffs'] = coeffs
    
    with open(psfFilepath, 'w+') as pfile:
        print "pickling...."
        pickle.dump(outputDict, pfile, protocol=-1)

