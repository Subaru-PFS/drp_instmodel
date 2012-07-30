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
    def __init__(self, band, detector, constructFromSpots=False):
        """ Read in our persisted form. Optionally that from the zemax spots. """

        self.band = band
        self.wave = []
        self.fiber = []
        self.spots = []
        self.coeffs = []
        self.xcCoeffs = []
        self.ycCoeffs = []
        self.detector = detector

        # Map self.band to filenames using some config file. For now, hardcode
        dataRoot = os.environ.get('PFS_INSTDATA_DIR', '.')
        psfFilepath = os.path.join(dataRoot, 'data', 'spots', band, 'psf.pck')

        if constructFromSpots:
            spotFilepath = os.path.join(dataRoot, 'data', 'spots', band, 'spots.fits')
            self.makeSplinesFromSpotFile(spotFilepath)
            self._saveToFile(psfFilepath)

        # We don't actually need to do this, but it seems like a good check.
        self.loadFromFile(psfFilepath)

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

    @property 
    def imshape(self):
        return self.spots[0].shape
    
    def makeSplinesFromSpotFile(self, filename):
        """ Load ourselves from the preprocessed spot FITS files the given directory. """
        
        sf = pyfits.open(filename, mode='readonly')
        s = sf[1].data

        self.wave = s['wavelength']
        self.fiber = s['fiberIdx']
        self.xc = s['spot_xc']
        self.yc = s['spot_yc']
        self.spots = s['spot'].astype('float64')
        self.spotScale = sf[0].header['pixscale']

        # Shift offsets to origin.
        self.xc += self.detector.config['ccdSize'][1] * self.detector.config['pixelScale'] / 2
        self.yc += self.detector.config['ccdSize'][0] * self.detector.config['pixelScale'] / 2
        
        self.coeffs = np.zeros(self.imshape,
                               dtype='O')

        # Make a spline for each pixel. The spots are centered on the output grid,
        # and we track the xc,yc offset separately.

        # XXX - This fiber/wave indexing scheme is not safe in general
        xx = np.unique(self.fiber)
        yy = np.unique(self.wave)

        for ix in range(self.imshape[0]):
            print "splining col %d" % (ix)
            for iy in range(self.imshape[1]):
                self.coeffs[iy, ix] = spInterp.RectBivariateSpline(xx, yy,
                                                                  self.spots[:, iy, ix].reshape(len(xx), len(yy)))
        self.xcCoeffs = spInterp.RectBivariateSpline(xx, yy, self.xc.reshape(len(xx), len(yy)))
        self.ycCoeffs = spInterp.RectBivariateSpline(xx, yy, self.yc.reshape(len(xx), len(yy)))


    def traceCenters(self, fibers, waves):
        """ Return the pixel centers for the given fibers and wavelengths """

        if waves[0] > waves[-1]:
            waveSlice = slice(None, None, -1)
        else:
            waveSlice = slice(None, None, 1)
            
        x = self.xcCoeffs(fibers, waves[waveSlice])
        y = self.ycCoeffs(fibers, waves[waveSlice])

        return x[waveSlice][0], y[waveSlice][0]
    
    def wavesForRows(self, fibers, rows=None):
        """ Return our best estimate for the wavelength at the given row centers """

        if rows == None:
            rows = np.arange(self.detector.config['ccdSize'][1])

        allWaves = np.linspace(self.wave.min(), self.wave.max(), self.detector.config['ccdSize'][1])
        allWaveRows = self.ycCoeffs(fibers, allWaves)[0] / self.detector.config['pixelScale']

        doReorder = allWaveRows[0] > allWaveRows[-1]
        if doReorder:
            allWaveRows = allWaveRows[::-1]
            allWaves = allWaves[::-1]
            
        waveFunc = spInterp.interp1d(allWaveRows, allWaves, 'linear', bounds_error=False)
        waves = waveFunc(rows)

        return waves
    
    def psfsAt(self, fibers, waves):
        """ Return a stack of PSFs, instantiated on a rectangular grid.

        Args:
           fibers, waves - two arrays of positions; we instatiate ourselves at the rectangular grid of these.

        Returns:
           - a 3D array of images [len(xs)*len(ys), imshape.x, imshape.y]
           - a list of the (fiberId, wavelength) positions the images are instantiated for.
           - a list of the (x,y) positions the images should be centered on.
        """

        newImages = np.zeros(shape=(len(fibers)*len(waves),
                                       self.imshape[0],
                                       self.imshape[1]))

        if waves[0] > waves[-1]:
            interpWaves = waves[::-1]
            interpSlice = slice(None, None, -1)
        else:
            interpWaves = waves
            interpSlice = slice(None)
            
        centers = [(x,y) for x in fibers for y in waves]
        for ix in range(self.imshape[0]):
            if ix % 16 == 0:
                print "fibers %s col %d" % (fibers, ix)
            for iy in range(self.imshape[1]):
                newImages[interpSlice, iy, ix] = self.coeffs[iy, ix](fibers, interpWaves).flat

        return newImages, centers, self.traceCenters(fibers, waves)

    def makeComb(self, ids, nth=100, hackScale=10000):
        """ Return a functor which returns True for every nth item in the initialization list. """
        class Comb(object):
            def __init__(self, ids, nth):
                self.waves = np.array([x[1] for x in ids])
                self.comb = self.waves[::nth].astype('i4').tolist()
                    
            def __call__(self, x):
                return hackScale * (int(x) in self.comb)

        return Comb(ids, nth=nth)
    
    def fiberImages(self, fibers, spectra=None, outImg=None):
        if outImg == None:
            outImg = self.detector.simBias().image

        if spectra == None:
            spectra = [None] * len(fibers)
        for i, fiber in enumerate(fibers):
            self.fiberImage(fiber, spectra[i], outImg=outImg)

        return outImg
    
    def fiberImage(self, fiber, spectrum=None, outImg=None):
        """ Return an interpolated image of a fiber """

        minX, maxX = self.xcCoeffs([fiber], [self.wave.min(), self.wave.max()])[0]
        minY, maxY = self.ycCoeffs([fiber], [self.wave.min(), self.wave.max()])[0]

        minRow = minY/self.detector.config['pixelScale']
        maxRow = maxY/self.detector.config['pixelScale']
        minCol = minX/self.detector.config['pixelScale']

        if minRow > maxRow:
            minRow, maxRow = maxRow, minRow
            
        # Get the wavelengths for the fiber pixels.
        allWaves = self.wavesForRows([fiber], np.arange(int(minRow+1), 
                                                        int(maxRow)))

        # Get the oversampled PSFs and their locations at the fiber pixels.
        fiberPsfs, psfIds, centers = self.psfsAt([fiber], allWaves)
        xCenters, yCenters = centers

        # With no input, use a comb spectrum
        if spectrum == None:
            spectrum = self.makeComb(psfIds, nth=20)

        traceWidth = int((xCenters.max() - xCenters.min())/self.detector.config['pixelScale'] + 0.5)
        traceHeight = int((yCenters.max() - yCenters.min())/self.detector.config['pixelScale'] + 0.5)

        psfPixelScale = int(self.detector.config['pixelScale'] / self.spotScale + 0.5)
        spotWidth = fiberPsfs[0].shape[-1] / psfPixelScale
        spotRadius = (spotWidth+1) / 2

        if outImg == None:
            outExp = self.detector.simBias(shape=(traceHeight, traceWidth + 2*spotWidth))
            outImg = outExp.image
            outImgOffset = (xCenters.min(), yCenters.min())
        else:
            outImgOffset = (0,0)
            
        for row in range(len(allWaves)):
            rawPsf = fiberPsfs[row]
            fiberID, wave = psfIds[row]

            # in mm
            xc = xCenters[row]
            yc = yCenters[row]

            # pix offset
            xPixOffset = (xc-outImgOffset[0]) / self.detector.config['pixelScale'] + spotRadius
            yPixOffset = (yc-outImgOffset[1]) / self.detector.config['pixelScale']

            # bin psf to ccd pixels, shift by fractional pixel only.
            psf = self.scalePsf(rawPsf, -self.frac(xPixOffset), -self.frac(yPixOffset))

            # Assumed to be a spline on the final resolution.
            flux = spectrum(wave)

            # When do we need to do this before the downsampling?
            spot = psf * flux

            self.placeSubimage(outImg, spot, int(xPixOffset), int(yPixOffset))

        return outImg, minRow, minCol+spotRadius

    def scalePsf(self, rawPsf, xc, yc, doRescaleFlux=False):
        """ Given an oversampled image scale it to the detector grid after shifting it. """

        psfPixelScale = self.detector.config['pixelScale'] / self.spotScale

        psf = scipy.ndimage.affine_transform(rawPsf, 
                                             psfPixelScale*np.eye(2),
                                             (yc*psfPixelScale, xc*psfPixelScale),
                                             output_shape=np.array(rawPsf.shape)/psfPixelScale)

        if doRescaleFlux:
            psf *= (psfPixelScale*psfPixelScale)
        return psf

    def placeSubimage(self, img, subImg, xc, yc):
        w, h = subImg.shape

        parentx, childx = self.trimSpan((0, img.shape[1]-1),
                                        (0, subImg.shape[1]-1),
                                        xc)
        
        parenty, childy = self.trimSpan((0, img.shape[0]-1),
                                        (0, subImg.shape[0]-1),
                                        yc)

        img[parenty, parentx] += subImg[childy, childx]
        
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

    def frac(self, n):
        """ Return the fractional part of a float. """
        return n-int(n)
    
    def _sincfunc(self, x, dx):
        dampfac = 3.25
        if dx != 0.0:
            return np.exp( -((x+dx)/dampfac)**2 ) * np.sin( np.pi*(x+dx) ) / (np.pi * (x+dx))
        else:
            xx = np.zeros(len(x))
            xx[len(x)/2] = 1.0
            return xx

    def _sincshift(self, image, dx, dy):
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

    def _saveToFile(self, filename):
        """ Saves our state (esp. the splines), currently to a pickle.  """

        with open(filename, 'w+') as pfile:
            d = {}

            d['coeffs'] = self.coeffs
            d['xcCoeffs'] = self.xcCoeffs
            d['ycCoeffs'] = self.ycCoeffs

            d['wave'] = self.wave
            d['fiber'] = self.fiber
            d['spots'] = self.spots

            d['spotScale'] = self.spotScale
            print "pickling...."
            pickle.dump(d, pfile, protocol=-1)

    def loadFromFile(self, filename):
        with open(filename, 'r') as pfile:
            d = pickle.load(pfile)

        self.spots = d['spots'] * 100
        self.wave = d['wave']
        self.fiber = d['fiber']
        self.coeffs = d['coeffs']
        self.spotScale = d['spotScale']
        
        print "WARNING:  forcing spotScale to 0.001"
        self.spotScale = 0.001
                
        print "WARNING: swapping x and y centers, cuz I think spectra disperse along columns"
        self.xcCoeffs = d['ycCoeffs']
        self.ycCoeffs = d['xcCoeffs']
