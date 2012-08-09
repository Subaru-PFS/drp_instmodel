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

def constructSplinesFromSpots(band, spotType='zemax'):
    psfFilepath = os.path.join(SplinedPsf.psfFileDir(band, spotType),
                               'psfSplines.pck')

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

    spots = s['spot'].astype('float64')
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

class SplinedPsf(psf.Psf):
    def __init__(self, band, detector, spotType='zemax'):
        """ Read in our persisted form. Optionally that from the zemax spots. """

        self.band = band
        self.wave = []
        self.fiber = []
        self.spots = []
        self.coeffs = []
        self.xcCoeffs = []
        self.ycCoeffs = []
        self.detector = detector

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
    
    def wavesForRows(self, fibers, rows=None, waveRange=None):
        """ Return our best estimate for the wavelength at the given row centers.

        Returns:
          rows   - the rows which we evaluated. [NROWS]
          waves  - for each fiber, the wavelengths at rows [NFIBERS, NROWS]
    
        """

        if waveRange == None:
            waveRange = self.wave.min(), self.wave.max()

        # Assume that the full spectrum fits on the detector.
        minY, maxY = self.ycCoeffs(fibers, waveRange)[0]
        doReorder = minY > maxY
        if doReorder:
            minY, maxY = maxY, minY
        assert minY > 0 and maxY < self.detector.config['ccdSize'][1]
        
        minRow = int(minY/self.detector.config['pixelScale']) + 1
        maxRow = int(maxY/self.detector.config['pixelScale'])

        if rows == None:
            rows = np.arange(minRow, maxRow+1)

        # Invert the spline into a row->wave map. 
        # Just use a linear interpolation based on the evaluation near the pixels.
        allWaves = np.linspace(waveRange[0], waveRange[1], maxRow-minRow)

        waves = []
        for f in fibers:
            allWaveRows = self.ycCoeffs([f], allWaves)[0] / self.detector.config['pixelScale']

            if doReorder:
                allWaveRows0 = allWaveRows[::-1]
                allWaves0 = allWaves[::-1]
            else:
                allWaveRows0 = allWaveRows
                allWaves0 = allWaves
                
            waveFunc = spInterp.interp1d(allWaveRows0, allWaves0, 'linear', bounds_error=False)
            fiberWaves = waveFunc(rows)
            waves.append(fiberWaves)

        return rows, waves
    
    def psfsAt(self, fibers, waves, everyNthPsf=1, usePsfs=None):
        """ Return a stack of PSFs, instantiated on a rectangular grid.

        Args:
           fibers, waves - two arrays of positions; we instatiate ourselves at the rectangular grid of these.

        Returns:
           - a 3D array of images [len(xs)*len(ys), imshape.x, imshape.y]
           - a list of the (fiberId, wavelength) positions the images are instantiated for.
           - a list of the (x,y) positions the images should be centered on.
        """

        waveSign = 1 if waves[-1] > waves[0] else -1
        interpWaves = waves[::waveSign*everyNthPsf]
        interpSlice = slice(None, None, waveSign*everyNthPsf)
            
        centers = [(x,y) for x in fibers for y in waves]
        if usePsfs != None:
            newImages = usePsfs
        else:
            newImages = np.zeros(shape=(len(fibers)*len(waves),
                                        self.imshape[0],
                                        self.imshape[1]))
            for ix in range(self.imshape[0]):
                if ix == 0:
                    print "fibers %s col %d" % (fibers, ix)
                for iy in range(self.imshape[1]):
                    newImages[interpSlice, iy, ix] = self.coeffs[iy, ix](fibers, interpWaves).flat

        if everyNthPsf > 1:
            interpSliceLen = newImages[interpSlice,:,:].shape[0]
            for i in range(1, everyNthPsf):
                if waves[0] > waves[-1]:
                    # a[-2::nth] and down
                    start = waveSign*(i+1)
                else:
                    # a[1::nth] and up
                    start = waveSign*i
                outSlice = slice(start, None, waveSign*everyNthPsf)
                if interpSliceLen > newImages[outSlice,:,:].shape[0]:
                    interpSlice = slice(None, waveSign*(interpSliceLen-1)*everyNthPsf, waveSign*everyNthPsf)
                
                newImages[outSlice,:,:] = newImages[interpSlice,:,:]
                    
        return newImages, centers, self.traceCenters(fibers, waves)

    def makeComb(self, ids, nth=100, hackScale=1000):
        """ Return a functor which returns True for every nth item in the initialization list. """
        class Comb(object):
            def __init__(self, ids, nth):
                self.waves = np.array([x[1] for x in ids])
                self.comb = self.waves[::nth].tolist()
                    
            def __call__(self, x):
                return hackScale * (x in self.comb)

        return Comb(ids, nth=nth)
    
    def fiberImages(self, fibers, spectra=None, outImg=None, waveRange=None, everyNthPsf=1):
        if outImg == None:
            outImg = self.detector.simBias().image

        if spectra == None:
            spectra = [None] * len(fibers)
        for i, fiber in enumerate(fibers):
            self.fiberImage(fiber, spectra[i], outImg=outImg, 
                            waveRange=waveRange, everyNthPsf=everyNthPsf)

        return outImg
    
    def fiberImage(self, fiber, spectrum=None, outImg=None, waveRange=None, everyNthPsf=1):
        """ Return an interpolated image of a fiber """

        if waveRange == None:
            waveRange = self.wave.min(), self.wave.max()
        minX, maxX = self.xcCoeffs([fiber], waveRange)[0]
        minY, maxY = self.ycCoeffs([fiber], waveRange)[0]

        minRow = minY/self.detector.config['pixelScale']
        maxRow = maxY/self.detector.config['pixelScale']
        minCol = minX/self.detector.config['pixelScale']

        # Generalize this... CPL
        if minRow > maxRow:
            minRow, maxRow = maxRow, minRow
            
        # Get the wavelengths for the fiber pixels.
        allWaveRows, allWaves = self.wavesForRows([fiber], waveRange=waveRange)
        allWaves = allWaves[0]
        
        # Get the oversampled PSFs and their locations at the fiber pixels.
        fiberPsfs, psfIds, centers = self.psfsAt([fiber], allWaves, everyNthPsf=everyNthPsf)
        xCenters, yCenters = [c[0] for c in centers]
        
        # With no input, use a comb spectrum
        if spectrum == None:
            spectrum = self.makeComb(psfIds, nth=20)

        traceWidth = int((xCenters.max() - xCenters.min())/self.detector.config['pixelScale'] + 0.5)
        traceHeight = int((yCenters.max() - yCenters.min())/self.detector.config['pixelScale'] + 0.5)

        psfPixelScale = int(self.detector.config['pixelScale'] / self.spotScale + 0.5)
        spotWidth = fiberPsfs[0].shape[-1] / psfPixelScale
        spotRadius = (spotWidth+1) / 2

        if outImg == None:
            outExp = self.detector.simBias(shape=(traceHeight + spotWidth, 
                                                  traceWidth + spotWidth))
            outImg = outExp.image
            outImgOffset = (xCenters.min(), yCenters.min())
        else:
            outImgOffset = (0,0)
            
        lasty = 0
        for i, row in enumerate(allWaveRows):
            rawPsf = fiberPsfs[i]
            fiberID, wave = psfIds[i]

            # in mm
            xc = xCenters[i]
            yc = yCenters[i]

            # pix offset
            xPixOffset = (xc-outImgOffset[0]) / self.detector.config['pixelScale'] + spotRadius
            yPixOffset = (yc-outImgOffset[1]) / self.detector.config['pixelScale'] + spotRadius

            # Keep the shift to the smallest fraction possible, or rather keep the integer steps 
            # exactly +/- 1.
            inty = round(yPixOffset)
            fracy = yPixOffset - inty

            intx = round(xPixOffset)
            fracx = xPixOffset - intx

            doDetails = False
            if row % 100 == 0:
                print("%5d %5d %6.1f yc: %3.4f %3.10f %d %d %3.10f " % (i, row, wave, yc, yPixOffset, lasty, inty, fracy)),
                doDetails = True
            lasty = inty
            
            # bin psf to ccd pixels, shift by fractional pixel only.
            psf = self.scalePsf(rawPsf, -fracx, -fracy, doDetails=doDetails)

            # Assumed to be a spline on the final resolution.
            flux = spectrum(wave)

            # When do we need to do this before the downsampling?
            spot = psf * flux

            self.placeSubimage(outImg, spot, intx, inty)

        return outImg, minRow, minCol

    def scalePsf(self, rawPsf, xc, yc, doRescaleFlux=False, doDetails=False):
        """ Given an oversampled image scale it to the detector grid after shifting it. """

        psfPixelScale = self.detector.config['pixelScale'] / self.spotScale
            
        psf = scipy.ndimage.affine_transform(rawPsf, 
                                             [psfPixelScale, psfPixelScale],
                                             (yc, xc),
                                             output_shape=np.array(rawPsf.shape)/psfPixelScale,
                                             order=1)

        if doDetails:
            print ("           shift=%d/%0.3f, %d/%0.3f flux=%0.2f,%0.2f" % 
                   (round(xc*psfPixelScale), xc, 
                    round(yc*psfPixelScale), yc,
                    np.sum(rawPsf), np.sum(psf) * (psfPixelScale * psfPixelScale)))

        if doRescaleFlux:
            psf *= (psfPixelScale*psfPixelScale)

        return psf

    def placeSubimage(self, img, subImg, xc, yc):
        parentx, childx = self.trimSpan((0, img.shape[1]-1),
                                        (0, subImg.shape[1]-1),
                                        xc)
        
        parenty, childy = self.trimSpan((0, img.shape[0]-1),
                                        (0, subImg.shape[0]-1),
                                        yc)

        # print("  %3.5f  parent %s  child %s" % (yc, parenty, childy))
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

    def frac(self, n):
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
            filename = os.path.join(SplinedPsf.psfFileDir(self.band, spotType),
                                    'psfSplines.pck')

        with open(filename, 'r') as pfile:
            d = pickle.load(pfile)

        self.spots = d['spots']
        self.wave = d['wave']
        self.fiber = d['fiber']
        self.coeffs = d['coeffs']
        self.spotScale = d['spotScale']

        # Make fiber -> (x,y) and wave->(x,y) splines. We have to do this here, because
        # the detector geometry matters.
        
        # Shift offsets to origin.
        xc = d['xc'] + self.detector.config['ccdSize'][1] * self.detector.config['pixelScale'] / 2
        yc = d['yc'] + self.detector.config['ccdSize'][0] * self.detector.config['pixelScale'] / 2

        # XXX - This fiber/wave indexing scheme is not safe in general
        xx = np.unique(self.fiber)
        yy = np.unique(self.wave)

        xcCoeffs = spInterp.RectBivariateSpline(xx, yy, xc.reshape(len(xx), len(yy)))
        ycCoeffs = spInterp.RectBivariateSpline(xx, yy, yc.reshape(len(xx), len(yy)))

        print "NOTE:  spotScale=%0.5f" % (self.spotScale)
                
        print "WARNING: swapping x and y centers, spectra actually disperse along columns"
        self.xc = d['yc']
        self.yc = d['xc']
        self.xcCoeffs = ycCoeffs
        self.ycCoeffs = xcCoeffs
