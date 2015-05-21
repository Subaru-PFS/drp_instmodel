#!/usr/bin/env python

import logging
import os
import time

from astropy.io import fits as pyfits
import numpy as np

import scipy.ndimage
import scipy.signal
import scipy.interpolate as spInterp

import psf
import pfs_tools
from pfs_tools import pydebug
import spotgames

class SplinedPsf(psf.Psf):
    def __init__(self, detector, spotType='jeg', spotID=None, logger=None):
        """ Create or read in our persisted form. By default use JEG's models. 

        Parameters
        ----------

        detector : a Detector object
           carries the properties of the detector (and camera).
        spotType : str, optional
           Whether to load 'jeg' or 'zemax' spot images.
        spotID : list/string/dict, optional
           If spotType is set, some dataset identifier. Opaque to us.
        """

        psf.Psf.__init__(self, detector, logger=logger)

        # The locations at which we have PSFs
        self.wave = []
        self.fiber = []
        self.spots = []

        self.coeffs = []
        self.xcCoeffs = []
        self.ycCoeffs = []
        self.spotID = spotID

        if spotType:
            self.loadFromSpots(spotType, spotID)

        
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
        dataRoot = os.environ.get('DRP_INSTDATA_DIR', '.')
        return os.path.join(dataRoot, 'data', 'spots', spotType, band)
        
    @staticmethod
    def psfSplinesFile(band, spotType):
        return os.path.join(SplinedPsf.psfFileDir(band, spotType),
                            'psfSplines.pck')

    @staticmethod
    def psfSpotsFile(band, spotType):
        return os.path.join(SplinedPsf.psfFileDir(band, spotType),
                            'spots.fits')
        
    @property 
    def imshape(self):
        return self.spots[0].shape

    def traceCenters(self, fibers, waves):
        """ Return the pixel centers for the given fibers and wavelengths """

        fibers = np.array(fibers)
        waves = np.array(waves)
        
        fiberIdx = np.argsort(fibers)
        waveIdx = np.argsort(waves)

        x = self.evalSpline(self.xcCoeffs, fibers[fiberIdx], waves[waveIdx])
        y = self.evalSpline(self.ycCoeffs, fibers[fiberIdx], waves[waveIdx])

        return (x[fiberIdx][:, waveIdx], 
                y[fiberIdx][:, waveIdx])
    
    def wavesForRows(self, fibers, rows=None, waveRange=None, pixelScale=None):
        """ Return our best estimate for the wavelength at the given row centers.

        Parameters
        ----------
        fibers : array_like
           the fibers we want solutions for.
        rows : array_like, optional
           the rows we want solutions for. If not set, all rows.
        waveRange : (low, high), optional
           limit the rows to those entirely within this range.  
        pixelScale : float, optional
           override our detector's mm/pixel scale
           
        Returns
        -------
        rows : float[NROWS]
            the rows which we evaluated. 
        waves : float[NFIBERS, NROWS]
            for each fiber, the wavelengths at rows
    
        """

        if pixelScale is None:
            pixelScale = self.detector.config['pixelScale']            
        
        if waveRange is None:
            waveRange = self.wave.min(), self.wave.max()

        # Assume that the full spectrum fits on the detector.
        minY, maxY = self.evalSpline(self.ycCoeffs, fibers, waveRange)[0]
        doReorder = minY > maxY
        if doReorder:
            minY, maxY = maxY, minY
        if minY < 0:
            print("one of wavelengths %s maps below the detector (%0.5f mm)" % (waveRange, minY))
        if maxY > self.detector.config['ccdSize'][1]:
            print("one of wavelengths %s maps above the detector (%0.5f mm)" % (waveRange, maxY))
        
        minRow = int((minY+1)/pixelScale)
        maxRow = int(maxY/pixelScale)

        if rows is None:
            rows = np.arange(minRow, maxRow+1, dtype='i4')

        # Invert the spline into a row->wave map. 
        # Just use a linear interpolation based on the evaluation near the pixels.
        allWaves = np.linspace(waveRange[0], waveRange[1], maxRow-minRow+1, dtype='f8')

        waves = []
        for f in fibers:
            allWaveRows = self.evalSpline(self.ycCoeffs, [f], allWaves)[0] / pixelScale

            if doReorder:
                allWaveRows0 = allWaveRows[::-1]
                allWaves0 = allWaves[::-1]
            else:
                allWaveRows0 = allWaveRows
                allWaves0 = allWaves
                
            waveFunc = spInterp.interp1d(allWaveRows0, allWaves0, 'linear', bounds_error=False)
            fiberWaves = waveFunc(rows)
            waves.append(fiberWaves.astype('f8'))

        return rows, waves
    
    def psfsAt(self, fibers, waves=None, usePsfs=None):
        """ Return a stack of PSFs, instantiated on a rectangular grid.

        Parameters
        ----------
        fibers : array_like, int
           the fiber IDs we want PSFs for
        waves : array_like, AA, optional
           the wavelengths we want PSFs as. If not set, uses the native wavelength grid.
        usePsfs : array_like , optional
           
        Returns
        -------
        psfImages : array
           a 3D array of images [len(fibers)*len(waves), imshape.x, imshape.y]
        psfIds : array
           the (fiberId, wavelength) positions the images are instantiated for.
        psfPos : array
           the (x,y) positions the PSF images should be centered on.

        """

        if waves is None:
            waves = np.unique(self.wave)
            
        waveSign = 1 if waves[-1] > waves[0] else -1
        interpWaves = waves[::waveSign]
        
        centers = [(x,y) for x in fibers for y in waves]
        if usePsfs is not None:
            newImages = usePsfs
        else:
            newImages = np.zeros(shape=(len(fibers)*len(interpWaves),
                                           self.imshape[0],
                                           self.imshape[1]), dtype='f4')
            for ix in range(self.imshape[0]):
                for iy in range(self.imshape[1]):
                    newImages[:, iy, ix] = self.evalSpline(self.coeffs[iy, ix], fibers, interpWaves).flat

        lo_w = newImages < 0
        if np.any(lo_w):
            minpix = newImages[lo_w].min()
            self.logger.warn("%d/%d low PSF pixels, min=%g" % (lo_w.sum(), newImages.size, minpix))
            newImages += minpix

        finalImages = newImages
        self.logger.info("psfsAt: fibers %s, for %d %s waves, returned %d %s unique psfs" % (fibers,
                                                                                             len(waves), waves[0].dtype,
                                                                                             len(interpWaves), finalImages.dtype))
        return finalImages, centers, self.traceCenters(fibers, waves)


    def fiberGeometry(self, fiber, waveRange=None):
        # Evaluate at image resolution
        pixelScale = self.detector.config['pixelScale']
        
        if waveRange is None:
            waveRange = self.wave.min(), self.wave.max()

        # Get the wavelengths for the fiber pixels.
        allPixelRows, allPixelWaves = self.wavesForRows([fiber], waveRange=waveRange, 
                                                        pixelScale=pixelScale)
        pixelWaves = allPixelWaves[0]

        centers = self.traceCenters([fiber], pixelWaves)
        
        return pixelWaves, centers

    def fiberImage(self, fiber, spectrum, outExp=None, waveRange=None, 
                   returnUnbinned=False,
                   shiftPsfs=True):
        """ Return an interpolated image of a fiber """

        # Evaluate at highest resolution
        pixelScale = self.spotScale
        
        if waveRange is None:
            waveRange = self.wave.min(), self.wave.max()

        # Get the wavelengths for the fiber pixels.
        self.logger.debug("waves range: %s", waveRange)
        allPixelRows, allPixelWaves = self.wavesForRows([fiber], waveRange=waveRange, 
                                                        pixelScale=pixelScale)
        self.logger.debug("allwaves: %s %d", allPixelWaves[0].dtype, len(allPixelWaves[0]))

        isLinelist = spectrum.__class__.__name__ == "CombSpectrum"
        waves, flux = spectrum.flux(allPixelWaves[0])

        self.logger.debug("waves: %s %d %d", waves.dtype, len(waves), waves.nbytes)
        self.logger.debug("flux:  %s %d %d", flux.dtype, len(flux), flux.nbytes)
        
        # Get the PSFs and their locations on the oversampled pixels.
        fiberPsfs, psfIds, centers = self.psfsAt([fiber], waves)

        xCenters, yCenters = [c[0] for c in centers]
        
        psfToSpotRatio = self.detector.config['pixelScale'] / pixelScale
        psfToSpotPixRatio = int(round(psfToSpotRatio))

        # pixels
        traceWidth = int((xCenters.max() - xCenters.min())/pixelScale) + 1
        traceHeight = int((yCenters.max() - yCenters.min())/pixelScale) + 1
        spotWidth = fiberPsfs[0].shape[-1]
        spotRad = spotWidth / 2
        self.logger.debug("spot size: %s %s %s %s" % (spotWidth, spotRad, psfToSpotRatio, psfToSpotPixRatio))
        self.logger.debug("trace    : %s %s" % (traceHeight, traceWidth))
        
        # We want out fiber image to be sized in units of outImgSpotPixelScale.
        fiHeight = traceHeight + spotWidth
        fiWidth = traceWidth + spotWidth
        chunkUp = np.array([fiHeight, fiWidth])
        chunkUp = psfToSpotPixRatio - chunkUp%psfToSpotPixRatio
        chunkUp[chunkUp == 10] = 0
        fiHeight += chunkUp[0]
        fiWidth += chunkUp[1]
        self.logger.debug("trace    : %s %s -> %s %s" % (traceHeight, traceWidth, fiHeight, fiWidth))

        # mm
        fiberImageOffset = np.asarray((yCenters.min(), xCenters.min()))

        # spot pixels; expand to include and fall on full image pixels
        fiberImagePixelOffset = (fiberImageOffset / pixelScale).astype('i4') - spotRad
        expandDown = fiberImagePixelOffset%psfToSpotPixRatio
        fiberImagePixelOffset -= expandDown
        
        # pixels
        outImgOffset = fiberImagePixelOffset / psfToSpotPixRatio
        self.logger.debug("fiber offset: pix=%s base=%s, mm=%s out=%s" % (fiberImagePixelOffset,
                                                                          fiberImageOffset/pixelScale, fiberImageOffset,
                                                                          outImgOffset))
        fiberImage = np.zeros((fiHeight, fiWidth), dtype='f4')

        if outExp is None:
            outExp = self.detector.makeExposure()

        # construct the oversampled fiber image
        geometry = np.zeros(len(waves), dtype=[('xc','f4'),('yc','f4'),
                                               ('intx','i4'),('inty','i4'),
                                               ('wavelength','f4'),('flux','f4')])
        for i in range(len(waves)):
            specWave = waves[i]
            specFlux = flux[i]

            if specFlux == 0.0:
                continue
            
            rawPsf = fiberPsfs[i]

            # in mm
            xc = xCenters[i]
            yc = yCenters[i]

            # pix offset
            xPixOffset = xc / pixelScale - fiberImagePixelOffset[1] - spotRad
            yPixOffset = yc / pixelScale - fiberImagePixelOffset[0] - spotRad

            # Keep the shift to the smallest fraction possible, or rather keep the integer steps 
            # exact.
            inty = int(yPixOffset)
            intx = int(xPixOffset)

            fracy = yPixOffset - inty
            fracx = xPixOffset - intx

            if shiftPsfs:
                shiftedPsf, kernels = spotgames.shiftSpot1d(rawPsf, fracx, fracy)
                spot = specFlux * shiftedPsf
            else:
                spot = specFlux * rawPsf
                
            if isLinelist or i % 1000 == 0 or i > len(waves)-2:
                self.logger.debug("%5d %6.1f (%3.3f, %3.3f) (%0.2f %0.2f) %0.2f %0.2f %0.2f" % (i, specWave, xc, yc, 
                                                                                                xPixOffset + spotRad,
                                                                                                yPixOffset + spotRad,
                                                                                                rawPsf.sum(), spot.sum(),
                                                                                                specFlux))
            self.placeSubimage(fiberImage, spot, (inty, intx))
            geometry[i] = (xc,yc,intx,inty,specWave,specFlux)


            if isLinelist:
                # mc1 = pfs_tools.centroid(fiberImage[inty-spotRad:inty+spotRad, intx-spotRad:intx+spotRad])
                mc2 = pfs_tools.centroid(spot)
                self.logger.debug("(%0.3f, %0.3f) (%0.3f, %0.3f)",
                                  xc / pixelScale, yc / pixelScale,
                                  mc2[0] + fiberImagePixelOffset[1] + intx,
                                  mc2[1] + fiberImagePixelOffset[0] + inty)
                
        # transfer flux from oversampled fiber image to final resolution output image
        resampledFiber = self.addOversampledImage(fiberImage, outExp, outImgOffset, psfToSpotPixRatio)

        if returnUnbinned:
            return outExp, fiberImagePixelOffset[0], fiberImagePixelOffset[1], geometry, fiberImage, resampledFiber
        else:
            return outExp, fiberImagePixelOffset[0], fiberImagePixelOffset[1], geometry

    def addOversampledImage(self, inImg, outExp, outOffset, outScale):
        """ Add the outScale-oversampled inImg to outImg at the given offset. 

        Must convolve with the pixel response of the oversampled pixel.
        """

        pixelKernel = np.ones((outScale, outScale), dtype='f4') / (outScale*outScale)
        pixImg = scipy.ndimage.convolve(inImg, pixelKernel, mode='constant', cval=0.0)
                                
        resampled = pfs_tools.rebin(pixImg, pixImg.shape[0]/outScale, pixImg.shape[1]/outScale)
        self.logger.debug("addOversampled kernel shape: %s, out type: %s",
                          pixelKernel.shape, resampled.dtype)
        parentIdx, childIdx = self.trimRect(outExp, resampled, outOffset)
        try:
            outExp.addFlux(resampled[childIdx], outSlice=parentIdx, addNoise=True)
        except Exception, e:
            self.logger.warn("failed to place child at %s): %s" % (outOffset, e))

        return resampled
    
    def placeSubimage(self, outImg, subImg, subOffset):
        parentIdx, childIdx = self.trimRect(outImg, subImg, subOffset)

        try:
            outImg[parentIdx] += subImg[childIdx]
            # exp.addFlux(subImg[childy, childx], outSlice=(parenty, parentx), addNoise=True)
        except Exception, e:
            self.logger.warn("failed to place child at %s: %s" % (subOffset, e))
    
    def trimRect(self, parent, child, childOffset=(0,0)):
        """ Given two images and an offset of the second, return the intersecting rectangles. 


        Arguments
        ---------
        parent, child : 2-d ndarray
        childOffset : offset of child on parent.

        Returns
        -------
        (parenty, parentx) - A pair of slice() objects.
        (childy, childx)   - A pair of slice() objects.

        """

        parentx, childx = self.trimSpan((0, parent.shape[1]-1),
                                        (0, child.shape[1]-1),
                                        childOffset[1])
        
        parenty, childy = self.trimSpan((0, parent.shape[0]-1),
                                        (0, child.shape[0]-1),
                                        childOffset[0])


        return (parenty, parentx), (childy, childx)

    def trimSpan(self, _parent, _child, childOffset=None):
        """ Return the final overlapping span as slices. We NEED to adopt some point/span/rect system. """ 

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
            return np.exp(-((x+dx)/dampfac)**2) * np.sin(np.pi*(x+dx)) / (np.pi * (x+dx))
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

    def loadFromSpots(self, spotType='jeg', spotIDs=None, spotArgs=None):
        """ Generate ourself from a semi-digested pfs_instdata spot file. 

        """
        
        self.logger.info("reading and interpolating %s PSF spots: %s..." % (spotType, spotIDs))
        if spotType == 'jeg':
            import jegSpots

            if spotArgs is None:
                spotArgs = dict()
            rawSpots, spotInfo = jegSpots.readSpotFile(spotIDs, verbose=True, **spotArgs)
            assert spotInfo['XPIX'] == spotInfo['YPIX']

            self.wave = rawSpots['wavelength']
            self.fiber = rawSpots['fiberIdx']
            self.spots = rawSpots['spot'].astype('float32')
            self.spotScale = spotInfo['XPIX']

        else:
            spotsFilepath = SplinedPsf.psfSpotsFile(self.detector.band, spotType, spotIDs)
            sf = pyfits.open(spotsFilepath, mode='readonly')
            rawSpots = sf[1].data

            self.spotScale = sf[0].header['pixscale']
            self.wave = rawSpots['wavelength']
            self.fiber = rawSpots['fiberIdx']
            self.spots = rawSpots['spot'].astype('float32')

        # Keep the raw info around, at least until production.
        self.rawSpots = rawSpots

        imshape = self.spots.shape[1:]

        # Make a spline for each pixel. The spots are centered on the output grid,
        # and we track the xc,yc offset separately.

        # XXX - Check that the input is properly sorted.
        xx = np.unique(self.fiber)
        yy = np.unique(self.wave)

        splineType = spInterp.RectBivariateSpline
        coeffs = np.zeros(imshape, dtype='O')
        for ix in range(imshape[0]):
            for iy in range(imshape[1]):
                coeffs[iy, ix] = self.buildSpline(splineType,
                                                  xx, yy, self.spots[:, iy, ix].reshape(len(xx), len(yy)))
        self.coeffs = coeffs

        # Shift offsets to origin.
        self.xc = rawSpots['spot_xc'] + self.detector.config['ccdSize'][1] * self.detector.config['pixelScale'] / 2
        self.yc = rawSpots['spot_yc'] + self.detector.config['ccdSize'][0] * self.detector.config['pixelScale'] / 2

        self.xcCoeffs = self.buildSpline(splineType,
                                         xx, yy, self.xc.reshape(len(xx), len(yy)))
        self.ycCoeffs = self.buildSpline(splineType,
                                         xx, yy, self.yc.reshape(len(xx), len(yy)))

    def buildSpline(self, splineType, x, y, z):
        if splineType is spInterp.RectBivariateSpline:
            return splineType(x, y, z)
        elif splineType is spInterp.RegularGridInterpolator:
            return splineType((x, y), z, method='linear')
        else:
            raise RuntimeError('evalSpline: unknown spline type: %s' % (splineType))

    def evalSpline(self, spline, x, y):
        if spline.__class__ is spInterp.RectBivariateSpline:
            return spline(x, y)
        elif spline.__class__ is spInterp.RegularGridInterpolator:
            x = np.asarray(x).flatten()
            y = np.asarray(y).flatten()
            xx,yy = np.meshgrid(x, y, indexing='ij')
    
            zz = np.array((xx,yy)).T.reshape(x.shape[0]*y.shape[0], 2)
    
            return np.atleast_2d(spline(zz))
        else:
            # print("warning: evalSpline: unknown spline type: %s'" % (spline))
            return spline(x, y)

    def spotGridImage(self):
        """ Return an array showing the rawspots in their given locations. """

        pass
        
    def setConstantSpot(self, spot=None, spotScale=None):
        """ Set ourself to use a single PSF over the whole focal plane. Uses the existing .fiber and .wave maps. """

        if spot is None:
            spotFiber = 0
            spotWave = sorted(self.wave)[len(self.wave)/2]
            spotId = np.where((self.fiber == spotFiber) & (self.wave == spotWave))[0]
            spot = self.spots[int(spotId)]

            self.logger.warn("set constant PSF to spot %s (fiber %d, wave %d)", spotId, spotFiber, spotWave)
            
        imshape = spot.shape
        self.spots = [spot]
        if spotScale is not None:
            self.spotScale = spotScale

        def freezeVal(val):
            return lambda x,y: np.ones((len(y), len(x))) * val
        
        coeffs = np.zeros(imshape, dtype='O')
        for ix in range(imshape[0]):
            for iy in range(imshape[1]):
                spotval = float(spot[iy,ix])
                coeffs[iy, ix] = freezeVal(spotval)
        self.coeffs = coeffs
        tx = ty = imshape[0]/2

    def setConstantX(self):
        """ Force spot x positions to be for the center spot on their trace. """
        
        xx = np.unique(self.fiber)
        yy = np.unique(self.wave)

        xc = self.xc.reshape(len(xx), len(yy))
        xc[:,:] = xc[:,len(yy)/2:len(yy)/2+1]
        self.xcCoeffs = self.buildSpline(self.xcCoeffs.__class__,
                                         xx, yy, xc)
        self.logger.warn("set constant X")
