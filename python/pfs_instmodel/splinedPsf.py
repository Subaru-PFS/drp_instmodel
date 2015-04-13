#!/usr/bin/env python

import os
import time
from astropy.io import fits as pyfits
import numpy

import scipy.ndimage
import scipy.signal
import scipy.interpolate as spInterp

import psf
import pfs_tools
from pfs_tools import pydebug
import spotgames

class SplinedPsf(psf.Psf):
    def __init__(self, detector, spotType='jeg', spotID=None):
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

        psf.Psf.__init__(self, detector)
        
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
                     len(numpy.unique(self.wave)), numpy.min(self.wave), numpy.max(self.wave),
                     len(numpy.unique(self.fiber))))
        else:
            return ("%s<0 spots>" %
                    (self.__class__.__name__))

    @staticmethod
    def psfFileDir(band, spotType):
        dataRoot = os.environ.get('PFS_INSTDATA_DIR', '.')
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

        fibers = numpy.array(fibers)
        waves = numpy.array(waves)
        
        fiberIdx = numpy.argsort(fibers)
        waveIdx = numpy.argsort(waves)

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
            rows = numpy.arange(minRow, maxRow+1)

        # Invert the spline into a row->wave map. 
        # Just use a linear interpolation based on the evaluation near the pixels.
        allWaves = numpy.linspace(waveRange[0], waveRange[1], maxRow-minRow+1)

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
            waves.append(fiberWaves)

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
            waves = numpy.unique(self.wave)
            
        waveSign = 1 if waves[-1] > waves[0] else -1
        interpWaves = waves[::waveSign]
        
        centers = [(x,y) for x in fibers for y in waves]
        if usePsfs is not None:
            newImages = usePsfs
        else:
            newImages = numpy.zeros(shape=(len(fibers)*len(interpWaves),
                                           self.imshape[0],
                                           self.imshape[1]))
            for ix in range(self.imshape[0]):
                for iy in range(self.imshape[1]):
                    newImages[:, iy, ix] = self.evalSpline(self.coeffs[iy, ix], fibers, interpWaves).flat

        lo_w = newImages < 0
        if numpy.any(lo_w):
            minpix = newImages[lo_w].min()
            print("%d/%d low PSF pixels, min=%g" % (lo_w.sum(), newImages.size, minpix))
            newImages += minpix

        finalImages = newImages
        print "psfsAt: fibers %s, for %d waves, returned %d unique psfs" % (fibers, len(waves),
                                                                            len(interpWaves))
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
        minX, maxX = self.evalSpline(self.xcCoeffs, [fiber], waveRange)[0]
        minY, maxY = self.evalSpline(self.ycCoeffs, [fiber], waveRange)[0]

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
        fiberPsfs, psfIds, centers = self.psfsAt([fiber], pixelWaves)

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
        fiberImage = numpy.zeros((fiHeight, fiWidth), dtype='f4')

        # mm
        fiberImageOffset = numpy.asarray((yCenters.min(), xCenters.min()))

        # pixels
        outImgOffset = fiberImageOffset / self.detector.config['pixelScale']
        print("fiber offset: pix=%s mm=%s" % (outImgOffset, fiberImageOffset))

        # Adjust offset by fractional detector pixel
        fiberImageOffset -= (outImgOffset - outImgOffset.astype('i4')) * self.detector.config['pixelScale']
        outImgOffset = outImgOffset.astype('i4')
        print("net fiber offset: pix=%s mm=%s" % (outImgOffset, fiberImageOffset))
      
        if outExp is None:
            outExp = self.detector.makeExposure()

        # construct the oversampled fiber image
        geometry = numpy.zeros(len(pixelWaves), dtype=[('xc','f4'),('yc','f4'),
                                                       ('intx','i4'),('inty','i4'),
                                                       ('wavelength','f4'),('flux','f4')])
        for i in range(len(pixelWaves)):
            specWave = pixelWaves[i]
            specFlux = pixelFlux[i]

            if specFlux == 0.0:
                continue
            
            rawPsf = fiberPsfs[i]

            # in mm
            xc = xCenters[i]
            yc = yCenters[i]
            xoffset = xc - fiberImageOffset[1]
            yoffset = yc - fiberImageOffset[0]

            # pix offset
            xPixOffset = xoffset / pixelScale
            yPixOffset = yoffset / pixelScale

            # Keep the shift to the smallest fraction possible, or rather keep the integer steps 
            # exact.
            inty = int(yPixOffset)
            intx = int(xPixOffset)

            fracy = yPixOffset - inty
            fracx = xPixOffset - intx

            if shiftPsfs:
                shiftedPsf, kernels = spotgames.shiftSpot1d(rawPsf, fracx, 0)
                spot = specFlux * shiftedPsf
            else:
                spot = specFlux * rawPsf
                
            if i % 1000 == 0 or i > len(pixelWaves)-2:
                print("%5d %6.1f (%3.3f, %3.3f) %0.2f %0.2f %0.2f shift=%s" % (i, specWave, xc, yc, 
                                                                               rawPsf.sum(), spot.sum(), specFlux,
                                                                               shiftPsfs))

            # This is unaccountably expensive.
            #if spot.sum() < 1:
            #    continue
            
            self.placeSubimage(fiberImage, spot, (inty, intx))
            geometry[i] = (xc,yc,intx,inty,specWave,specFlux)

            # bin psf to ccd pixels, shift by fractional pixel only.
            #psf = self.scalePsf(rawPsf, -fracx, -fracy, doDetails=doDetails)
            #self.placeSubimage(outImg, spot, (inty, intx))

        # transfer flux from oversampled fiber image to final resolution output image
        resampledFiber = self.addOversampledImage(fiberImage, outExp, outImgOffset, outImgSpotPixelScale)

        if returnUnbinned:
            return outExp, minRow, minCol, geometry, fiberImage, resampledFiber
        else:
            return outExp, minRow, minCol, geometry

    def addOversampledImage(self, inImg, outExp, outOffset, outScale):
        """ Add the outScale-oversampled inImg to outImg at the given offset. """

        resampled = pfs_tools.rebin(inImg, inImg.shape[0]/outScale, inImg.shape[1]/outScale)

        parentIdx, childIdx = self.trimRect(outExp, resampled, outOffset)
        try:
            outExp.addFlux(resampled[childIdx], outSlice=parentIdx, addNoise=True)
        except Exception, e:
            print "failed to place child at %s): %s" % (outOffset, e)

        return resampled
    
    def placeSubimage(self, outImg, subImg, subOffset):
        parentIdx, childIdx = self.trimRect(outImg, subImg, subOffset)

        try:
            outImg[parentIdx] += subImg[childIdx]
            # exp.addFlux(subImg[childy, childx], outSlice=(parenty, parentx), addNoise=True)
        except Exception, e:
            print "failed to place child at %s: %s" % (subOffset, e)
    
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

        return numpy.modf(n)[0]
    
    def _sincfunc(self, x, dx):
        dampfac = 3.25
        if dx != 0.0:
            return numpy.exp(-((x+dx)/dampfac)**2) * numpy.sin(numpy.pi*(x+dx)) / (numpy.pi * (x+dx))
        else:
            xx = numpy.zeros(len(x))
            xx[len(x)/2] = 1.0
            return xx

    def _sincshift(self, image, dx, dy):
        """ UNUSED & UNTESTED from boss -- for comparison only, etc. etc. """
        
        sincrad = 10
        s = numpy.arange(-sincrad, sincrad+1)
        sincx = self._sincfunc(s, dx)

        #- If we're shifting just in x, do faster 1D convolution with wraparound                                
	    #- WARNING: can introduce edge effects if PSF isn't nearly 0 at edge                                    
        if abs(dy) < 1e-6:
            newimage = scipy.signal.convolve(image.ravel(), sincx, mode='same')
            return newimage.reshape(image.shape)

        sincy = self._sincfunc(s, dy)
        kernel = numpy.outer(sincy, sincx)
        newimage = scipy.signal.convolve2d(image, kernel, mode='same')
        return newimage

    def loadFromSpots(self, spotType='jeg', spotIDs=None, spotArgs=None):
        """ Generate ourself from a semi-digested pfs_instdata spot file. 

        """
        
        print "reading and interpolating %s PSF spots: %s..." % (spotType, spotIDs)
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
        xx = numpy.unique(self.fiber)
        yy = numpy.unique(self.wave)

        splineType = spInterp.RectBivariateSpline
        coeffs = numpy.zeros(imshape, dtype='O')
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
            x = numpy.asarray(x).flatten()
            y = numpy.asarray(y).flatten()
            xx,yy = numpy.meshgrid(x, y, indexing='ij')
    
            zz = numpy.array((xx,yy)).T.reshape(x.shape[0]*y.shape[0], 2)
    
            return numpy.atleast_2d(spline(zz))
        else:
            raise RuntimeError('evalSpline: unknown spline type: %s' % (spline))
        

    def spotGridImage(self):
        """ Return an array showing the rawspots in their given locations. """

        pass
        
    def setConstantSpot(self, spot, spotScale):
        """ Set ourself to use a single PSF over the whole focal plane. Uses the existing .fiber and .wave maps. """

        imshape = spot.shape
        self.spots = [spot]
        self.spotScale = spotScale

        def freezeVal(val):
            return lambda x,y: numpy.ones((len(y), len(x))) * val
        coeffs = numpy.zeros(imshape, dtype='O')
        for ix in range(imshape[0]):
            for iy in range(imshape[1]):
                spotval = float(spot[iy,ix])
                coeffs[iy, ix] = freezeVal(spotval)
        self.coeffs = coeffs
        tx = ty = 128
        print "%d,%d: %s)" % (tx, ty, coeffs[ty, tx]([0],[0]))
        


