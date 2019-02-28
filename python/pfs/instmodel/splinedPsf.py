#!/usr/bin/env python

from importlib import reload
import logging
import os
import time

from astropy.io import fits as pyfits
import numpy as np

import pickle
import gzip

import scipy.ndimage
import scipy.signal
import scipy.interpolate as spInterp
import scipy.ndimage
import scipy.ndimage.interpolation

from . import spotgames
from .spectrum import LineSpectrum, ArcSpectrum

from . import psf
reload(psf)


class SpotCoeffs(object):
    """ Trivial wrapper around pixel and position spline coefficients, so that we can 
        treat per-fiber (1-d) splines the same as sparse (2-d) splines. 
    """
    def __init__(self, pixelCoeffs, xcCoeffs, ycCoeffs, spots=None):
        self.pixelCoeffs = pixelCoeffs
        self.xcCoeffs = xcCoeffs
        self.ycCoeffs = ycCoeffs
        self.spots = spots

class SpotCache(object):
    def __init__(self, dir, logLevel=logging.DEBUG):
        self.dir = dir
        self.logger = logging.getLogger('spotCache')
        self.lastFiber = None, None

    def getInfoPath(self):
        return os.path.join(self.dir, 'info.pkl.gz')
    
    def getFiberPath(self, fiber):
        return os.path.join(self.dir, "fiber-%03d.pkl.gz" % (fiber))

    def create(self, info, force=False):
        if self.exists() and not force:
            raise RuntimeError("spot cache %s already exists" % (self.dir))

        if not os.path.isdir(self.dir):
            os.mkdir(self.dir)
        with gzip.open(self.getInfoPath(), 'wb+') as f:
            pickle.dump(info, f, -1)

    def getInfo(self):
        with gzip.open(self.getInfoPath(), 'rb') as f:
            info = pickle.load(f)
        return info

    def exists(self):
        return os.path.isfile(self.getInfoPath()) and os.path.isfile(self.getFiberPath(650))

    def __getitem__(self, fiber):
        self.logger.debug("fetching %s", self.getFiberPath(fiber))
        t0 = time.time()
        if self.lastFiber[0] == fiber:
            data = self.lastFiber[1]
        else:
            with gzip.open(self.getFiberPath(fiber), 'rb') as f:
                data = pickle.load(f)
        t1 = time.time()
        self.logger.debug("fetched in %0.2fs", (t1-t0))
        self.lastFiber = (fiber, data)

        return data

    def __setitem__(self, fiber, coeffs):
        with gzip.open(self.getFiberPath(fiber), 'wb+') as f:
            pickle.dump(coeffs, f, -1)

class SplinedPsf(psf.Psf):
    def __init__(self, detector, spotType='jeg', spotID=None,
                 logger=None, logLevel=logging.INFO,
                 slitOffset=(0.0, 0.0), everyNth=20,
                 doTrimSpots=True, doRebin=False):
        """ Create or read in our persisted form. By default use JEG's models. 

        Parameters
        ----------

        detector : a Detector object
           carries the properties of the detector (and camera).
        spotType : str, optional
           Whether to load 'jeg' or 'zemax' spot images.
        spotID : list/string/dict, optional
           If spotType is set, some dataset identifier. Opaque to us.
        doRebin : bool/integer
           How much to unbin raw spots by.
        slitOffset : pair
           X, wavelength offset of the slit. Used for dithered flats.
        """

        psf.Psf.__init__(self, detector, logger=logger)

        self.logger.setLevel(logLevel)
        
        # The locations at which we have PSFs
        self.wave = []
        self.fiber = []
        self.spots = []

        self.spotID = spotID
        self.slitOffset = slitOffset
        self.everyNth = everyNth
        self.perFiberCoeffs = None

        self.spotID['band'] = self.detector.armName
        if spotType:
            self.loadFromSpots(spotType, spotID, spotArgs=dict(doTrimSpots=doTrimSpots, doRebin=doRebin))

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
    def psfFileDir(armName, spotType):
        dataRoot = os.environ.get('DRP_INSTDATA_DIR', '.')
        return os.path.join(dataRoot, 'data', 'spots', spotType, armName)
        
    @staticmethod
    def psfSplinesFile(armName, spotType):
        return os.path.join(SplinedPsf.psfFileDir(armName, spotType),
                            'psfSplines.pck')

    @staticmethod
    def psfSpotsFile(armName, spotType):
        return os.path.join(SplinedPsf.psfFileDir(armName, spotType),
                            'spots.fits')
        
    @property 
    def imshape(self):
        return self.spotsInfo['YPIX'], self.spotsInfo['XPIX']

    def getCoeffs(self, fibers):
        if self.perFiberCoeffs is None:
            return self.allFiberCoeffs
        else:
            try:
                if len(fibers) != 1:
                    raise ValueError('only one fiber can be requested for per-fiber models')
                return self.perFiberCoeffs[fibers[0]]
            except TypeError:
                return self.perFiberCoeffs[fibers]

    def getCards(self):
        cards = [('HIERARCH sim.slit.xoffset', self.slitOffset[0], 'slit fiber offset, mm'),
                 ('HIERARCH sim.slit.yoffset', self.slitOffset[1], 'slit wavelength offset, mm')]

        return cards

    def _focalPlaneXToDetectorX(self, x, correctToLL=False):
        """ Given focal plane positions, return detector positions.

        Args
        ----
        x : ndarray
            mm from center
        correctToLL : bool
            whether to shift to left edge.

        Returns
        -------
        x : adjusted coordinates.

        """
        halfGap = self.detector.config['interCcdGap'] / 2
        halfPoint = 0.0

        # Chip gap
        neg_w = x < halfPoint
        pos_w = x > halfPoint
        x[neg_w] += halfGap
        x[pos_w] -= halfGap

        if correctToLL:
            x += self.detector.xcMmOffset

        return x

    def _focalPlaneYToDetectorY(self, y, correctToLL=False):
        """ Given focal plane positions, return detector positions.

        Args
        ----
        y : ndarray
            mm from center
        correctToLL : bool
            whether to shift to bottom edge.

        Returns
        -------
        y : adjusted coordinates.

        """

        if correctToLL:
            y += self.detector.ycMmOffset

        return y

    def traceCenters(self, fibers, waves):
        """ Return the pixel centers in mm for the given fibers and wavelengths """

        fibers = np.array(fibers)
        waves = np.array(waves)

        fiberIdx = np.argsort(fibers)
        waveIdx = np.argsort(waves)

        self.logger.debug("tracing fiber %s", fibers[fiberIdx])

        coeffs = self.getCoeffs(fibers)
        x = self.evalSpline(coeffs.xcCoeffs, fibers[fiberIdx], waves[waveIdx])
        y = self.evalSpline(coeffs.ycCoeffs, fibers[fiberIdx], waves[waveIdx])

        x = self._focalPlaneXToDetectorX(x, correctToLL=True)
        y = self._focalPlaneYToDetectorY(y, correctToLL=True)

        return (x[waveIdx], y[waveIdx])

    def wavesForRow(self, fiber, waveRange=None, pixelScale=None):
        """ Return our best estimate for the wavelength at the given row centers.

        Parameters
        ----------
        fiber : int
           the fiber we want a solution for.
        waveRange : (low, high), optional
           limit the rows to those entirely within this range.  
        pixelScale : float, optional
           override our detector's mm/pixel scale

        Returns
        -------
        rows : float[NROWS]
            the rows which we evaluated. 
        waves : float[NROWS]
            the wavelengths at the rows
        """

        if pixelScale is None:
            pixelScale = self.detector.config['pixelScale']

        if waveRange is None:
            waveRange = self.wave.min(), self.wave.max()

        coeffs = self.getCoeffs(fiber)

        # Assume that the full spectrum fits on the detector.
        minY, maxY = self.evalSpline(coeffs.ycCoeffs, fiber, waveRange)

        doReorder = minY > maxY
        if doReorder:
            minY, maxY = maxY, minY
        if minY + self.detector.ycMmOffset < 0:
            self.logger.info("one of wavelengths %s maps below the detector (%0.5f mm)" % (waveRange,
                                                                                           minY + self.detector.ycMmOffset))
        if maxY > self.detector.config['ccdSize'][1]:
            self.logger.info("one of wavelengths %s maps above the detector (%0.5f mm)" % (waveRange,
                                                                                           maxY + self.detector.ycMmOffset))

        minRow = int((minY + self.detector.ycMmOffset)/pixelScale)
        maxRow = int((maxY + self.detector.ycMmOffset)/pixelScale)
        rows = np.arange(minRow, maxRow+1, dtype='i4')

        # Invert the spline into a row->wave map. 
        # Just use a linear interpolation based on the evaluation near the pixels.
        rowWaves = np.linspace(waveRange[0], waveRange[1], maxRow-minRow+1, dtype='f8')

        return rows, rowWaves

    def psfsAt(self, fibers, waves=None, usePsfs=None, everyNth=None, doFill=True):
        """ Return a stack of PSFs, instantiated on a rectangular grid.

        Parameters
        ----------
        fibers : array_like, int
           the fiber IDs we want PSFs for
        waves : array_like, AA, optional
           the wavelengths we want PSFs as. If not set, uses the native wavelength grid.
        everyNth :integer
           How many psf instances to reuse
        usePsfs : array_like , optional
           existing PSF images to reuse.

        Returns
        -------
        psfImages : array
           a 3D array of images [len(fibers)*len(waves), imshape.x, imshape.y]
        psfIds : array
           the (fiberId, wavelength) positions the images are instantiated for.
        psfPos : array
           the (x,y) positions the PSF images should be centered on.

        """

        if isinstance(fibers, int):
            fibers = [fibers]
        coeffs = self.getCoeffs(fibers)

        if everyNth is None:
            everyNth = self.everyNth
        self.logger.debug('everyNth: %s', everyNth)

        if waves is None:
            waves = np.unique(self.wave)

        # Why did I think I needed to do this? Check. XXXCPL
        waveSign = 1 if waves[-1] > waves[0] else -1

        allWaves = waves[::waveSign]
        psfWaves = allWaves[::everyNth]

        centers = [(x,y) for x in fibers for y in waves]
        if usePsfs is not None:
            newImages = usePsfs
            t0 = t1 = t2 = 0.0
        else:
            t0 = time.time()
            newImages = np.zeros(shape=(len(fibers)*len(allWaves),
                                        self.imshape[0],
                                        self.imshape[1]), dtype='f8')
            self.logger.debug("psfsAt: fibers %s for %d/%d waves, in %s %s array" % (fibers, len(psfWaves),
                                                                                     len(allWaves),
                                                                                     newImages.shape,
                                                                                     newImages.dtype))
            for ix in range(self.imshape[0]):
                for iy in range(self.imshape[1]):
                    newImages[::everyNth, iy, ix] = self.evalSpline(coeffs.pixelCoeffs[iy, ix],
                                                                    fibers, psfWaves).flat
                    if everyNth > 1:
                        newImages[-1:, iy, ix] = self.evalSpline(coeffs.pixelCoeffs[iy, ix],
                                                                 fibers, allWaves[-1]).flat
            t1 = time.time()
            if everyNth > 1 and doFill:
                nSpots = len(newImages)
                nBlocks = (nSpots + everyNth - 1)//everyNth
                for block_i in range(nBlocks):
                    i0 = block_i * everyNth
                    i1 = i0 + everyNth
                    if i1 > nSpots-1:
                        i1 = nSpots-1
                    if i0 == i1:
                        dimg = 0
                    else:
                        dimg = (newImages[i1] - newImages[i0]) / (i1-i0)
                    # print("%d/%d %d: %d %d" % (block_i, nBlocks, nSpots, i0, i1))
                    for i in range(i0+1, i1):
                        newImages[i,:,:] = newImages[i-1,:,:] + dimg
            t2 = time.time()

        lo_w = newImages < -0.01
        if np.any(lo_w):
            minpix = newImages[lo_w].min()
            self.logger.warn("%d/%d low PSF pixels, min=%g" % (lo_w.sum(), newImages.size, minpix))
            newImages += minpix

        finalImages = newImages
        self.logger.debug("psfsAt: fibers %s, for %d %s waves, returned %d unique psfs. r=%g,%g" % (fibers,
                                                                                                    len(waves),
                                                                                                    waves[0].dtype,
                                                                                                    len(psfWaves),
                                                                                                    t1-t0, t2-t1))
        return finalImages, centers, self.traceCenters(fibers, waves)


    def fiberGeometry(self, fiber, waveRange=None):
        """ Return the wavelengths and positions of the pixels of the given fiber.
        """
        
        # Evaluate at image resolution
        pixelScale = self.detector.config['pixelScale']
        
        if waveRange is None:
            waveRange = self.wave.min(), self.wave.max()

        # Get the wavelengths for the fiber pixels.
        pixelRows, pixelWaves = self.wavesForRow(fiber, waveRange=waveRange,
                                                 pixelScale=pixelScale)

        centers = self.traceCenters([fiber], pixelWaves)

        return pixelWaves, centers

    def genFiberImage(self, fiber, spectrum, everyNth=1):
        """ Return an oversampled fiber image.

        Args:
         fiber : integer
           The holeId.
         spectrum: Spectrum
           The lambda -> flux map
         everyNth : integer
           How tightly to fully interpolate spots. 

        Returns:
          fiberImage : ndarray
            the oversampled fiberimage, trimmed to the full-pixel grid
          imageOffset: (x, y)
            the offsets into the full detector image.
          geometry : tuple
            stuff

        """

        pass

    def placeFiberImage(self):
        pass

    def addFibers(self, fibers, spectra, everyNth=1):
        pass
    
    def fiberImage(self, fiber, spectrum, waveRange=None,
                   shiftPsfs=True, everyNth=None, applyThroughput=True):
        """ Return an interpolated image of a fiber """

        # Evaluate at highest resolution
        pixelScale = self.spotScale
        
        if waveRange is None:
            waveRange = self.wave.min(), self.wave.max()

        # Get the wavelengths for the fiber pixels.
        self.logger.debug("waves range: %s", waveRange)

        pixelRows, pixelWaves = self.wavesForRow(fiber, waveRange=waveRange,
                                                 pixelScale=pixelScale)
        lower = np.empty_like(pixelWaves)
        upper = np.empty_like(pixelWaves)
        middle = 0.5*(pixelWaves[:-1] + pixelWaves[1:])
        # The bounds for most are easy
        lower[1:] = middle
        upper[:-1] = middle
        # The ends have to be handled with a bit more care
        lower[0] = pixelWaves[0] - 0.5*(pixelWaves[1] - pixelWaves[0])
        upper[-1] = pixelWaves[-1] + 0.5*(pixelWaves[-1] - pixelWaves[-2])

        flux = spectrum.integrate(lower, upper)  # W/m^2
        # Convert flux to photons
        # E = n.h.nu = F.T.A
        # n = F.T.A.lambda/h/c
        #   = (F / W/m^2) . (T / sec) . (A / m^2) . 10^-9 . (lambda/nm) / 6.63e-34 J.s / 3e8 m/s
        # aperture is 8.2m (effective; https://subarutelescope.org/Introduction/telescope.html)
        # Integration time gets factored in elsewhere.
        area = np.pi*(0.5*8.2)**2
        planck = 6.63e-34  # J.s
        speedOfLight = 3.0e8*1.0e9  # nm/s
        flux *= area*pixelWaves/planck/speedOfLight

        self.logger.info("Spectral min/max: %f %f", flux.min(), flux.max())

        waves = pixelWaves  # Change name for historical reasons

        # This is temporary: the NIST line strengths were adjusted to match the observed
        # intensities. We need to correct the line list instead of doing that here. XXXCPL
        if applyThroughput:
            self.logger.debug('applying throughput')
            tp = self.detector.throughput(waves)
            flux *= tp

        self.logger.debug("waves: %s %d %d", waves.dtype, len(waves), waves.nbytes)
        self.logger.debug("flux:  %s %d %d", flux.dtype, len(flux), flux.nbytes)

        # Get the PSFs and their locations on the oversampled pixels.
        fiberPsfs, psfIds, centers = self.psfsAt([fiber], waves, everyNth=everyNth)

        xCenters, yCenters = centers

        psfToSpotRatio = self.detector.config['pixelScale'] / pixelScale
        psfToSpotPixRatio = int(round(psfToSpotRatio))

        # pixels
        traceWidth = int((xCenters.max() - xCenters.min())/pixelScale) + 1
        traceHeight = int((yCenters.max() - yCenters.min())/pixelScale) + 1
        spotWidth = fiberPsfs[0].shape[-1]
        spotRad = spotWidth // 2
        self.logger.debug("spot size: %s %s %s %s" % (spotWidth, spotRad, psfToSpotRatio, psfToSpotPixRatio))
        self.logger.debug("trace    : %s %s" % (traceHeight, traceWidth))
        
        # We want out fiber image to be sized in units of outImgSpotPixelScale.
        fiHeight = traceHeight + spotWidth
        fiWidth = traceWidth + spotWidth
        chunkUp = np.array([fiHeight, fiWidth])
        chunkUp = psfToSpotPixRatio - chunkUp%psfToSpotPixRatio
        chunkUp[chunkUp == psfToSpotPixRatio] = 0
        fiHeight += chunkUp[0] + psfToSpotPixRatio
        fiWidth += chunkUp[1] + psfToSpotPixRatio
        self.logger.info("fiber %-3d   : %s %s %s -> %s %s" % (fiber, spectrum,
                                                               traceHeight, traceWidth,
                                                               fiHeight, fiWidth))

        # mm
        fiberImageOffset = np.asarray((yCenters.min(), xCenters.min()))

        # spot pixels; expand to include and fall on full image pixels
        fiberImagePixelOffset = (fiberImageOffset / pixelScale).astype('i4') - spotRad
        expandDown = fiberImagePixelOffset%psfToSpotPixRatio
        fiberImagePixelOffset -= expandDown
        
        # pixels
        outImgOffset = fiberImagePixelOffset // psfToSpotPixRatio
        self.logger.debug("fiber offset: pix=%s base=%s, mm=%s out=%s" % (fiberImagePixelOffset,
                                                                          fiberImageOffset/pixelScale,
                                                                          fiberImageOffset,
                                                                          outImgOffset))
        fiberImage = np.zeros((fiHeight, fiWidth), dtype='f4')

        # construct the oversampled fiber image
        geometry = np.zeros(len(waves), dtype=[('xc','f8'),('yc','f8'),
                                               ('dxc','f8'),('dyc','f8'),
                                               ('intx','i4'),('inty','i4'),
                                               ('wavelength','f8'),('flux','f4')])

        ymin = 1
        ymax = -1
        
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
            xPixOffset = xc/pixelScale - fiberImagePixelOffset[1] - spotRad + psfToSpotPixRatio/2
            yPixOffset = yc/pixelScale - fiberImagePixelOffset[0] - spotRad + psfToSpotPixRatio/2

            # XXX Shift the slit
            # This assumes everything moves linearly, which isn't right, but we don't currently
            # have an alternative.
            xPixOffset += self.slitOffset[0]*psfToSpotPixRatio
            yPixOffset += self.slitOffset[1]*psfToSpotPixRatio

            # Keep the shift to the smallest fraction possible, or rather keep the integer steps 
            # exact.
            inty = int(np.rint(yPixOffset))
            intx = int(np.rint(xPixOffset))

            fracy = yPixOffset - inty
            fracx = xPixOffset - intx
            if fracy < ymin:
                ymin = fracy
            elif fracy > ymax:
                ymax = fracy

            if shiftPsfs:
                if True:
                    shiftedPsf, kernels = spotgames.shiftSpot(rawPsf, fracx, fracy)
                elif False:
                    shiftedPsf, kernels = spotgames.shiftSpot(rawPsf, fracx, fracy, method='1d',
                                                              kargs=dict(n=spotRad))
                else:
                    shiftedPsf = scipy.ndimage.interpolation.shift(rawPsf,(fracy,fracx))

                _c0x, _c0y = spotgames.centroid(rawPsf)
                dxc = _c0x - spotRad
                dyc = _c0y - spotRad
                spot = specFlux * shiftedPsf
            else:
                spot = specFlux * rawPsf

            if (i % 1000 == 0 or i > len(waves)-2):
                self.logger.debug("%5d %6.1f (%3.3f, %3.3f) (%0.2f %0.2f) %0.2f %0.2f %0.2f",
                                  i, specWave,
                                  xc/pixelScale, yc/pixelScale, 
                                  xPixOffset, yPixOffset,
                                  rawPsf.sum(), spot.sum(),
                                  specFlux)

            self.placeSubimage(fiberImage, spot, (inty, intx))
            geometry[i] = (xc,yc,dxc,dyc,intx,inty,specWave,specFlux)

        return fiberImage, outImgOffset, psfToSpotPixRatio, geometry

    def addOversampledImage(self, inImg, outExp, outOffset, outScale, skyImage, skyOffset):
        """ Add the outScale-oversampled inImg to outImg at the given offset. 
        """

        resampled = spotgames.rebinBy(inImg, outScale)
        parentIdx, childIdx = self.trimRect(outExp, resampled, outOffset)

        try:
            outExp.addFlux(resampled[childIdx], outSlice=parentIdx)
        except Exception as e:
            self.logger.warn("failed to place child at %s, slices=%s,%s): %s" %
                             (outOffset, childIdx, parentIdx, e))

        if skyImage is not None and skyOffset is not None:
            skyResampled = spotgames.rebinBy(skyImage, outScale)
            skyParent, skyChild = self.trimRect(outExp, skyResampled, skyOffset)
            try:
                outExp.addSky(skyResampled[skyChild], outSlice=skyParent)
            except Exception as e:
                self.logger.warn("failed to place sky child at %s, slices=%s,%s): %s" %
                                 (skyOffset, skyChild, skyParent, e))

        return resampled
    
    def placeSubimage(self, outImg, subImg, subOffset):
        parentIdx, childIdx = self.trimRect(outImg, subImg, subOffset)

        try:
            outImg[parentIdx] += subImg[childIdx]
        except Exception as e:
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

        assert parent[1]-parent[0] == child[1]-child[0]

        return (slice(parent[0], parent[1]+1),
                slice(child[0], child[1]+1))

    def _spotCacheDir(self, spotIds=None):
        from . import jegSpots
        reload(jegSpots)

        rawSpotPath, _ = jegSpots.resolveSpotPathSpec(spotIds)
        
        spotDir = os.path.dirname(rawSpotPath)
        cacheDir = os.path.join(spotDir, '_cache')

        return cacheDir
    
    def primeSpotCache(self, spotIDs, rawSpotPath):
        cacheDir = self._spotCacheDir(rawSpotPath)
        cache = SpotCache(cacheDir)
        if cache.exists():
            self.perFiberCoeffs = cache
            spotsInfo = cache.getInfo()
            self.wave = spotsInfo['wave']
            self.fiber = spotsInfo['fiber']
            self.spotScale = spotsInfo['spotScale']
            self.spotsInfo = spotsInfo['spotsInfo']
            return True
        return False
    
    def createSpotCache(self, spotIDs, rawSpotPath, doRebin=None, doTrimSpots=False, spotArgs=None):
        from . import jegSpots
        reload(jegSpots)

        if spotArgs is None:
            spotArgs = {}
        rawSpots, spotsInfo = jegSpots.readSpotFile(rawSpotPath, verbose=True, **spotArgs)
        assert spotsInfo['XPIX'] == spotsInfo['YPIX']
        self.logger.info('loaded spot data from %s' % rawSpotPath)

        self.wave = rawSpots['wavelength']
        self.fiber = rawSpots['fiberIdx']
        self.spotScale = spotsInfo['PIXSIZE']
        self.spotsInfo = spotsInfo

        imshape = rawSpots['spot'][0].shape
        self.perFiberCoeffs = SpotCache(self._spotCacheDir(rawSpotPath))
        self.perFiberCoeffs.create(dict(wave=self.wave,
                                        fiber=self.fiber,
                                        spotScale=self.spotScale,
                                        spotsInfo=spotsInfo))

        splineType = spInterp.InterpolatedUnivariateSpline
        spotWaves = np.unique(self.wave)
        spots = rawSpots['spot']
        for fid in np.unique(self.fiber):
            self.logger.info('constructing splines for fiber %d', fid)
            fib_w = np.where(self.fiber == fid)[0]
            spotPixelCoeffs = np.zeros(imshape, dtype='O')
            fiberSpots = spots[fib_w]
            assert len(spotWaves) == len(fiberSpots)

            for ix in range(imshape[0]):
                for iy in range(imshape[1]):
                    spotPixelCoeffs[iy, ix] = self.buildSpline(splineType,
                                                               None, spotWaves, fiberSpots[:, iy, ix])
            xcCoeffs = self.buildSpline(splineType,
                                        None, spotWaves,
                                        rawSpots['spot_xc'][fib_w])
            ycCoeffs = self.buildSpline(splineType,
                                        None, spotWaves,
                                        rawSpots['spot_yc'][fib_w])

            self.perFiberCoeffs[fid] = SpotCoeffs(spotPixelCoeffs, xcCoeffs, ycCoeffs,
                                                  spots=fiberSpots)

    def loadFromSpots(self, spotType='jeg', spotIDs=None, spotArgs=None):
        """ Generate ourself from a semi-digested pfs_instdata spot file. 

        """

        self.logger.info("reading and interpolating %s PSF spots: %s..." % (spotType, spotIDs))

        if spotArgs is None:
            spotArgs = dict()

        if spotType == 'jeg':
            from . import jegSpots
            reload(jegSpots)

            rawSpotPath, _ = jegSpots.resolveSpotPathSpec(spotIDs)

            cached = self.primeSpotCache(spotIDs, rawSpotPath)
            if not cached:
                self.logger.warn('creating per-fiber spot cache for %s. This takes several minutes...' %
                                 rawSpotPath)
                self.createSpotCache(spotIDs, rawSpotPath, spotArgs)
        else:
            spotsFilepath = SplinedPsf.psfSpotsFile(self.detector.armName, spotType, spotIDs)
            sf = pyfits.open(spotsFilepath, mode='readonly')
            rawSpots = sf[1].data

            self.spotScale = sf[0].header['pixscale']
            self.wave = rawSpots['wavelength']
            self.fiber = rawSpots['fiberIdx']
            self.spots = rawSpots['spot'].astype('float32')

            self.buildAllSplines(perFiberCoeffs=False)

            # Keep the raw info around, at least until production.
            self.rawSpots = rawSpots

    def buildAllSplines(self, perFiberCoeffs):
        """ Build the spot-to-spot and spot center splines. """

        imshape = self.spots.shape[1:]

        # Make a spline for each pixel. The spots are centered on the output grid,
        # and we track the xc,yc offset separately.

        # XXX - Check that the input is properly sorted.
        xx = np.unique(self.fiber)
        yy = np.unique(self.wave)

        splineType = spInterp.RectBivariateSpline
        pixelCoeffs = np.zeros(imshape, dtype='O')
        for ix in range(imshape[0]):
            for iy in range(imshape[1]):
                pixelCoeffs[iy, ix] = self.buildSpline(splineType,
                                                       xx, yy,
                                                       self.spots[:, iy, ix].reshape(len(xx), len(yy)))
        xc = self.rawSpots['spot_xc']
        yc = self.rawSpots['spot_yc']

        # Add slit shift offsets.
        slitOffset = self.calcSlitOffset()
        xc += slitOffset[0]
        yc += slitOffset[1]

        xcCoeffs = self.buildSpline(splineType,
                                    xx, yy, xc.reshape(len(xx), len(yy)))
        ycCoeffs = self.buildSpline(splineType,
                                    xx, yy, yc.reshape(len(xx), len(yy)))
        self.allFiberCoeffs = SpotCoeffs(pixelCoeffs, xcCoeffs, ycCoeffs)

    def calcSlitOffset(self):
        """ Given .slitOffset, calculate spot offsets for .xc. """

        if self.slitOffset[1] != 0:
            raise RuntimeError("Cannot yet shift in wavelength!")

        if self.slitOffset[0] == 0 and self.slitOffset[1] == 0:
            return 0.0, 0.0

        offsetData = np.genfromtxt('slitmove.dat',dtype='f4', names=True)
        offsetWaves = np.unique(offsetData['lambda'])
        offsetShifts = np.unique(offsetData['yslt'])

        offsetDetX = offsetData['xdet']
        offsetDetY = offsetData['ydet']

        # Totally cheat for now -- CPL
        return (self.slitOffset[0] * offsetData['dyd_dys'].mean(),
                self.slitOffset[1] * offsetData['dxd_dxs'].mean())
        
        # x-y swap here: the physical detectors have the spectra going l-r.
        xInterp = spInterp.RectBivariateSpline(offsetShifts, offsetWaves,
                                               offsetDetY.reshape(len(offsetShifts),
                                                                  len(offsetWaves)))
        yInterp = spInterp.RectBivariateSpline(offsetShifts, offsetWaves,
                                               offsetDetX.reshape(len(offsetShifts),
                                                                  len(offsetWaves)))
        
    def buildSpline(self, splineType, x, y, z):
        if splineType is spInterp.RectBivariateSpline:
            return splineType(x, y, z)
        elif splineType is spInterp.InterpolatedUnivariateSpline:
            return splineType(y, z)
        elif splineType is spInterp.RegularGridInterpolator:
            return splineType((x, y), z, method='linear')
        elif splineType is float:
            return z[0]
        else:
            raise RuntimeError('buildSpline: unknown spline type: %s' % (splineType))

    def evalSpline(self, spline, x, y):
        splineType = spline.__class__
        if splineType is spInterp.RectBivariateSpline:
            return spline(x, y)
        elif splineType is spInterp.InterpolatedUnivariateSpline:
            return spline(y)
        elif splineType is spInterp.RegularGridInterpolator:
            x = np.asarray(x).flatten()
            y = np.asarray(y).flatten()
            xx,yy = np.meshgrid(x, y, indexing='ij')
    
            zz = np.array((xx,yy)).T.reshape(x.shape[0]*y.shape[0], 2)
    
            return np.atleast_2d(spline(zz))
        elif splineType is float:
            return np.asarray(spline)
        else:
            # print("warning: evalSpline: unknown spline type: %s'" % (spline))
            return spline(x, y)

    def spotGridImage(self):
        """ Return an array showing the rawspots in their given locations. """

        pass
        
    def setConstantSpot(self, spot=None, spotScale=None):
        """ Set ourself to use a single PSF over the whole focal plane. 

        Uses the existing .fiber and .wave maps. 
        """

        if spot is None:
            spotFiber = 0
            spotWave = sorted(self.wave)[len(self.wave)//2]
            spotId = np.where((self.fiber == spotFiber) & (self.wave == spotWave))[0]
            spot = self.spots[int(spotId)]

            self.logger.warn("set constant PSF to spot %s (fiber %d, wave %d)", spotId, spotFiber, spotWave)
        else:
            spotShape = self.spots[0].shape
            spotRad = spotShape[0]/2
            alpha = spot
            rad = spotgames.radToImageIndices(spotRad, sampling=0.1)
            spot = spotgames.gaussianFiber(rad, sigma=1.0, alpha=alpha)
            spot /= spot.sum()
            self.logger.warn("set constant PSF to given spot, with alpha=%s. shape=%s from %s (%s..%s)" %
                             (alpha,
                              spot.shape,
                              spotShape,
                              rad.min(), rad.max()))
            
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
                coeffs[iy, ix] = spotval
        self.coeffs = coeffs

    def setConstantX(self):
        """ Force spot x positions to be for the center spot on their trace. """

        xx = np.unique(self.fiber)
        yy = np.unique(self.wave)

        xc = self.xc.reshape(len(xx), len(yy))
        xc[:,:] = xc[:,len(yy)//2:len(yy)//2+1]
        self.xcCoeffs = self.buildSpline(self.xcCoeffs.__class__,
                                         xx, yy, xc)
        self.logger.warn("set constant X")

    def makeDetectorMap(self, fname, obsdate=None):
        """ Create a DetectorMap file for DRP. 

        Why here? We know both about the millimeters from the optical model and
        about the pixels of the detector. 

        """
        import lsst.afw.geom as afwGeom
        import lsst.daf.base as dafBase
        import pfs.drp.stella as drpStella

        pixelScale = self.detector.config['pixelScale']
        ccdSize = self.detector.config['ccdSize']

        # Per RHL, we want the _detector_ geometry here.
        bbox = afwGeom.BoxI(afwGeom.PointI(0,0), afwGeom.PointI(ccdSize[1]-1, ccdSize[0]-1))

        fiberIds = np.unique(self.fiber)
        fiber0_w = np.where(self.fiber == min(self.fiber))
        nKnots = len(fiber0_w[0])

        lam0 = self.spotsInfo['LAM0']
        dlam = self.spotsInfo['LAMINC']
        lams = np.linspace(lam0, lam0+(nKnots-1)*dlam, nKnots, dtype=np.float32)

        allYKnots = []
        allXCenters = []
        for holeId in fiberIds:
            coeffs = self.getCoeffs(holeId)

            yKnot = self._focalPlaneYToDetectorY(coeffs.ycCoeffs.get_coeffs(), correctToLL=True) / pixelScale
            xc = self._focalPlaneXToDetectorX(coeffs.xcCoeffs.get_coeffs(), correctToLL=True) / pixelScale
            allYKnots.append(yKnot.astype(np.float32))
            allXCenters.append(xc.astype(np.float32))

            assert len(lams) == len(yKnot)
            assert len(lams) == len(xc)

            midY = len(lams)//2
            self.logger.info("hole %d: xcKnot, xc, wl: %s %s %s" % (holeId,
                                                                    yKnot[midY], xc[midY], lams[midY]))

        detMap = drpStella.DetectorMap(bbox, fiberIds.astype('i4'), allYKnots, allXCenters,
                                       allYKnots, [lams]*len(fiberIds))

        if obsdate is None:
            now = time.gmtime()
        else:
            try:
                now = time.strptime(obsdate, '%Y-%m-%dT%H:%M:%S')
            except:
                now = time.strptime(obsdate, '%Y-%m-%d')

        obsdate = time.strftime('%Y-%m-%dT%H:%M:%S', now)
        calibDate = time.strftime('%Y-%m-%d', now)
        ccd = 3*(self.detector.spectrograph - 1) + dict(b=0, r=1, m=1, n=2)[self.detector.arm]
        metadata = detMap.getMetadata()
        metadata.addString('DATE-OBS', obsdate)
        metadata.addString('CALIB_ID', 'arm=%s spectrograph=%d ccd=%d filter=%s calibDate=%s visit0=0' %
                           (self.detector.arm, self.detector.spectrograph, ccd, self.detector.arm, calibDate))

        detMap.writeFits(fname)

        return detMap, fiberIds, allYKnots, allXCenters
