from importlib import reload

import logging
import numpy

from concurrent.futures import ProcessPoolExecutor
from concurrent.futures import as_completed

import pfs.instmodel.detector as pfsDet
import pfs.instmodel.sky as pfsSky
import pfs.instmodel.splinedPsf as pfsPsf
from pfs.instmodel.spectrum import NullSpectrum

"""
Example
_______
>>> simg = SimImage('IR')
>>> fibers = numpy.concatenate([numpy.arange(5),
                                numpy.arange(5) + 100,
                                numpy.arange(5) + 290])
>>> irImage = irPsf.fiberImages(fibers,
                                spectra=[irSky]*len(fibers))
"""
class SimImage(object):
    def __init__(self, detector, sky, psf=None, simID=None,
                 dtype='i4',
                 everyNth=20,
                 constantPsf=False, constantX=False,
                 slitOffset=(0.0, 0.0),
                 logger=None):

        if logger is None:
            logger = logging.getLogger()
        self.logger = logger

        self.detector = pfsDet.Detector(detector)
        self.sky = sky
        self.psf = psf if psf else pfsPsf.SplinedPsf(self.detector, spotID=simID,
                                                     everyNth=everyNth,
                                                     slitOffset=slitOffset,
                                                     logger=logger)
        if constantPsf:
            self.psf.setConstantSpot(constantPsf if constantPsf else None)
        if constantX:
            self.psf.setConstantX()

        self.exposure = self.detector.makeExposure(dtype=dtype)
        self.fibers = {}

    def clone(self):
        """Return a copy of ``self``"""
        new = type(self)(self.detector.detectorName, self.sky, self.psf, self.exposure.dtype)
        new.exposure = self.exposure.clone()
        new.fibers = self.fibers.copy()
        return new

    @property
    def image(self):
        return self.exposure.image

    def addFibers(self, fibers, spectra, doSkyForFiber, waveRange=None,
                  shiftPsfs=True):
        """ Add images of the given fibers. 

        Parameters
        ----------
        fibers : array_like
            the fiber IDs to add images of
        spectra : array_like
            the spectra for the given fibers.
        waveRange : (minWave, maxWave), optional
            limit the spectra to the given inclusive wavelength range.

        Returns
        ------
        image - a full image for the given detector.

        Notes
        -----
        The fiber images are added to our internal image, so multiple calls should just add flux.

        The spectra arg is currently something which has a Spectrum signature (flux =  __call__(wave)),
        but should hoisted up to have the full probe schema.
        
        I believe that sky spectra could be added entirely differently from object spectra. So for
        the following 
          * 
        """

        for i, (fiber, doSky) in enumerate(zip(fibers, doSkyForFiber)):
            self.addSingleFiber(fiber, spectra[i], doSky, waveRange, shiftPsfs)
        return self.exposure

    def addSingleFiber(self, fiber, spectrum, doSky, waveRange=None, shiftPsfs=True):
        """Add image of the given fiber"""
        parts = self.psf.fiberImage(fiber, spectrum*self.sky.getExtinction(),
                                    waveRange=waveRange, shiftPsfs=shiftPsfs)
        fiberImage, outImgOffset, psfToSpotPixRatio, geometry = parts
        if fiberImage is not None:
            # transfer flux from oversampled fiber image to final resolution output image
            self.psf.addOversampledImage(fiberImage, self.exposure, outImgOffset, psfToSpotPixRatio)

        self.fibers[fiber] = dict(spectrum=spectrum, geometry=geometry)

        if doSky:
            skyParts = self.psf.fiberImage(fiber, self.sky.getSky(), waveRange=waveRange, shiftPsfs=shiftPsfs)
            skyImage = skyParts[0]
            skyOffset = skyParts[1]
            if skyImage is not None:
                self.psf.addOversampledImage(skyImage, self.exposure, skyOffset, psfToSpotPixRatio)

    def waveImage(self):
        """ """

        nFibers = len(self.fibers)
        f_i = sorted(self.fibers)
        
        waveArr = numpy.zeros((nFibers, self.detector.config['ccdSize'][0]), dtype='f4')
        waveArr[:] = numpy.nan
        for i in range(nFibers):
            rows, waves = self.psf.wavesForRow([f_i[i]])
            select = (rows >= 0) & (rows < self.detector.config['ccdSize'][0])
            waveArr[i] = waves[select]

        return waveArr
        
    def writeTo(self, outputFile=None, addNoise=True,
                exptime=1.0, pfsDesignId=0x0, timestamp="2020-01-01T00:00:00.0",
                compress='RICE', allOutput=False,
                imagetyp=None, realBias=None, realFlat=None,
                addCards=None):
        if addCards is None:
            addCards = []
        if outputFile is None:
            from .utils import SeqPath

            baseTemplate = '%(filePrefix)s%(seqno)06d'
            self.fileMgr = SeqPath.NightFilenameGen('/data/pfsSim',
                                                    filePrefix='PFFA',
                                                    filePattern="%s%s.fits" % (baseTemplate,
                                                                               self.detector.detectorName))
            outputFile = self.fileMgr.getNextFileset()[0]
        if realBias is True:
            realBias = int(outputFile[-8])

        print("output to %s, addNoise=%s, realBias=%s" %
              (outputFile, addNoise, realBias))

        addCards += self.psf.getCards()
        
        self.exposure.writeto(outputFile, addNoise=addNoise,
                              exptime=exptime, pfsDesignId=pfsDesignId, timestamp=timestamp,
                              realBias=realBias, realFlat=realFlat,
                              imagetyp=imagetyp,
                              addCards=addCards,
                              compress=compress, allOutput=allOutput)
