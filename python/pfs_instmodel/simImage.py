from importlib import reload

import logging
import numpy

from concurrent.futures import ProcessPoolExecutor
from concurrent.futures import as_completed

import pfs_instmodel.detector as pfsDet
import pfs_instmodel.sky as pfsSky
import pfs_instmodel.splinedPsf as pfsPsf
reload(pfsPsf)

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
    def __init__(self, detector, sky=None, psf=None, simID=None,
                 addNoise=True, dtype='i4',
                 everyNth=20,
                 constantPsf=False, constantX=False,
                 slitOffset=(0.0, 0.0),
                 logger=None):

        if logger is None:
            logger = logging.getLogger()
        self.logger = logger

        self.detector = pfsDet.Detector(detector)
        self.sky = sky if sky else pfsSky.StaticSkyModel(self.detector.armName)
        self.psf = psf if psf else pfsPsf.SplinedPsf(self.detector, spotID=simID,
                                                     everyNth=everyNth,
                                                     slitOffset=slitOffset,
                                                     logger=logger)
        if constantPsf:
            self.psf.setConstantSpot(constantPsf if constantPsf else None)
        if constantX:
            self.psf.setConstantX()

        self.exposure = self.detector.makeExposure(dtype=dtype, addNoise=addNoise)
        self.fibers = {}

    @property
    def image(self):
        return self.exposure.image

    def addFibers(self, fibers, spectra, doSkyForFiber, skySpectrum, waveRange=None,
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
            parts = self.psf.fiberImage(fiber, spectra[i],
                                        waveRange=waveRange,
                                        shiftPsfs=shiftPsfs)

            fiberImage, outImgOffset, psfToSpotPixRatio, geometry = parts

            skyImage = None
            skyOffset = None
            if doSky:
                skyParts = self.psf.fiberImage(fiber, skySpectrum, waveRange=waveRange, shiftPsfs=shiftPsfs)
                skyImage = skyParts[0]
                skyOffset = skyParts[1]

            # transfer flux from oversampled fiber image to final resolution output image
            self.psf.addOversampledImage(fiberImage, self.exposure,
                                         outImgOffset, psfToSpotPixRatio,
                                         skyImage, skyOffset)

            self.fibers[fiber] = dict(spectrum=spectra[i],
                                      geometry=geometry)

        return self.exposure

    def waveImage(self):
        """ """

        nFibers = len(self.fibers)
        f_i = sorted(self.fibers)
        
        waveArr = numpy.zeros((nFibers, self.detector.config['ccdSize'][0]), dtype='f4')
        waveArr[:] = numpy.nan
        for i in range(nFibers):
            rows, waves = self.psf.wavesForRows([f_i[i]])
            waveArr[i][rows] = waves[0]

        return waveArr

    def lineGeometry(self):
        """ Return a recarray describing all the lines in the spectra with line lists. """
        
        combFibers = set()
        for f in self.fibers:
            if self.fibers[f]['spectrum'].__class__.__name__ == 'CombSpectrum':
                combFibers.add(f)

        nFibers = len(combFibers)
        fiberIds = sorted(combFibers)

        geomType = numpy.dtype([('fiberId','i2'),
                                ('xc','f8'), ('yc','f8'),
                                ('dxc','f4'), ('dyc','f4'),
                                ('wavelength', 'f8'),
                                ('flux', 'f4')])
        npoints = 0
        for f in fiberIds:
            g = self.fibers[f]['geometry']
            f_npoints = numpy.sum(g['flux'] > 0)
            npoints += f_npoints

        geomArr = numpy.zeros(npoints, dtype=geomType)
        geomDone = 0
        for i, f_i in enumerate(fiberIds):
            fiberGeom = self.fibers[f_i]['geometry']
            hasFlux = fiberGeom['flux'] > 0
            geomLen = numpy.sum(hasFlux)

            g_i = slice(geomDone, geomDone+geomLen)
            geomArr[g_i]['fiberId'] = f_i
            geomArr[g_i]['xc'] = fiberGeom['xc'][hasFlux] / self.detector.config['pixelScale']
            geomArr[g_i]['yc'] = fiberGeom['yc'][hasFlux] / self.detector.config['pixelScale']
            geomArr[g_i]['dxc'] = fiberGeom['dxc'][hasFlux]
            geomArr[g_i]['dyc'] = fiberGeom['dyc'][hasFlux]
            geomArr[g_i]['wavelength'] = fiberGeom['wavelength'][hasFlux]
            geomArr[g_i]['flux'] = fiberGeom['flux'][hasFlux]

            geomDone += geomLen
                
        return geomArr
        
    def writeTo(self, outputFile=None, addNoise=True, skySwindle=True,
                exptime=1.0, pfiDesignId=0x0,
                compress='RICE', allOutput=False,
                imagetyp=None, realBias=None, realFlat=None):
        import fitsio

        # import pyfits
        
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

        addCards = self.psf.getCards()
        
        self.exposure.writeto(outputFile, addNoise=addNoise, skySwindle=skySwindle,
                              exptime=exptime, pfiDesignId=pfiDesignId,
                              realBias=realBias, realFlat=realFlat,
                              imagetyp=imagetyp,
                              addCards=addCards,
                              compress=compress, allOutput=allOutput)

        if allOutput:
            waveImage = self.waveImage()
            fitsio.write(outputFile, waveImage, extname='wavelengths', compress='RICE')

            lineGeometry = self.lineGeometry()
            fitsio.write(outputFile, lineGeometry, extname='lines')
