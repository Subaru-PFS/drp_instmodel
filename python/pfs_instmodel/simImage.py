import logging
import numpy

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
    def __init__(self, band, sky=None, psf=None, simID=None,
                 addNoise=True, dtype='i4',
                 constantPsf=False, constantX=False,
                 logger=None):

        if logger is None:
            logger = logging.getLogger()
        self.logger = logger

        self.detector = pfsDet.Detector(band)
        self.sky = sky if sky else pfsSky.StaticSkyModel(band)
        self.psf = psf if psf else pfsPsf.SplinedPsf(self.detector, spotID=simID,
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

    def addFibers(self, fibers, spectra, waveRange=None,
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
        for i, fiber in enumerate(fibers):
            parts = self.psf.fiberImage(fiber, spectra[i], outExp=self.exposure,
                                        waveRange=waveRange,
                                        shiftPsfs=shiftPsfs)
            self.fibers[fiber] = dict(spectrum=spectra[i],
                                      geometry=parts[3])

        return self.exposure

    def waveImage(self):
        """ """

        nFibers = len(self.fibers)
        f_i = sorted(self.fibers)
        
        waveArr = numpy.zeros((nFibers, self.detector.config['ccdSize'][1]), dtype='f4')
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
        
    def writeTo(self, outputFile, addNoise=True,
                compress='RICE', allOutput=False,
                imagetyp=None, realBias=None):
        import fitsio

        print("output to %s, addNoise=%s, realBias=%s, dtype=self.exposure" %
              (outputFile, realBias, addNoise))
        
        self.exposure.writeto(outputFile, addNoise=addNoise,
                              realBias=realBias, imagetyp=imagetyp,
                              compress=compress, allOutput=allOutput)
        
        waveImage = self.waveImage()
        fitsio.write(outputFile, waveImage, extname='wavelengths', compress='RICE')

        lineGeometry = self.lineGeometry()
        fitsio.write(outputFile, lineGeometry, extname='lines')
        
def fiberInfo(self):
    """ Return a single numpy array containing what we know about the fibers. """

    nFibers = len(self.fibers)
    s_i = sorted(self.fibers)

    dtype = numpy.dtype([('id','i2')])

