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
    def __init__(self, band, sky=None, psf=None, simID=None, addNoise=True, dtype='i4'):
        self.detector = pfsDet.Detector(band)
        self.sky = sky if sky else pfsSky.StaticSkyModel(band)
        self.psf = psf if psf else pfsPsf.SplinedPsf(self.detector, spotID=simID)
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

    def writeTo(self, outputFile, addNoise=True):
        import fitsio

        print("output to %s, addNoise=%s, dtype=self.exposure" % (outputFile, addNoise))

        self.exposure.writeto(outputFile, addNoise=addNoise)
        
        waveImage = self.waveImage()
        fitsio.write(outputFile, waveImage, extname='wavelengths', compress='RICE')
        

def fiberInfo(self):
    """ Return a single numpy array containing what we know about the fibers. """

    nFibers = len(self.fibers)
    s_i = sorted(self.fibers)

    dtype = numpy.dtype([('id','i2')])

