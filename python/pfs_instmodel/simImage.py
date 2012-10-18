import numpy

import pfs_instmodel.detector as pfsDet
import pfs_instmodel.sky as pfsSky
import pfs_instmodel.splinedPsf as pfsPsf

"""
Example
_______
>>> simg = SimImage('IR')
>>> fibers = numpy.concatenate([numpy.arange(5),
                                numpy.arange(5) + 100,
                                numpy.arange(5) + 290])
>>> irImage = irPsf.fiberImages(fibers,
                                spectra=[irSky]*len(fibers),
                                everyNthPsf=50)
"""

class SimImage(object):
    def __init__(self, band, sky=None, psf=None):
        self.detector = pfsDet.Detector(band)
        self.sky = sky if sky else pfsSky.StaticSkyModel(band)
        self.psf = psf if psf else pfsPsf.SplinedPsf(self.detector)
        self.image = None
        
    def addFibers(self, fibers, spectra, waveRange=None, everyNthPsf=1):
        """ Add images of the given fibers. 

        Parameters
        ----------
        fibers : array_like
            the fiber IDs to add images of
        spectra : array_like
            the spectra for the given fibers.
        waveRange : (minWave, maxWave), optional
            limit the spectra to the given inclusive wavelength range.
        everyNthPsf : int, optional
            only require the PSFs to vary on every Nth pixel. default=1

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
        if self.image == None:
            self.image = self.detector.simBias().image

        for i, fiber in enumerate(fibers):
            self.psf.fiberImage(fiber, spectra[i], outImg=self.image,
                                waveRange=waveRange, everyNthPsf=everyNthPsf)

        return self.image

