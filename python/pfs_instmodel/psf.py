import numpy
import scipy

class Psf(object):
    def __init__(self, detector):
        self.detector = detector
    
    def psfAt(self, fiber, wave):
        """ Instantiate a single PSF at the given position.
        """
        
        ims, ids, ctrs =  self.psfsAt([fiber], [wave])
        return ims[0]
    
    def fiberImages(self, fibers, spectra=None, outImg=None, waveRange=None, everyNthPsf=1):
        """ Return and/or place an image of the given spectra through the given fibers. 

        Parameters
        ----------
        fibers : int[NFIBERS] XXX - should be some useful fiber objects.
           The IDs of the fibers we want images through. 
        spectra : Spectrum[NFIBERS]
           
        """

        if outImg is None:
            outImg = self.detector.simBias().image

        if spectra is None:
            spectra = [None] * len(fibers)
        for i, fiber in enumerate(fibers):
            self.fiberImage(fiber, spectra[i], outImg=outImg, 
                            waveRange=waveRange, everyNthPsf=everyNthPsf)

        return outImg

    def scalePsf(self, rawPsf, xc, yc, doRescaleFlux=False, doDetails=False):
        """ Given an oversampled image scale it to the detector grid after shifting it. """

        psfPixelScale = self.detector.config['pixelScale'] / self.spotScale
            
        psf = scipy.ndimage.affine_transform(rawPsf, 
                                             [psfPixelScale, psfPixelScale],
                                             (yc, xc),
                                             output_shape=numpy.array(rawPsf.shape)/psfPixelScale,
                                             order=1)

        if doDetails:
            print ("           shift=%d/%0.3f, %d/%0.3f flux=%0.2f,%0.2f" % 
                   (round(xc*psfPixelScale), xc, 
                    round(yc*psfPixelScale), yc,
                    numpy.sum(rawPsf), numpy.sum(psf) * (psfPixelScale * psfPixelScale)))

        if doRescaleFlux:
            psf *= (psfPixelScale*psfPixelScale)

        return psf

    
    
