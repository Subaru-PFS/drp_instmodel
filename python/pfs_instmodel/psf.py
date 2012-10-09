import numpy
import scipy

class Psf(object):
    def __init__(self, band, detector):
        self.band = band
        self.detector = detector
    
    def psfAt(self, fiber, wave):
        """ Instantiate a single PSF at the given position.
        """
        
        ims, ctrs =  self.psfsAt([fiber], [wave])
        return ims[0]
    
    def fiberImages(self, fibers, spectra=None, outImg=None, waveRange=None, everyNthPsf=1):
        """ Return an image of the """

        if outImg == None:
            outImg = self.detector.simBias().image

        if spectra == None:
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

    
    
