import numpy

class Exposure(object):
    def __init__(self, detector, doNew=True, dtype='u2'):
        """ """
        
        self.detector = detector

        if doNew:
            self._image = numpy.zeros(self.detector.config['ccdSize'], dtype=dtype)
            self._mask = numpy.zeros(self.detector.config['ccdSize'], dtype='u4')
            self._ivar = numpy.zeros(self.detector.config['ccdSize'], dtype='f4')
            self.planes = [('flux', 'image'),
                           ('ivar', 'ivar'),
                           ('mask', 'mask')]

    def __str__(self):
        return "Exposure(rows=%s, cols=%s, dtype=%s)" % (self._image.shape[0],
                                                         self._image.shape[1],
                                                         self._image.dtype)
    @property 
    def image(self):
        return self._image
    
    @property 
    def ivar(self):
        return self._ivar
    
    @property 
    def mask(self):
        return self._mask

    @property 
    def shape(self):
        return self._image.shape
    
    def addFlux(self, addIm, outSlice=None, addNoise=True, addPlane=None):
        """ Add some flux to ourselves.

        """
        
        if outSlice is not None:
            outIm = self._image[outSlice]
        else:
            outIm = self._image

        self._addFlux(outIm, addIm, addNoise=addNoise, addPlane=addPlane)

    def _addFlux(self, outIm, inIm, addNoise=True, addPlane=None):
        outIm += inIm

        if addNoise:
            noise = numpy.rint(numpy.sqrt(inIm) * numpy.random.normal(size=inIm.shape)).astype('i4')

            # Do we need a soft bias?
            w = numpy.where((outIm + noise) < 0)
            outIm += noise
            outIm[w] = 0


    def setImage(self, im, ivar):
        # assert im.shape == self._image.shape, "cannot overwrite with image of new shape"
        
        self._image = im
        self._ivar = ivar
