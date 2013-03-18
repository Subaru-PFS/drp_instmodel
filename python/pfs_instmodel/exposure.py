import numpy
import pyfits

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
    
    def addFlux(self, im, ivar=None, mask=None, addPlane=None):
        """ Add some flux to ourselves.

        addPlane - optionally track the new flux as a separate image plane. This is
                   mostly for simulations, where we want to save the components of an 
                   image individually.
        """
        
        self.image += im
        if ivar:
            self.ivar += ivar
        if mask:
            self.mask += mask
        if addPlane:
            assert isinstance(addPlane, basestring), "addPlane argument must be a string"
            if hasattr(self.planes, addPlane):
                raise KeyError("plane named %s already exists" % (addPlane)

            # Go back and be fancy with properties -- CPL
            setattr(self, addPlane, im)
        
    def setImage(self, im, ivar):
        # assert im.shape == self._image.shape, "cannot overwrite with image of new shape"
        
        self._image = im
        self._ivar = ivar
