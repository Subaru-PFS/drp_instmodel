import numpy

class Exposure(object):
    def __init__(self, detector, 
                 doNew=True, 
                 addBias=False,
                 addNoise=False,
                 dtype='u2'):
        """ """
        
        self.detector = detector
        self.addNoise = addNoise
        self.dtype = dtype

        if doNew:
            self._flux = numpy.zeros(self.detector.config['ccdSize'], dtype='f4')
            self._mask = numpy.zeros(self.detector.config['ccdSize'], dtype='u2')
            self.planes = dict()
            self.addPlane('flux', self._flux)
            self.addPlane('mask', self._mask)
            self.pixelImage = None

        if addBias:
            bias = self.detector.addBias(self)
            self.addPlane('bias', bias)

    def __str__(self):
        return "Exposure(rows=%s, cols=%s, dtype=%s, planes=%s)" % (self._flux.shape[0],
                                                                    self._flux.shape[1],
                                                                    self._flux.dtype,
                                                                    self.planes.keys())
    @property 
    def mask(self):
        return self._mask

    @property 
    def shape(self):
        return self._flux.shape
    
    def addFlux(self, addIm, outSlice=None, addNoise=True, addPlane=None):
        """ Add some flux to ourselves. Optionally 

        """
        
        if outSlice is not None:
            outIm = self._flux[outSlice]
        else:
            outIm = self._flux

        outIm += addIm

        if addPlane is not None:
            self.addPlane(addPlane, addIm)

        if addNoise and self.addNoise:
            noise = numpy.rint(numpy.sqrt(addIm.clip(0,addIm.max())) * numpy.random.normal(size=addIm.shape)).astype('i4')

            # Do we need a soft bias?
            w = numpy.where((outIm + noise) < 0)
            outIm += noise

            if len(w[0] > 0):
                minLow = min((outIm + noise)[w])
                print("%d low pixels (min=%0.4f,%0.4f,%0.4f) at %s" % (len(w[0]),
                                                                       minLow,
                                                                       min(addIm[w]),
                                                                       min(noise[w]),
                                                                       outSlice))

    def writeto(self, outputFile, doCombine=True, doWriteAll=True, addNoise=True, compress='RICE'):
        import fitsio

        print("writing to %s, addNoise=%s" % (outputFile, addNoise))

        self.readout(addNoise=addNoise)
        fitsio.write(outputFile, self.pixelImage, clobber=True, compress=compress)
        fitsio.write(outputFile, self.planes['mask'], extname='mask', compress=compress)
        fitsio.write(outputFile, self.planes['bias'], extname='bias', compress=compress)

    def readout(self, addNoise=True):
        if self.pixelImage is not None:
            return

        if addNoise:
            print("adding noise=%s" % (addNoise))
            lo_w = numpy.where(self._flux < 0)
            if len(lo_w[0]) > 0:
                print("%d low pixels, min=%f" % (len(lo_w[0]), 
                                                 self._flux.min()))
                self._flux[lo_w] = 0

            noisyFlux = numpy.random.poisson(self._flux)
            noise = noisyFlux - self._flux
            self.addPlane('shotnoise', noise)
        else:
            noisyFlux = self._flux
            
        self.detector.readout(self, noisyFlux)

    def XXXXaddNoise(self):
        noise = numpy.rint(numpy.sqrt(self.flux * numpy.random.normal(size=self.shape))).astype('i2')

        self.addPlane('shotnoise', noise)

    def addPlane(self, name, im):
        if name in self.planes:
            raise KeyError("exposure image plane %s already exists!" % (name))

        self.planes[name] = im

    def setImage(self, im, ivar):
        # assert im.shape == self._image.shape, "cannot overwrite with image of new shape"
        
        assert self._image.dtype == im.dtype

        self._image[:] = im
        self._ivar[:] = ivar

class SimExposure(Exposure):
    def __init__(self, detector, 
                 doNew=True, 
                 addBias=False,
                 addNoise=False,
                 dtype='u2'):
        """ """

        super(Exposure, self).__init__(doNew=True,
                                       addBias=False,
                                       addNoise=addNoise,
                                       dtype='f4')

