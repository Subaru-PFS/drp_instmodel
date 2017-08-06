import glob
import numpy
import os
import time

import astropy.io.fits as pyfits
from fpga import geom

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
    
    def addFlux(self, addIm, outSlice=None, addNoise=False, addPlane=None):
        """ Add some flux to ourselves. 

        """
        
        if outSlice is not None:
            outIm = self._flux[outSlice]
        else:
            outIm = self._flux

        outIm += addIm

        if addPlane is not None:
            self.addPlane(addPlane, addIm)

    def ts(self, t=None):
        if t is None:
            t = time.time()
            return time.strftime("%Y-%m-%dT%H:%M:%S", time.gmtime(t))
                                    
    def writeto(self, outputFile, doCombine=True, doWriteAll=True,
                addNoise=True, compress='RICE',
                realBias=None, realFlat=None,
                addCards=(),
                imagetyp=None, allOutput=False):

        self.readout(addNoise=addNoise, realBias=realBias, realFlat=realFlat)
        if realBias is not None:
            outIm = self.biasExp.replaceActiveFlux(self.pixelImage, leadingRows=False)
            hdr = self.biasExp.header
        else:
            outIm = self.pixelImage
            hdr = pyfits.Header()

        print("writing to %s, addNoise=%s, imagetyp=%s, %s, dtype=%s" % (outputFile, addNoise,
                                                                         imagetyp, type(outIm), outIm.dtype))
            
        hdulist = pyfits.HDUList()
        if imagetyp is not None:
            imageTyp, expTime = imagetyp.split(',')
            
            hdr.set('IMAGETYP', imageTyp)
            hdr.set('EXPTIME', float(expTime))
            hdr.set('DATE-OBS', self.ts(), 'Crude time')

        for c in addCards:
            hdr.set(*c)
            
        hdu0 = pyfits.CompImageHDU(outIm, header=hdr, name='image')
        # hdu0.data = outIm
        hdulist.append(hdu0)
            
        if allOutput:
            hdulist.append(pyfits.CompImageHDU(self.planes['mask'], name='mask'))
            hdulist.append(pyfits.CompImageHDU(self.planes['bias'], name='bias'))
            hdulist.append(pyfits.CompImageHDU(self.planes['shotnoise'], name='shotnoise'))
            if realBias:
                hdulist.append(pyfits.CompImageHDU(self.pixelImage, name='active'))
            
        hdulist.writeto(outputFile, checksum=True, clobber=True)

    def loadBias(self, biasID):
        """ Load a real detector bias, return its active image. """

        dataRoot = os.environ.get('DRP_INSTDATA_DIR', '.')
        filepath = os.path.join(dataRoot, 'data', 'pfs', 'PFSA00715%d%d%d.fits' %
                                (biasID, 9, 2 if self.detector.band == 'Red' else 1))

        print("loading bias %s" % (filepath))
        self.biasExp = geom.Exposure(obj=filepath)
        print("  bias geom: %s" % (self.biasExp))

        if self._flux.shape[0] == 4174:
            self.biasExp.leadinRows = 50
        
        return self.biasExp.finalImage(leadingRows=False)
        
    def loadFlat(self):
        """ Load a real detector flat, return its active image. """

        dataRoot = os.environ.get('DRP_INSTDATA_DIR', '.')
        fileglob = os.path.join(dataRoot, 'data', 'pfs', 'flats', 'pfsFlat-*-%d%s.fits' %
                                (1, self.detector.band[0].lower()))

        print("looking for flats %s" % (fileglob))
        filepaths = glob.glob(fileglob)
        filepath = filepaths[0]
        print("fetching flat %s" % (filepath))
        flat = pyfits.getdata(filepath)

        return flat
    
    def readout(self, addNoise=True,
                realBias=None, realFlat=None):
        
        if self.pixelImage is not None:
            return

        if addNoise:
            print("adding noise=%s" % (addNoise))
            lo_w = numpy.where(self._flux < 0)
            if len(lo_w[0]) > 0:
                print("%d low pixels, min=%f" % (len(lo_w[0]), 
                                                 self._flux.min()))

            noisyFlux = numpy.random.poisson(self._flux)
            noise = noisyFlux - self._flux
            self.addPlane('shotnoise', noise)
        else:
            noisyFlux = self._flux
            
        self.detector.readout(self, noisyFlux,
                              ontoBias=realBias, applyFlat=realFlat)

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

