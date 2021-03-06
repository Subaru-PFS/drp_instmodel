import glob
import logging
import numpy
import os
import time
import warnings

import astropy.io.fits as pyfits
from .utils import geom
from .arm import Arm


class Exposure(object):
    def __init__(self, detector, 
                 dtype='u2'):
        """ """
        
        self.detector = detector
        self.dtype = dtype
        self.logger = logging.getLogger('exposure')
        self.logger.setLevel(logging.INFO)

        self._flux = numpy.zeros(self.detector.config['ccdSize'], dtype='f4')
        self._mask = numpy.zeros(self.detector.config['ccdSize'], dtype='u2')
        self._sky = numpy.zeros(self.detector.config['ccdSize'], dtype='f4')
        self.planes = dict()
        self.addPlane('flux', self._flux)
        self.addPlane('mask', self._mask)
        self.addPlane('sky', self._sky)
        self.pixelImage = None

    def clone(self):
        """Return a copy of ``self``

        This is probably most useful *before* calling the ``readout``,
        ``writeto``, ``loadBias`` or ``loadFlat`` methods, which add details
        specific to one particular instance to this object.
        """
        new = type(self)(self.detector, dtype=self.dtype)
        for attr in ("_flux", "_mask", "_sky"):
            getattr(new, attr)[:] = getattr(self, attr)
        for plane in self.planes:
            if plane not in new.planes:
                new.addPlane(self.planes[plane])
        return new

    def __str__(self):
        return "Exposure(rows=%s, cols=%s, dtype=%s, planes=%s)" % (self._flux.shape[0],
                                                                    self._flux.shape[1],
                                                                    self._flux.dtype,
                                                                    list(self.planes.keys()))
    @property 
    def mask(self):
        return self._mask

    @property 
    def shape(self):
        return self._flux.shape
    
    def addFlux(self, addIm, outSlice=None, addPlane=None):
        """ Add some flux to ourselves. 

        """
        
        if outSlice is not None:
            outIm = self._flux[outSlice]
        else:
            outIm = self._flux

        outIm += addIm

        if addPlane is not None:
            self.addPlane(addPlane, addIm)

    def addSky(self, skyImage, outSlice=None):
        outSky = self._sky[outSlice] if outSlice is not None else self._sky
        outSky += skyImage

    def writeto(self, outputFile, doCombine=True, doWriteAll=True,
                exptime=1.0, pfsDesignId=0x0, timestamp="2020-01-01T00:00:00.0",
                addNoise=True, compress='RICE',
                realBias=None, realFlat=None,
                addCards=(),
                imagetyp=None, allOutput=False):

        hdulist = pyfits.HDUList()

        self.readout(addNoise=addNoise,
                     realBias=realBias, realFlat=realFlat, exptime=exptime)
        if realBias is not None:
            outIm = self.biasExp.replaceActiveFlux(self.pixelImage, leadingRows=True)
            hdr = self.biasExp.header
        else:
            outIm = self.pixelImage
            hdr = pyfits.Header()

        print("writing to %s, addNoise=%s, imagetyp=%s, %s, dtype=%s" % (outputFile, addNoise,
                                                                         imagetyp, type(outIm), outIm.dtype))

        hdr.set('EXPTIME', float(exptime))
        hdr.set('DATE-OBS', timestamp)
        hdr.set('W_SIMBIA', realBias)
        hdr.set('W_SIMFLA', realFlat)
        hdr.set('W_PFDSGN', pfsDesignId)

        for c in addCards:
            hdr.set(*c)

        hdu0 = pyfits.CompImageHDU(outIm, name='image')
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", pyfits.verify.VerifyWarning)  # creating HIERARCH keys
            hdu0.header.extend(hdr)

        # hdu0.data = outIm
        hdulist.append(hdu0)
            
        if allOutput:
            hdulist.append(pyfits.CompImageHDU(self.planes['mask'], name='mask'))
            hdulist.append(pyfits.CompImageHDU(self.planes['bias'], name='bias'))
            if "shotnoise" in self.planes:
                hdulist.append(pyfits.CompImageHDU(self.planes['shotnoise'], name='shotnoise'))
            hdulist.append(pyfits.CompImageHDU(self.planes['sky'], name='sky'))
            if realBias:
                hdulist.append(pyfits.CompImageHDU(self.pixelImage, name='active'))

        hdulist.update_extend()
        hdulist.writeto(outputFile, checksum=True, overwrite=True)

        return hdulist
    
    def loadBias(self, biasID):
        """ Load a real detector bias, return its active image. """

        if biasID is None:
            return None

        if isinstance(biasID, str):
            filepath = biasID
        else:
            dataRoot = os.environ.get('DRP_INSTDATA_DIR', '.')
            fileglob = os.path.join(dataRoot, 'data', 'pfs', 'biases', 'PF?A0*%d%d%d.fits' %
                                    (biasID, 1, 2 if self.detector.arm == Arm.RED else 1))

            print("looking for biases %s" % (fileglob))
            filepaths = glob.glob(fileglob)
            filepath = filepaths[0]

        print("loading bias %s" % (filepath))
        self.biasExp = geom.Exposure(obj=filepath)
        print("  bias geom: %s" % (self.biasExp))

        # Paper over a temporary geometry botch in the DA
        if self._flux.shape[0] == 4174:
            self.biasExp.leadinRows = 50

        return self.biasExp.finalImage(leadingRows=False)

    def loadFlat(self):
        """ Load a real detector flat, return its active image. """

        dataRoot = os.environ.get('DRP_INSTDATA_DIR', '.')
        fileglob = os.path.join(dataRoot, 'data', 'pfs', 'flats', 'pfsFlat-*-%d%s.fits' %
                                (1, self.detector.arm.value))

        print("looking for flats %s" % (fileglob))
        filepaths = glob.glob(fileglob)
        filepath = filepaths[0]
        print("fetching flat %s" % (filepath))
        flat = pyfits.getdata(filepath)

        return flat
    
    def readout(self, exptime=1.0, addNoise=True,
                realBias=None, realFlat=None):
        
        if self.pixelImage is not None:
            return

        self._flux *= exptime
        self._sky *= exptime
        if addNoise:
            print("adding noise=%s" % (addNoise))
            lo_w = numpy.where(self._flux < 0)
            if len(lo_w[0]) > 0:
                print("%d low image pixels, min=%f" % (len(lo_w[0]), 
                                                 self._flux.min()))
                self._flux += -self._flux.min()
            numLowSky = (self._sky < 0).sum()
            if numLowSky > 0:
                print("%d low sky pixels, min=%f" % (numLowSky, self._sky.min()))
                self._sky += -self._sky.min()

            noisyFlux = numpy.random.poisson(self._flux + self._sky)
            noise = noisyFlux - self._flux
            self.addPlane('shotnoise', noise)
        else:
            noisyFlux = self._flux

        self.detector.readout(self, noisyFlux, exptime=exptime,
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

