from __future__ import print_function
from __future__ import division
from builtins import range
from past.builtins import basestring
from builtins import object
import itertools
import logging
import numpy as np

import astropy.io.fits as pyfits

class Exposure(object):
    def __init__(self, obj=None, dtype=None, nccds=2, copyExposure=False,
                 simpleGeometry=False, logLevel=logging.WARN):
        self.logger = logging.getLogger('geom.Exposure')
        self.logger.setLevel(logLevel)
        
        self.nccds = nccds
        self._setDefaultGeometry()

        if obj is None:
            self.image = None
            self.header = dict()
        elif isinstance(obj, Exposure):
            self.image = obj.image
            self.header = obj.header
            self.deduceGeometry()
            
            if copyExposure:
                self.image = self.image.copy()
                self.header = self.header.copy()
        elif isinstance(obj, basestring):
            ffile = pyfits.open(obj)
            self.image = ffile[0].data
            self.header = ffile[0].header
            self.deduceGeometry(simpleGeometry)
            self.image = self.fixEdgeColsBug(self.image)

        elif isinstance(obj, np.ndarray):
            self.image = obj.copy() if copyExposure else obj
            self.header = dict()
            self.deduceGeometry()
        else:
            raise RuntimeError("do not know how to construct from a %s" % (type(obj)))

        if dtype is not None:
            self.image = self.image.astype(dtype, copy=False)

    def __str__(self):
        return "Exposure(shape=%s, dtype=%s, rows=(%d,%d,%d) cols=(%d,%d,%d)*%d)" % (self.image.shape, self.image.dtype,
                                                                                     self.leadinRows, self.activeRows,
                                                                                     self.overRows,
                                                                                     self.leadinCols, self.activeCols,
                                                                                     self.overCols,
                                                                                     self.namps)
    
    def fixEdgeColsBug(self, image):
        """ Fix extra 0th pixel in early raw images.
        """
        
        try:
            vers = self.header['versions.FPGA']
            if vers >= 0xa071:
                return image
        except:
            pass
        
        try:
            flag = self.header.get('geom.edgesOK')
        except:
            return image

        if not flag:
            self.logger.info('fixing corner pixel')
            if False:
                fixedImage1 = np.ndarray(shape=image.shape, dtype=image.dtype)
                fixedImage1[:,:-1] = image[:,1:]
                fixedImage1[:,-1] = image[:,0]

                fixedImage = np.ndarray(shape=image.shape, dtype=image.dtype)
                fixedImage[:-1,:] = fixedImage1[1:,:]
                fixedImage[-1,:] = fixedImage1[0,:]
            else:
                fixedImage = np.ndarray(shape=image.size, dtype=image.dtype)
                fixedImage[:-1] = image.flat[1:]
                fixedImage[-1] = fixedImage[-2]
                fixedImage.shape = image.shape
                
            return fixedImage
        else:
            return image

    def parseHeader(self):
        try:
            self.namps = self.header['geom.namps']
            self.nccds = self.namps//2
            
            self.leadinRows = self.header['geom.rows.leadin']
            self.overRows = self.header['geom.rows.overscan']
            self.leadinCols = self.header['geom.cols.leadin'] 
            self.overCols = self.header['geom.cols.overscan']
            self.readDirection = self.header['geom.readDirection']

            self.ccdRows = self.image.shape[0] - self.overRows
            self.ampCols = self.image.shape[1]//self.namps - self.overCols

            return True
        except Exception as e:
            self.logger.warn("FAILED to parse geometry cards: %s", e)
            self.logger.warn(("header: %s", self.header))
            return False
        
    def deduceGeometry(self, simple=False):
        """ Use .image to generate geometry for an unbinned full-frame. 

        Args
        ----
        simple : bool
           If we do not have geometry in the header, generate one from the image shape,
           assuming that it has .namps
        """

        hdrOk = self.parseHeader()
        if not hdrOk:
            self.header['geom.edgesOK'] = True
            self._setDefaultGeometry()
            if simple:
                self.overCols = self.overRows = 0
                self.leadinCols = self.leadinRows = 0
                self.readDirection = 0
                self.ccdRows = self.image.shape[0]
                self.ampCols = self.image.shape[1]//self.namps
        
        imh,imw = self.image.shape
        if (self.ampCols + self.overCols)*self.namps != imw:
            raise RuntimeError("Strange geometry: %d amps * (%d cols + %d overscan cols) != image width %d)" %
                               (self.namps, self.ampCols, self.overCols, imw))
                               
    def _setDefaultGeometry(self):
        self.ampCols = 520
        self.ccdRows = 4224
        self.overCols = 32
        self.overRows = 76
        self.leadinCols = 8
        self.leadinRows = 48
        
        self.readDirection = 0b10101010
        
        self.namps = 4 * self.nccds

    def forceLeadinRows(self, newLeadinRows):
        self.leadinRows = newLeadinRows
        if 'geom.rows.leadin' in self.header:
            self.header['geom.rows.leadin'] = newLeadinRows


    @property
    def expType(self):
        return self.header.get('IMAGETYP', 'unknown').strip()
    
    @property
    def expTime(self):
        return self.header.get('EXPTIME', -1.0)

    @property
    def activeRows(self):
        return self.ccdRows - self.leadinRows
        
    @property
    def activeCols(self):
        return self.ampCols - self.leadinCols
        
    @property
    def nrows(self):
        return self.ccdRows + self.overRows

    @property
    def ncols(self):
        return self.ampCols + self.overCols

    def ampExtents(self, ampId, leadingCols=False, leadingRows=False, overscan=False):
        """ Return the y,x slices covering the given amp. 
        
        Args
        ----
        ampId : int
           The index of the amp. 0..self.namps
        leadingCols : bool
           Whether to include the leadin columns.
        leadingRows : bool
           Whether to include the leadin rows
        overscan : bool
           Whether to include the overscan pixels.

        """
        x0 = ampId*self.ncols + self.leadinCols*(not leadingCols)
        x1 = (ampId + 1)*self.ncols
        if not overscan:
            x1 -= self.overCols

        if not self.readDirection or (self.readDirection & (1 << ampId)) == 0:
            xr = slice(x0, x1)
        else:
            xr = slice(x1-1, x0-1, -1)

        yr = slice(self.leadinRows*(not leadingRows), self.ccdRows + overscan*self.overRows)

        return yr, xr

    def finalAmpExtents(self, ampId, leadingRows=True):
        """ Return the y,x slices for the  
        
        Args
        ----
        ampId : int
           The index of the amp. 0..self.namps
        leadingRows : bool
           Whether to include the leadin rows.
        """
        
        x0 = ampId*self.activeCols
        x1 = x0 + self.activeCols

        xr = slice(x0, x1)

        y0 = (not leadingRows)*self.leadinRows
        y1 = y0 + self.activeRows + leadingRows*self.leadinRows
        yr = slice(y0, y1)
        
        return yr, xr

    def ampImage(self, ampId, im=None, leadingCols=False, leadingRows=False):
        """ Return the image for a single amp.

        Args
        ----
        ampId : int
            The index of the desired amp. 0..7
        im : ndarray of image, optional
            If set, return the amp from the given image. Must have a geometry 
            compatible with self.
        leadingCols : bool
            If set, include the leadin columns
        leadingRows : bool
            If set, include the leadin rows
        """
        
        if im is None:
            im = self.image
            
        yr, xr = self.ampExtents(ampId, leadingRows=leadingRows, leadingCols=leadingCols)

        return im[yr, xr]

    def overscanCols(self, ampId, leadingRows=False, overscanRows=False):
        x0 = ampId*self.ncols + self.ampCols
        x1 = x0 + self.overCols

        xr = slice(x0, x1)
        yr = slice(self.leadinRows*(not leadingRows),
                   self.ccdRows + self.overRows*(overscanRows))

        return yr, xr

    def overscanRows(self, ampId, leadingCols=False, overscanCols=False):
        x0 = ampId*self.ncols + self.leadinCols*(not leadingCols)
        x1 = ampId*self.ncols + self.ampCols + self.overCols*(overscanCols)

        xr = slice(x0, x1)
        yr = slice(self.ccdRows, self.ccdRows + self.overRows)

        return yr, xr

    def overscanRowImage(self, ampId, im=None, leadingCols=False, overscanCols=False):
        if im is None:
            im = self.image
        yr, xr = self.overscanRows(ampId, leadingCols=leadingCols, overscanCols=overscanCols)
    
        return im[yr, xr]

    def overscanColImage(self, ampId, im=None, leadingRows=False, overscanRows=False):
        if im is None:
            im = self.image
        yr, xr = self.overscanCols(ampId, leadingRows=leadingRows, overscanRows=overscanRows)
    
        return im[yr, xr]

    def coreAmpImage(self, ampId, im=None, offset=(0,0)):
        if im is None:
            im = self.image
        ampImg = self.ampImage(ampId, im=im)

        amph, ampw = ampImg.shape
        hslice = slice(int(amph/2.0 - 200 + offset[0]),
                       int(amph/2.0 + 200 + offset[0]))
        wslice = slice(int(ampw/2.0 - 200 + offset[1]),
                       int(ampw/2.0 + 200 + offset[1]))
    
        return ampImg[hslice, wslice]

    def coreOverscanColImage(self, ampId, im=None):
        if im is None:
            im = self.image
        ampImg = self.overscanColImage(ampId, im=im)

        amph, ampw = ampImg.shape
        hslice = slice(int(amph/2.0 - 50), int(amph/2.0 + 50))
        wslice = slice(5, -1)
    
        return ampImg[hslice, wslice]

    def coreOverscanRowImage(self, ampId, im=None):
        if im is None:
            im = self.image
        ampImg = self.overscanRowImage(ampId, im=im)

        amph, ampw = ampImg.shape
        hslice = slice(5, -1)
        wslice = slice(int(ampw/2.0 - 50), int(ampw/2.0 + 50))
    
        return ampImg[hslice, wslice]

    def allAmpsImages(self, leadingCols=False, leadingRows=False):
        ampIms = []
        
        for i in range(self.namps):
            ampIms.append(self.ampImage(i, leadingCols=leadingCols,
                                        leadingRows=leadingRows))

        return ampIms
    
    def allOverscanColImages(self, leadingRows=False, overscanRows=False):
        osIms = []
        
        for i in range(self.namps):
            osIms.append(self.overscanColImage(i, leadingRows=leadingRows,
                                               overscanRows=overscanRows))

        return osIms
    
    def allOverscanRowImages(self, leadingCols=False, overscanCols=False):
        osIms = []
        
        for i in range(self.namps):
            osIms.append(self.overscanRowImage(i, leadingCols=leadingCols,
                                               overscanCols=overscanCols))

        return osIms

    def splitImage(self, doTrim=False, doFull=False):

        if doTrim and doFull:
            raise RuntimeError("cannot accept both doFull and doTrim")
        
        ampIms = []
        osColIms = []
        osRowIms = []

        for a_i in range(self.namps):
            if doTrim:
                ampIm = self.coreAmpImage(a_i)
                osColIm = self.coreOverscanColImage(a_i)
                osRowIm = self.coreOverscanRowImage(a_i)
            elif doFull:
                ampIm = self.ampImage(a_i, leadingCols=True, leadingRows=True)
                osColIm = self.overscanColImage(a_i, leadingRows=True)
                osRowIm = self.overscanRowImage(a_i, leadingCols=True, overscanCols=True)
            else:
                ampIm = self.ampImage(a_i)
                osColIm = self.overscanColImage(a_i)
                osRowIm = self.overscanRowImage(a_i)

            ampIms.append(ampIm)
            osColIms.append(osColIm)
            osRowIms.append(osRowIm)

        return ampIms, osColIms, osRowIms

    def replaceActiveFlux(self, newFlux, leadingRows=True):
        newImage = self.image.copy()

        for a_i in range(self.namps):
            yslice, xslice = self.ampExtents(a_i, leadingRows=False)
            inYslice, inXslice = self.finalAmpExtents(a_i, leadingRows=leadingRows)
            self.logger.debug("leadingRows=%s inYslice=%s", leadingRows, inYslice)
            newImage[yslice, xslice] = newFlux[inYslice,inXslice]

        return newImage
    
    def biasSubtract(self, bias=None):
        if bias is None:
            bias = Exposure(self.image, copyExposure=True, dtype='i4')
            biasParts = bias.splitImage()
            coreBiasParts = bias.splitImage(doTrim=True)
            for a_i in range(self.namps):
                m = np.median(coreBiasParts[1][a_i]) 
                coreBiasParts[1][a_i][:] = m
                m = np.median(biasParts[0][a_i]) 
                biasParts[0][a_i][:] = m
        else:
            biasParts = bias.splitImage()
            coreBiasParts = bias.splitImage(doTrim=True)
        
        parts = self.splitImage()
        coreParts = self.splitImage(doTrim=True)

        newAmps = []
        newOverCols = []
        
        for a_i in range(self.namps):
            biasOffset = int(np.median(coreBiasParts[1][a_i]) - np.median(coreParts[1][a_i]))

            newAmps.append(parts[0][a_i] - (biasParts[0][a_i] - biasOffset))
            newOverCols.append(parts[1][a_i] - (biasParts[1][a_i] - biasOffset))

        print(newAmps[0].shape, newOverCols[0].shape)
        return newAmps, newOverCols
        
    def biasSubtractOne(self, im=None, byRow=False):
        if im is None:
            im = self.image
        amps = []

        for a_i in range(8):
            amps.append(self.biasSubtractAmp(a_i, im=im, byRow=byRow))
        
        return np.hstack(amps)

    
    def biasSubtractAmp(self, ampId, im=None, byRow=False):
        if im is None:
            im = self.image

        ampIm = self.ampImage(ampId, im=im)
        osYr, osXr = self.overscanCols(ampId)

        if byRow:
            imMed = np.median(im[osYr, osXr],
                              axis=1, keepdims=True).astype('i4')
        else:
            osYr = slice(osYr.start + 500,
                         osYr.stop - 500)
            osXr = slice(osXr.start+3, osXr.stop)
        
            imMed = int(np.median(im[osYr, osXr]))

        return ampIm.astype('i4') - imMed
    
    def finalImage(self, leadingRows=False):
        ampImages = self.allAmpsImages(leadingRows=leadingRows)
        return np.hstack(ampImages)
    
def clippedStats(a, nsig=3.0, niter=20):
    a = a.reshape(-1)
    keep = np.ones(a.size, dtype=np.bool)

    for i in range(niter):
        mn = a[keep].mean()
        sd = a[keep].std()
        nkeep0 = keep.sum()
        keep &= np.abs(a - mn) < nsig*sd
        nkeep1 = keep.sum()

        if nkeep1 == nkeep0:
            return mn, sd, float(nkeep1)/a.size

    print("too many iterations (%d): fullsize=%d clipped=%d lastclipped=%d" % (i, a.size,
                                                                               float(nkeep0)/a.size,
                                                                               float(nkeep1)/a.size))
    return a[keep].mean(), a[keep].std(), float(nkeep1)/a.size

def clippedStack(flist, dtype='i4'):
    """ Return the median of a stack of images.

    Notes:
    
    The individual images are median-subtracted before medianed together, and 
    the median offset is added back to the result.

    All images must match the images type and exposure time of the first image.
    
    """

    exp = Exposure(flist[0], dtype=dtype)
    imshape = exp.image.shape
    imtype = exp.header['IMAGETYP']
    exptime = exp.header['EXPTIME']

    stack = np.zeros((len(flist), imshape[0], imshape[1]), dtype=dtype)
    meds = np.zeros(len(flist))

    for i, f in enumerate(flist):
        exp = Exposure(f)
        if exp.header['IMAGETYP'] != imtype or exp.header['EXPTIME'] != exptime:
            print("%s: unexpected imagetyp: %s(%0.2f) vs %s(%0.2f)"
                  % (f, 
                     exp.header['IMAGETYP'], exp.header['EXPTIME'],
                     imtype, exptime))
            
        meds[i] = np.median(exp.image)
        stack[i] = exp.image - int(meds[i])

    medimg = np.median(stack, axis=0)

    return Exposure(medimg)
                                                                                                            
def superBias(flist):
    biasParts = bias.splitImage(doTrim=False)
    coreBiasParts = bias.splitImage(doTrim=True)

    newAmps = []
    newOverCols = []
        
    for a_i in range(self.namps):
        biasOffset = int(np.median(coreBiasParts[1][a_i]) - np.median(coreParts[1][a_i]))
        
        newAmps.append(parts[0][a_i] - (biasParts[0][a_i] - biasOffset))
        newOverCols.append(parts[1][a_i] - (biasParts[1][a_i] - biasOffset))

    print(newAmps[0].shape, newOverCols[0].shape)
    return newAmps, newOverCols
        
def finalImage(exp, bias=None, dark=None, flat=None):
    exp = Exposure(exp)
    return exp.biasSubtract(bias)
    
def constructImage(ampIms, osColIms=None, osRowIms=None):
    orderedIms = []

    if osColIms is not None:
        for i in range(len(ampIms)):
            orderedIms.append(np.concatenate((ampIms[i], osColIms[i]), axis=1))
    else:
        orderedIms = ampIms[:]

    if osRowIms is not None:
        for i in range(len(orderedIms)):
            orderedIms[i] = np.concatenate((orderedIms[i], osRowIms[i]), axis=0)

    ret = np.concatenate(orderedIms, axis=1)
    return ret

def normAmpLevels(exp, fullCol=False):
    """ Using the overscan region, crudely normalize all amps to min=0. """

    retexp = Exposure(exp, copyExposure=True, dtype='f4')

    for i in range(retexp.namps):
        yr,xr = retexp.ampExtents(i, leadingCols=True, leadingRows=True, overscan=True)

        overscanImage = retexp.overscanColImage(i, leadingRows=True, overscanRows=True)

        if fullCol:
            overscanMed = np.mean(overscanImage, axis=1)[:, np.newaxis]
            retexp.image[:,xr] -= overscanMed
        else:
            overscanLevel = int(np.rint(np.median(overscanImage)))
            retexp.image[yr,xr] -= overscanLevel

    return retexp.image
