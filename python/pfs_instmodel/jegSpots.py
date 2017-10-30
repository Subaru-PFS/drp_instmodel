from __future__ import division
from __future__ import print_function
from __future__ import absolute_import

from builtins import str
from builtins import range
from past.builtins import basestring
import glob
import gzip
import logging
import os
import re
import sys
import time

import numpy
import scipy.signal
import scipy.ndimage

import astropy.io.fits as fits

from . import spotgames

jegLogger = logging.getLogger('jegSpots')

from pfs_tools import pydebug
"""
Read the raw PSF spot images from JEG's optical design. 

The raw file 
OFFSETS :
    HEADER: 0
    DESIGN FILE: 2K
    X,Y,FOC (NIMAGE*3 floats): 8K
    DATA (NIMAGE*XSIZE*YSIZE u-shorts) 32K

"""

def getDataPath(date='2016-10-26', band='Red', frd=23, focus=0, slitFocus=0, fieldAngle=0, spotDir=None):
    """ Return complete data path for the specified spots. 

    The file format changed on 2017-10-30 from a homebrew binary to FITS.
    """
    if not spotDir:
        spotDir = os.path.join(os.environ['DRP_INSTDATA_DIR'], 'data/spots/jeg')

    if date > '2017-10-29':
        spotFile = os.path.join(spotDir, date, band, 
                                "PFSsim*_f%3d_*.fits" % (focus))
    elif date > '2016-10-01':
        spotFile = os.path.join(spotDir, date, band, 
                                "*.dat_foc%d_frd%d_sfld%02d.imgstk" % (focus, frd, fieldAngle))
    else:
        spotFile = os.path.join(spotDir, date, band, 
                                "*.dat_foc%d_frd%d.imgstk" % (focus, frd))

    files = glob.glob(spotFile) + glob.glob(spotFile + '.gz')
    if len(files) != 1:
        raise RuntimeError("There is not a unique JEG spotfile matching %s{.gz}; found %d" % 
                           (spotFile, len(files)))

    return files[0]

def resolveSpotPathSpec(pathSpec):
    if isinstance(pathSpec, basestring):
        return pathSpec
    if pathSpec is None:
        return getDataPath()
    
    return getDataPath(**pathSpec)

def makeFiberImage(fiberRadius=28, shape=(64,64), dtype='f4'):
    """ Return the image we convolve the spots with. """
    
    im = numpy.zeros(shape, dtype=dtype)
    spotRadius = shape[0]//2

    x,y = numpy.meshgrid(numpy.arange(-spotRadius,spotRadius+1),
                         numpy.arange(-spotRadius,spotRadius+1))
    d = numpy.sqrt(x**2 + y**2)

    im[numpy.where(d <= fiberRadius)] = 1

    return im
    
def convolveWithFiber(spot, fiberImage=None):
    if fiberImage is None:
        fiberImage = makeFiberImage()
    
    return scipy.ndimage.convolve(spot, fiberImage, mode='constant', cval=0.0)

def getFiberIds(header, useArrayKeys):
    """ Given the header, return the fiberIDs. We assume that the spacing is proportional to the SINANGs. """

    nang = header['NANG']
    maxFiber = header['MAXFIBER']
    if useArrayKeys is False:
        angs = []
        for i in range(nang):
            ang = float(header['SINANG[%d]' % (i)])
            angs.append(ang)
    elif useArrayKeys is True:
        angs = header['SINANG[]']
        assert len(angs) == nang
    else:
        raise ValueError("useArrayKeys must be True or False.")
    
    angs = numpy.array(angs)
    normAngs = angs / angs.max()
    fibers = (normAngs * (maxFiber-1)).astype('i2')

    jegLogger.info("fiber angles: %s" % (angs))
    jegLogger.info("fiber IDs   : %s" % (fibers))
    
    return fibers

_spotCache = dict()
def clearSpotCache():
    global _spotCache
    _spotCache.clear()


def readSpotFile(pathSpec, doConvolve=None, doRebin=False, 
                 doNorm=True, 
                 doSwapaxes=True, doTrimSpots=True,
                 doRecenter=None,
                 verbose=False, clearCache=False):
    """ Directly read some recent version of JEGs spots.

    Args
    ----
    pathSpec : dict
       Enough info to uniquely identify a jeg spot dataset.
       The keys are: (date, band, frd, focus). 
    doConvolve : bool or None
       Whether to convolve spots with a circular fiber tophat.
       By default, set iff spot data version == 1.
    doNorm : bool
       Whether to normalize flux to the brightest spot's total.
    doSwapaxes: bool
       Whether to swap X&Y, so that the traces and CCD columns 
       run vertically.
    doTrimSpots : bool
       Whether to trim 0-pixel borders from the input spots.
    doRecenter : bool
       Whether to "fix" the spot centers to the measured moments.
       By default, set iff spot data version == 2.
    doRebin : bool/integer
       If not False, attempt to rebin pixels by the given factor.

    Returns
    -------
    arr : the table of spots, with wavelength(AA), fiberID, xc, yc, image
    metadata : the partially converted header from the .imgstk file.

    The .imgstk file contains a few 1k-aligned sections, described by
    a commented-out header section. For version 1:
    
    OFFSETS :
        HEADER: 0
        DESIGN FILE: 2K
        X,Y,FOC (NIMAGE*3 floats): 8K
        DATA (NIMAGE*XSIZE*YSIZE u-shorts) 32K
        
    """

    path = resolveSpotPathSpec(pathSpec)
    if clearCache:
        clearSpotCache()

    if path in _spotCache:
        return _spotCache[path]
    
    if path.endswith('.gz'):
        fopen = gzip.open
    else:
        fopen = open

    with fopen(path, 'r') as f:
        rawHeader = f.read(2*1024)

    rawHeader = rawHeader.rstrip('\0\n')
    header = rawHeader.split('\n')

    def k2Decimal(s):
        """ Convert '8K' -> 8192 """
        if s[-1] in {'K', 'k'}:
            return int(s[:-1]) * 1024
        else:
            return int(s)
        
    headerDict = {'TITLE':header.pop(0)}
    converters = {'SLIT_RADIUS':float,
                  'SLIT-VPUPIL_DISTANCE':float,
                  'FOFFSET':float,
                  'FSLOPE':float,
                  'FRD_SIGMA':float,
                  'FOFFSET':float,
                  'XPIX':float,
                  'YPIX':float,
                  'IMPEAK':int,
                  
                  'FRD_ON':int,

                  'NIMAGE':int,
                  'NLAM':int,
                  'NANG':int,
                  'XSIZE':int,
                  'YSIZE':int,

                  'HDROFFS':k2Decimal,
                  'DESOFFS':k2Decimal,
                  'SBOFFS':k2Decimal,
                  'XYFOFFS':k2Decimal,
                  'DATOFFS':k2Decimal,
                  
                  'LAM[]': lambda _s : [float(x) for x in _s.split(' ')],
                  'SINANG[]': lambda _s : [float(x) for x in _s.split(' ')]
    }
    
    for i, h in enumerate(header):
        m = re.match('(^[A-Z][-A-Z_ 0-9\[\]]*) = (.*)', h)
        if m:
            key, value = m.groups()
            convert = converters.get(key, str)
            headerDict[key] = convert(value.strip())
        else:
            headerDict['LINE%03d'%(i)] = h

    if 'LAM[1]' in headerDict:
        dataVersion = 1
    elif 'SBOFFS' in headerDict:
        dataVersion = 3
    else:
        dataVersion = 2
        
    headerDict['DATA_VERSION'] = dataVersion
    headerDict['FILENAME'] = os.path.basename(path)
    
    with fopen(path, 'r') as f:
        _ = f.read(headerDict['DESOFFS'])
        if dataVersion < 3:
            rawDesign = f.read(headerDict['XYFOFFS']-headerDict['DESOFFS'])    # unused, so far.
            rawPositions = f.read(headerDict['DATOFFS']-headerDict['XYFOFFS'])
            rawData = f.read()
        else:
            rawDesign = f.read(headerDict['SBOFFS']-headerDict['DESOFFS'])     # unused, so far.
            rawBrightness = f.read(headerDict['XYFOFFS']-headerDict['SBOFFS']) # unused, so far.
            rawPositions = f.read(headerDict['DATOFFS']-headerDict['XYFOFFS'])
            rawData = f.read()

    # Add in a temporary fiber range, until the .imgstk file specifies one
    headerDict['MAXFIBER'] = 325
    jegLogger.debug("  XXX: hardwired max(fiberid)==%d logic" % (headerDict['MAXFIBER']))

    if doConvolve is None:
        doConvolve = dataVersion == 1
    if doRecenter is None:
        doRecenter = dataVersion == 2
        
    nimage = headerDict['NIMAGE']
    xsize = headerDict['XSIZE']
    ysize = headerDict['YSIZE']
    positions = numpy.fromstring(rawPositions, dtype='3f4', count=nimage)

    wavelengths = []
    nlam = headerDict['NLAM']
    useArrayKeys = dataVersion >= 2
    if not useArrayKeys:
        for i in range(nlam):
            waveKey = "LAM[%d]" % (i)
            wavelength = float(headerDict[waveKey])
            wavelengths.append(wavelength)
    else:
        wavelengths = headerDict['LAM[]']
        assert len(wavelengths) == nlam
                        
    fiberIDs = getFiberIds(headerDict, useArrayKeys)
    nfibers = len(fiberIDs)
    
    if dataVersion < 3:
        data = numpy.fromstring(rawData, dtype='(%d,%d)u2' % (xsize,ysize), count=nimage).astype('f4')
    elif dataVersion == 3:
        data = numpy.fromstring(rawData, dtype='(%d,%d)f4' % (xsize,ysize), count=nimage).astype('f4')
    else:
        raise ValueError("unknown spot file version: %s" % (dataVersion))
    jegLogger.info("raw spot version %d data type %s, range: %g..%g",
                   dataVersion, data.dtype, data.min(), data.max())

    rawspots = data.copy()
    if doRebin is not False:
        newSize = xsize//doRebin
        if newSize*doRebin != xsize:
            raise ValueError('doRebin must evenly divide the raw spot size')

        newData = numpy.empty(shape=(data.shape[0], newSize, newSize), dtype=data.dtype)
        for i in range(data.shape[0]):
            tspot = spotgames.rebinBy(data[i], doRebin)
            newData[i,:,:] = tspot
            
        data = newData
        xsize = ysize = newSize
        headerDict['XPIX'] *= doRebin
        headerDict['YPIX'] *= doRebin
        
    if doTrimSpots:
        if xsize % 2 == 0 and dataVersion < 3:
            jegLogger.info("trimming outer pixel of %s spots" % (xsize))
            data = data[:,:-1,:-1]
        extents = dataWidth(data)
        jegLogger.info("spot extents: %s" % (str(extents)))
        
        if not isinstance(doTrimSpots, bool):
            tryFor = doTrimSpots
        else:
            tryFor = None
            
        data, trimmedPixels = trimSpots(data, tryFor=tryFor, leaveBorder=3)
        if doTrimSpots > 1:
            assert data.shape[-1] == tryFor, ("found trimmed spots to be %d pixels, wanted %d" %
                                              (data.shape[-1], tryFor))
        jegLogger.info("trimmed spots by %d pixels to %s pixels (%g mm)" % (2*trimmedPixels, data.shape[-1],
                                                                            data.shape[0] * headerDict['XPIX']))

    xsize, ysize = data.shape[1:]
    headerDict['XSIZE'] = xsize
    headerDict['YSIZE'] = ysize

    if doConvolve:
        jegLogger.info("convolving with fiber image: %s" % (doConvolve))
        fiberImage = makeFiberImage()

    spots = []
    symSpots = []
    maxFlux = max([data[d].sum() for d in range(data.shape[0])])
    jegLogger.info("max spot = %0.2f" % (maxFlux))

    for i in range(data.shape[0]):
        fiberIdx = fiberIDs[i // nlam]
        wavelength = wavelengths[i % nlam]
        xc = positions[i,0]
        yc = -positions[i,1]
        focus = positions[i,2]

        if doConvolve:
            spot = convolveWithFiber(data[i,:,:], fiberImage)
        else:
            spot = data[i,:,:]

        assert spot.shape[0] == spot.shape[1]
        spotw = spot.shape[0]//2
        ctr0 = spotgames.centroid(spot)
        #assert numpy.abs(ctr0[0] - spotw) < 0.1, "centroid Y too far from center (%g %g)" % (ctr0[0], spotw)
        #assert numpy.abs(ctr0[1] - spotw) < 0.1, "centroid X too far from center (%g %g)" % (ctr0[1], spotw)

        if doRecenter:
            pspot, spotSlice = spotgames.padArray(spot, padTo=spot.shape[0]*2)
            pspot, _ = spotgames.shiftSpot1d(pspot, spotw-ctr0[0], spotw-ctr0[1], kargs=dict(n=spotw*3//2))
            spot = pspot[spotSlice, spotSlice]
            ctr2 = spotgames.centroid(spot)

            fluxMin = spot.min()
            if fluxMin < 0:
                spot -= fluxMin
                
            jegLogger.debug("recentered spot %i by (%0.5f, %0.5f) to (%0.5f, %0.5f) [min=%0.5f]" %
                            (i,
                             spotw-ctr0[0], spotw-ctr0[1],
                             ctr2[0], ctr2[1],
                             fluxMin))

        # Rotate x-up mechanical view to y-up detector view (dispersing along columns)
        if doSwapaxes:
            spot = numpy.swapaxes(spot,0,1)
            xc, yc = yc, xc

            rawSpot = numpy.swapaxes(rawspots[i,:,:],0,1)
            
        if doNorm:
            if doNorm == 'peak':
                # Normalize the flux of each spot w.r.t. this spot peak.
                spot /= spot.max()
            else:
                # Normalize the flux of each spot, w.r.t. the brightest spot.
                spot /= maxFlux
            jegLogger.debug("normalized to %g..%g sum=%g...." % (spot.min(), spot.max(), spot.sum()))
        
        ctr = spotgames.centroid(spot)
        dCtr = numpy.abs(ctr - numpy.array(spot.shape).T/2.0)
        if numpy.any(dCtr > 0.01):
            jegLogger.warn("spot  %d (%d, %0.2f) at (%0.2f %0.2f %0.3f) (%0.4f,%0.4f) (%0.4f, %0.4f), sum=%0.3f" % 
                           (i, fiberIdx, wavelength, xc, yc, focus,
                            ctr[0], ctr[1], dCtr[0], dCtr[1],
                            spot.sum()))

        spots.append((fiberIdx, wavelength, xc, yc, focus, rawSpot, spot, ctr))
            
        if fiberIdx != 0:
            flipCtr = ctr.copy()
            flipCtr[0] = spot.shape[1]/2 - (flipCtr[0] - spot.shape[1]/2)
            symSpots.append((-fiberIdx, wavelength, -xc, yc, focus, rawSpot[:,::-1], spot[:,::-1], flipCtr))

    allSpots = spots + symSpots    
    spotw = spots[0][-2].shape[0]
    spotDtype = numpy.dtype([('fiberIdx','i2'),
                             ('wavelength','f8'),
                             ('spot_xc','f8'), ('spot_yc','f8'), ('spot_focus','f4'),
                             ('rawspot', '(%d,%d)%s' % (rawSpot.shape[0],
                                                        rawSpot.shape[1],
                                                        rawSpot.dtype)),
                             ('spot', '(%d,%d)f4' % (spotw,spotw)),
                             ('ctr', '2f4')])

    tarr = numpy.array(allSpots, dtype=spotDtype)

    # Now clean up... later steps expect to have fiber IDs and wavelengths in order
    arr = numpy.sort(tarr, order=('fiberIdx','wavelength'))
    assert(numpy.all(arr['spot'] >= 0)), "spots must be non-negative."
    
    _spotCache[path] = (arr, headerDict)

    return arr, headerDict

def dataWidth(spots):
    """ Return the x/y radius of the furthest nonzero pixel. """
    
    assert spots.shape[-2] == spots.shape[-1], "input spots must be square"
    
    startWidth = spots.shape[-1]
    nz = numpy.where(spots > 0)
    mn = min(nz[1].min(), nz[2].min())
    mx = max(nz[1].max(), nz[2].max())

    mxt = startWidth - mx
    return max(mn, mxt), min(mn, mxt)

def borderWidth(spots):
    """ Return the number of all-zero pixels around the edge of the image. """
    
    assert spots.shape[-2] == spots.shape[-1], "input spots must be square"
    
    startWidth = spots.shape[-1]
    nz = numpy.where(spots > 0)
    mn = min(nz[1].min(), nz[2].min())
    mx = max(nz[1].max(), nz[2].max())
    mxt = startWidth - mx - 1

    return min(mn, mxt)

def trimSpots(spots, tryFor=None, leaveBorder=1):
    """ Return a spot array with the outer 0-level pixels trimmed off. """

    startWidth = spots.shape[-1]
    trimPix = borderWidth(spots)
    trimPix -= leaveBorder

    if trimPix <= 0:
        raise RuntimeError("cannot trim spots by %d pixels" % (trimPix))
        
    if tryFor is not None:
        if tryFor < startWidth - 2*trimPix:
            raise RuntimeError("cannot safely trim spots to %s pixels" % (tryFor))
        trimPix = (startWidth-tryFor)//2

    jegLogger.debug("trimming border by %d pixels; image goes from %d to %d pixels" %
                    (trimPix, startWidth, startWidth-trimPix*2)) 
    nspots = spots[:,trimPix:startWidth-trimPix,trimPix:startWidth-trimPix]

    return nspots, trimPix

def oversampleSpots(spots, factor):
    """ Oversample the given spots by factor. """

    if factor == 1:
        return spots

    assert spots.shape[1] == spots.shape[2]

    t0 = time.time()
    print("oversampling %s (%d bytes) by %d" % (spots.shape, 
                                                spots.nbytes,
                                                factor))

    newWidth = spots.shape[1]*factor
    newSpots = numpy.zeros(shape=(spots.shape[0], newWidth, newWidth),
                           dtype=spots.dtype)

    ix1 = (numpy.arange(newWidth,dtype='f4')//factor).astype('i4')
    ix = numpy.tile(ix1, (newWidth,1))
    ixy = (ix, ix.T)

    for i in range(spots.shape[0]):
        newSpots[i,:,:] = spots[i][ixy] / (factor*factor)

    t1 = time.time()
    print("oversampled to %s (%d bytes) in %0.4fs" % (newSpots.shape, 
                                                      newSpots.nbytes,
                                                      t1-t0))

    return newSpots

def writeSpotFITS(pathSpec, outDir, data, headerDict):
#    raise NotImplementedError("writeSpotFITS() no longer needed or tested")

    phdu = fits.PrimaryHDU()
    phdr = phdu.header
    phdr.update('pixscale', 0.003, 'mm/pixel')

    cols = []
    cols.append(fits.Column(name='fiberIdx',
                            format='I',
                            array=data['fiberIdx']))
    cols.append(fits.Column(name='wavelength',
                            format='D',
                            array=data['wavelength']))
    cols.append(fits.Column(name='spot_xc',
                            format='D',
                            array=data['spot_xc']))
    cols.append(fits.Column(name='spot_yc',
                            format='D',
                            array=data['spot_yc']))
    spots = data['spot'][:]
    spotw = spots[0].shape[0]
    spots.shape = (len(spots), spotw*spotw)
    cols.append(fits.Column(name='spot',
                            format='%dE' % (spotw*spotw),
                            dim='(%d,%d)' % (spotw, spotw),
                            array=spots))
    colDefs = fits.ColDefs(cols)

    thdu = fits.new_table(colDefs)
    hdulist = fits.HDUList([phdu, thdu])

    filename = "%(date)s_%(band)s_%(frd)02d_%(focus)+02d.fits" % pathSpec
    
    hdulist.writeto(os.path.join(outDir, filename), 
                    checksum=True, clobber=True)

    
def main(band, spotDir=None, filePattern=None):
    """ Convert a directory of zemax spot files into a slightly more convenient FITS table. """

    if not spotDir:
        spotDir = os.path.join(os.environ['DRP_INSTDATA_DIR'], 'data/spots/jeg', band)

    data = readSpotDir(spotDir, filePattern)
    writeSpotFITS(spotDir, data)

if __name__ == "__main__":
    main(sys.argv[1])
