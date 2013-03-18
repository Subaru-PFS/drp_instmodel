import glob
import os
import re
import sys

import numpy
import scipy.signal
import scipy.ndimage

import pyfits

import splinedPsf
import pfs_tools

from pfs_tools import pydebug
"""
Read the raw PSF spot images from JEG's optical design, and convert
them to some slightly more efficient intermediate FITS file..

The raw file 
OFFSETS :
    HEADER: 0
    DESIGN FILE: 2K
    X,Y,FOC (NIMAGE*3 floats): 8K
    DATA (NIMAGE*XSIZE*YSIZE u-shorts) 32K

"""

def displaySpot(spot, scale=15.0):
    import xplot

def makeFiberImage(fiberRadius=28, shape=(64,64), dtype='f4'):
    """ Return the image we convolve the spots with. """
    
    im = numpy.zeros(shape, dtype=dtype)
    spotRadius = shape[0]/2

    x,y = numpy.meshgrid(numpy.arange(-spotRadius,spotRadius+1),
                         numpy.arange(-spotRadius,spotRadius+1))
    d = numpy.sqrt(x**2 + y**2)

    im[numpy.where(d <= fiberRadius)] = 1

    return im
    
def convolveWithFiber(spot, fiberImage=None):
    if fiberImage == None:
        fiberImage = makeFiberImage()
    
    return scipy.ndimage.convolve(spot, fiberImage, mode='constant', cval=0.0)

def getFiberIds(header, headerVersion):
    """ Given the header, return the fiberIDs. We assume that the spacing is proportional to the SINANGs. """

    nang = int(header['NANG'])
    if headerVersion == 1:
        angs = []
        for i in range(nang):
            ang = float(header['SINANG[%d]' % (i)])
            angs.append(ang)
    else:
        angs = [float(ang) for ang in header['SINANG[]'].split()]
        assert len(angs) == nang
                        
    angs = numpy.array(angs)
    normAngs = angs / angs.max()
    fibers = (normAngs * 300).astype('i2')

    print "fiber angles: %s" % (angs)
    print "fiber IDs   : %s" % (fibers)
    
    return fibers
    
def readSpotDir(path, doConvolve=None, filePattern='*.imgstk'):
    """ Read slightly cleaned up versions of JEGs spots.

    The .imgstk file contains a few 1k-aligned sections:
    
    OFFSETS :
        HEADER: 0
        DESIGN FILE: 2K
        X,Y,FOC (NIMAGE*3 floats): 8K
        DATA (NIMAGE*XSIZE*YSIZE u-shorts) 32K
        
    """

    names = glob.glob(os.path.join(path, '*.imgstk'))
    assert len(names) == 1, "path %s does not have exactly one .imgstk file" % (path)
    name = names[0]
    if name.endswith('.gz'):
        fopen = gzip.open
    else:
        fopen = open

    with fopen(names[0], 'r') as f:
        rawHeader = f.read(2*1024)
        rawDesign = f.read((8-2)*1024)    # unused, so far.
        rawPositions = f.read((32-8)*1024)
        rawData = f.read()

    rawHeader = rawHeader.rstrip('\0\n')
    header = rawHeader.split('\n')

    headerDict = {'TITLE':header.pop(0)}
    for i, h in enumerate(header):
        m = re.match('(^[A-Z][-A-Z_ 0-9\[\]]*) = (.*)', h)
        if m:
            key, value = m.groups()
            headerDict[key] = value
        else:
            headerDict['LINE%03d'%(i)] = h

    headerVersion = 1 if 'LAM[1]' in headerDict else 2
    if doConvolve == None:
        doConvolve = headerVersion == 1
        
    nimage = int(headerDict['NIMAGE'])
    xsize = int(headerDict['XSIZE'])
    ysize = int(headerDict['YSIZE'])
    positions = numpy.fromstring(rawPositions, dtype='3f4', count=nimage)
    
    wavelengths = []
    nlam = int(headerDict['NLAM'])
    if headerVersion == 1:
        for i in range(nlam):
            waveKey = "LAM[%d]" % (i)
            wavelength = float(headerDict[waveKey])
            wavelengths.append(wavelength)
    else:
        wavelengths = [float(lam) for lam in headerDict['LAM[]'].split()]
        assert len(wavelengths) == nlam
                        
    fiberIDs = getFiberIds(headerDict, headerVersion)
        
    data = numpy.fromstring(rawData, dtype='(%d,%d)u2' % (xsize,ysize), count=nimage).astype('f4')
    fiberImage = makeFiberImage()

    spots = []
    symSpots = []
    for i in range(data.shape[0]):
        fiberIdx = fiberIDs[i / nlam]
        wavelength = wavelengths[i % nlam]
        xc = positions[i,0]
        yc = positions[i,1]
        if doConvolve:
            spot = convolveWithFiber(data[i,:,:], fiberImage)
        else:
            spot = data[i,:,:]

        # bin from 256x256 1um pixels to 85x85 3um pixels.
        # Still oversampled by 5.
        spot = spot[:-1,:-1]
        spot = pfs_tools.rebin(spot, 85,85)

        # Require total flux to be 1.0
        spot /= spot.sum()
        
        # Rotate x-up mechanical view to y-up detector view (dispersing along columns)
        spot = numpy.swapaxes(spot,0,1)
        xc, yc = yc, xc
        spots.append((fiberIdx, wavelength, xc, yc, spot))
        print("spot  %d (%d, %0.2f) at (%0.1f %0.1f)" % (i, fiberIdx, wavelength, xc, yc))
        if fiberIdx != 0:
            symSpots.append((-fiberIdx, wavelength, -xc, yc, spot[::-1,:]))

    allSpots = spots + symSpots    
    spotw = spots[0][-1].shape[0]
    spotDtype = numpy.dtype([('fiberIdx','i2'),
                             ('wavelength','f4'),
                             ('spot_xc','f4'), ('spot_yc','f4'),
                             ('spot', '(%d,%d)f4' % (spotw,spotw))])

    tarr = numpy.array(allSpots, dtype=spotDtype)

    # Now clean up... later steps expect to have fiber IDs and wavelengths in order
    arr = numpy.sort(tarr, order=('fiberIdx','wavelength'))

    return arr

def writeSpotFITS(spotDir, data):

    phdu = pyfits.PrimaryHDU()
    phdr = phdu.header
    phdr.update('pixscale', 0.003, 'mm/pixel')

    cols = []
    cols.append(pyfits.Column(name='fiberIdx',
                              format='I',
                              array=data['fiberIdx']))
    cols.append(pyfits.Column(name='wavelength',
                              format='D',
                              array=data['wavelength']))
    cols.append(pyfits.Column(name='spot_xc',
                              format='D',
                              array=data['spot_xc']))
    cols.append(pyfits.Column(name='spot_yc',
                              format='D',
                              array=data['spot_yc']))
    spots = data['spot'][:]
    spotw = spots[0].shape[0]
    spots.shape = (len(spots), spotw*spotw)
    cols.append(pyfits.Column(name='spot',
                              format='%dE' % (spotw*spotw),
                              dim='(%d,%d)' % (spotw, spotw),
                              array=spots))
    colDefs = pyfits.ColDefs(cols)

    thdu = pyfits.new_table(colDefs)
    hdulist = pyfits.HDUList([phdu, thdu])

    hdulist.writeto(os.path.join(spotDir, 'spots.fits'), 
                    checksum=True, clobber=True)

    
def main(band, spotDir=None, filePattern=None):
    """ Convert a directory of zemax spot files into a slightly more convenient FITS table. """

    if not spotDir:
        spotDir = os.path.join(os.environ['PFS_INSTDATA_DIR'], 'data/spots/jeg', band)

    data = readSpotDir(spotDir, filePattern)
    writeSpotFITS(spotDir, data)

if __name__ == "__main__":
    main(sys.argv[1])
