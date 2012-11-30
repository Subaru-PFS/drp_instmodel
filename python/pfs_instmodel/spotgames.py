import pyfits
import numpy
import scipy.ndimage
import scipy.signal
import matplotlib.pyplot as plt

import pfs_tools

def toyline(at, width, sampling=0.1, sigma=1.0, ongrid=False, donorm=False):
    x = numpy.arange(-width, width, step=sampling)
    y = numpy.exp(-0.5*(x/sigma)**2)
    if donorm:
        y /= sigma * numpy.sqrt(2*numpy.pi)
        
    return x+at, y

def readSpot(filename):
    spot = pyfits.getdata(filename)
    return spot

def imfreq(im):
    """ return the 2d frequency map for an image. """
    
    imfft = numpy.abs(numpy.fft.fftshift(numpy.fft.fftn(im)))
    return imfft

def imfreq1d(im):
    """ return the 1d frequency vector for an image, where the 2d frequencies are azimuthally flattened. """

    imfreqImg = imfreq(im)
    fmap = imfreqImg.flatten()
    dmap = distmap(im).flatten()

    sort_i = numpy.argsort(dmap)

    return dmap[sort_i], fmap[sort_i], imfreqImg

def distmap(arr, x0=None, y0=None):
    """ return the pixel distance map for a given array. """

    if x0 == None:
        x0 = arr.shape[1]/2
    if y0 == None:
        y0 = arr.shape[0]/2
    
    yd = numpy.linspace(0,arr.shape[0]-1,arr.shape[0]) - y0
    xd = numpy.linspace(0,arr.shape[1]-1,arr.shape[1]) - x0
    yyd,xxd = numpy.meshgrid(yd,xd)
    dmap = numpy.sqrt(yyd**2 + xxd**2)

    return dmap
    
def shiftSpot(spot, dx, dy, kernels=None, kargs=None):
    """ shift an image using seperable x&y lanczos kernels. """
    
    if not kargs:
        kargs = {}
    if kernels == None:
        xkernel = make1dKernel(offset=dx, **kargs)[1]
        ykernel = make1dKernel(offset=dy, **kargs)[1]
        kernels = (xkernel, ykernel)

    sspot = scipy.ndimage.convolve1d(spot, ykernel, axis=0)
    sspot = scipy.ndimage.convolve1d(sspot, xkernel, axis=1)
    return sspot, kernels
    
def shiftSpot2d(spot, dx, dy, kernels=None, kargs=None):
    """ shift an image using a 2d kernel built from seperable 1d lanczos kernels. """

    if not kargs:
        kargs = {}
    if kernels == None:
        xkernel = make1dKernel(offset=dx, **kargs)[1]
        ykernel = make1dKernel(offset=dy, **kargs)[1]
        kernel = numpy.outer(ykernel, xkernel)
        kernels = (xkernel, ykernel, kernel)

    sspot = scipy.signal.convolve(spot, kernel, mode='same')
    return sspot, kernels

def shiftSpotSpline(spot, dx, dy, kernels=None, kargs=None):
    """ shift an image using the scipy ndimage interpolation routines. """
    
    if not kargs:
        kargs = {}

    kernels = []
    sspot = scipy.ndimage.interpolation.shift(spot, (dy,dx), **kargs)
    return sspot, kernels
    
    
def lanczosWindow(x, n):
    if n > 0: 
        w = numpy.sinc(x/n)
        w[abs(x)>n] = 0
    else:
        w = x*0 + 1
    return w

def lanczos2(x):
    return lanczosWindow(x, 2)
    
def lanczos3(x):
    return lanczosWindow(x, 3)
    
def lanczos4(x):
    return lanczosWindow(x, 4)
    
def sincKernel(offset=0.0, padding=0, window=lanczosWindow, n=3, doNorm=True):
    if offset < -1 or offset > 1:
        raise ValueError("sinc offset must be in [-1,1], not %0.4f" % (offset))
        
    cnt = 2*(n+padding) + 1
    left = offset - n - padding
    right = offset + n + padding

    #print("n,offset=%d,%0.2f %d %0.2f %0.2f" % (n, offset,
    #                                            cnt,
    #                                            left, right))
        
    x = numpy.linspace(left, right, cnt)
    y = numpy.sinc(x)
    
    if window and n>0:
        w = window(x, n=n)
        y *= w
    
    if doNorm:
        y /= numpy.sum(y)
        
    return x, y
    
def make1dKernel(n=3, offset=0.0, padding=0, doNorm=True):
    """ Construct a centered 1-d lanczos-windowed sinc kernel. """

    x, y = sincKernel(n=n, offset=-offset, padding=padding)
    
    if doNorm:
        y /= numpy.sum(y)
        
    return x, y

def padArray(arr, padTo):
    assert arr.shape[0] == arr.shape[1]
    assert arr.shape[0] < padTo

    #print "old size=%d, new=%d" % (arr.shape[0], padTo)
    newSize = numpy.array([padTo, padTo])
    offset = (padTo-arr.shape[0])/2
    parr = numpy.zeros(newSize, dtype=arr.dtype)
    coreSlice = slice(offset, offset+arr.shape[0])
    parr[offset:offset+arr.shape[0],
         offset:offset+arr.shape[1]] = arr

    return parr, coreSlice

def unpadArray(arr, slice):
    return arr[slice, slice]

def applyPixelResponse(arr, pixelSize):
    kernelSize = pixelSize * 3
    kern = numpy.zeros((kernelSize, kernelSize), dtype=arr.dtype)
    kern[pixelSize:2*pixelSize, pixelSize:2*pixelSize] = 1
    kern /= numpy.sum(kern)
    
    out = scipy.signal.fftconvolve(arr, kern, mode='same')

    return out
    
def poo(arr, dx, dy, binFactor=10, padTo=0, applyPixelResp=False, kargs=None):
    assert dx>=0 and dy>=0

    # Trim raw image to multiples of binFactor pixels.
    maxSize = binFactor*(numpy.array(arr.shape,dtype='i4')/binFactor)
    newSize = (maxSize / binFactor).tolist()
    arr = arr[:maxSize[0],:maxSize[1]].copy()

    # Get our unshifted, binned, reference image.
    if applyPixelResp:
        arrs = applyPixelResponse(arr, 10)
        arr00 = pfs_tools.rebin(arrs, *newSize)
    else:
        arr00 = pfs_tools.rebin(arr, *newSize)

    if padTo:
        arr00unpadded = arr00.copy()
        arr00, padSlice = padArray(arr00, padTo)
        
    # Put the rest of this in loop....
    
    # interpolation-shift the binned reference image
    arrShifted, kernels = shiftSpot(arr00, float(dx)/binFactor, float(dy)/binFactor, kargs=kargs)

    # And pixel-shift (and optionally apply pixel response) then pixel-bin a comparison image
    arrPlaced = arr * 0
    arrPlaced[dy:,dx:] = arr[slice(None,-dy if dy else None),
                             slice(None,-dx if dx else None)]
    if applyPixelResp:
        arrPlaced = applyPixelResponse(arrPlaced, 10)
    arrPlaced = pfs_tools.rebin(arrPlaced, *newSize)

    if padTo:
        arrShifted = unpadArray(arrShifted, padSlice)
        arr00up = unpadArray(arr00, padSlice)
        assert numpy.all(numpy.equal(arr00unpadded, arr00up))
    
    return arr00, arrPlaced.copy(), arrShifted.copy(), kernels

def dispSpot(spotDict, key):
    bin, pad, shift = key
    res = spotDict[key]
    images = res[0]
    freqs = res[1]

    f = plt.figure()
    
def gatherPoo(spot, applyPixelResp=False, kargs=None):
    tries = ([1,0],
             [2,0],
             [3,0],
             [4,0],
             [4,128],
             [5,0],
             [6,0],
             [7,0],
             [8,0],
             [8,128],
             [9,0],
             [10,128],
             [10,0])

    if not kargs:
        kargs = dict(n=3)
        
    all = {}
    for bin,pad in tries:
        for shift in range(bin):
            print "processing %s with kargs=%s" % ((bin,pad,shift), kargs)
            ret = poo(spot, shift, 0, 
                      applyPixelResp=applyPixelResp,
                      binFactor=bin, 
                      padTo=pad,
                      kargs=kargs)
            freq = imfreq1d(ret[2])
            all[(bin, pad, shift)] = ret, freq

    return all, kargs
        
            
