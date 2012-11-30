import pyfits
import numpy
import scipy.ndimage
import scipy.signal

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
        kernels = (kernel, xkernel, ykernel)

    sspot = scipy.signal.convolve(spot, kernel, mode='same')
    return sspot, kernels
    
def poo(arr, dx, dy, binFactor=10, kargs=None):
    assert dx>=0 and dy>=0
    
    # Trim raw image to multiples of binFactor pixels.
    maxSize = binFactor*(numpy.array(arr.shape,dtype='i4')/binFactor)
    newSize = (maxSize / binFactor).tolist()
    arr = arr[:maxSize[0],:maxSize[1]].copy()

    # Get our unshifted, binned, reference image.
    arr00 = pfs_tools.rebin(arr, *newSize)

    # Put the rest of this in loop....
    
    # interpolation-shift the binned reference image
    arrShifted, kernels = shiftSpot2d(arr00, float(dx)/binFactor, float(dy)/binFactor, kargs=kargs)

    # And pixel-shift then pixel-bin a comparison image
    arrPlaced = arr * 0
    arrPlaced[dy:,dx:] = arr[slice(None,-dy if dy else None),
                             slice(None,-dx if dx else None)]
    arrPlaced = pfs_tools.rebin(arrPlaced, *newSize)
    
    return arr00, arrPlaced, arrShifted, kernels
    
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
    
def sincKernel(offset=0.0, n=3, padding=0, window=lanczos3, doNorm=True):
    if offset < -1 or offset > 1:
        raise ValueError("sinc offset must be in [-1,1], not %0.4f" % (offset))
        
    cnt = 2*(n+padding) + 1
    left = offset - n - padding
    right = offset + n + padding

    print("n,offset=%d,%0.2f %d %0.2f %0.2f" % (n, offset,
                                                cnt,
                                                left, right))
        
    x = numpy.linspace(left, right, cnt)
    y = numpy.sinc(x)
    
    if window:
        w = window(x)
        y *= w
    
    if doNorm:
        y /= numpy.sum(y)
        
    return x, y
    
def make1dKernel(n=3, offset=0.0, padding=0, doNorm=True):
    """ Construct a centered 1-d lanczos-windowed sinc kernel. """

    x, y = sincKernel(n=n, offset=-offset, padding=padding)
    w = lanczosWindow(x, n=n)
    y *= w
    
    if doNorm:
        y /= numpy.sum(y)
        
    return x, y
