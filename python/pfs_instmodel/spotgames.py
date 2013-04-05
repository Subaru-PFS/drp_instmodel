import pyfits
import numpy
import scipy.ndimage
import scipy.signal
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

import pfs_tools
import plotutils

import pfs_instmodel.jegSpots as jegSpots

def radToIndices(rad, sampling=1.0, offset=0.0):
    """ Return the indices for an odd-width vector. """
    
    if rad < 1:
        raise RuntimeError("radius must be positive.")
    x = numpy.linspace(-(rad-1), rad-1, 
                       (2*rad-1)/sampling).astype('f4')
    return x+offset

def radToImageIndices(rad, sampling=1.0, offset=None):
    """ Return the indices for an odd-width image. """
    
    x = radToIndices(rad, sampling=sampling)
    xx, yy = numpy.meshgrid(x, x)

    if offset != None:
        xx += offset[0]
        yy += offset[1]
        
    r = numpy.sqrt(xx**2 + yy**2)

    return r
    
def radFuncToVector(func, rad, sampling=1.0):
    """ Take a symmetric function and expand it into an array. 
    0 is placed at the center of an odd-width vector.
    """
    
    x = radToIndices(rad, sampling=sampling)
    return x, func(x)
    
def radFuncToImage(func, rad, sampling=1.0):
    """ Take a circularily symmetric function and expand it into an image. 
    0,0 is placed at the center of an odd-width image.
    """

    r = radToImageIndices(rad, sampling)
    im = func(r)

    return r, im


def tophat(rad, width):
    return (rad < width).astype('f4')
        
def gaussian(rad, at=0.0, sigma=1.0, donorm=False):
    y = numpy.exp(-0.5*(rad/sigma)**2)
    if donorm:
        y /= sigma * numpy.sqrt(2*numpy.pi)
        
    return rad+at, y

def toyspot(at, width, sampling=0.1, sigma=1.0, ongrid=False, donorm=False):
    x = numpy.arange(-width, width, step=sampling)
    y = numpy.exp(-0.5*(x/sigma)**2)
    if donorm:
        y /= sigma * numpy.sqrt(2*numpy.pi)
        
    return x+at, y

def readSpot(filename):
    spot = pyfits.getdata(filename)
    return spot

def imfreq(im, doShift=True, doAbs=False):
    """ return the 2d frequency map for an image. """
    
    fftim = numpy.fft.fftn(im)
    if doAbs:
        fftim = numpy.abs(fftim)
    if doShift:
        fftim = numpy.fft.fftshift(fftim)

    return fftim


def imfreq1d(im, doAbs=False, doAverage=True, sampling=None):
    """ return the 1d frequency vector for an image, where the 2d frequencies are azimuthally flattened. """

    imfreqImg = imfreq(im, doAbs=doAbs)
    fmap = imfreqImg.flatten()
    dmap = distmap(im).flatten()

    sort_i = numpy.argsort(dmap)

    freqs = dmap[sort_i]
    ampl = fmap[sort_i] 
    if doAverage:
        uux, uui = numpy.unique(freqs, return_index=True)
        for i in range(len(uui)-1):
            ampl[uui[i]:uui[i+1]] = numpy.mean(numpy.abs(ampl[uui[i]:uui[i+1]]))

    if sampling:
        freqs = freqs / (numpy.sqrt(2) * im.shape[0]*sampling)
    return freqs, ampl, imfreqImg

def freq2im(fftim, doAbs=True, doShift=True):
    """ return the image from an fft. """
    
    if doShift:
        fftim = numpy.fft.ifftshift(fftim)

    im = numpy.fft.ifftn(fftim)
    if doAbs:
        im = numpy.abs(im)

    return im

def imextent(im, scale=1.0, doCenter=False):
    """ return (l,r,b,t) indices for the pixel centers. """

    rows,cols = im.shape
    rows /= scale
    cols /= scale
    x0 = 0 if doCenter==False else -cols/2.0
    y0 = 0 if doCenter==False else -rows/2.0
    
    return (x0, x0+cols, y0, y0+rows)

def distmap(arr, x0=None, y0=None):
    """ return the pixel distance map for a given array. """

    if x0 == None:
        x0 = (arr.shape[1]-1)/2.0
    if y0 == None:
        y0 = (arr.shape[0]-1)/2.0
    
    yd = numpy.linspace(0,arr.shape[0]-1,arr.shape[0]) - y0
    xd = numpy.linspace(0,arr.shape[1]-1,arr.shape[1]) - x0
    yyd,xxd = numpy.meshgrid(yd,xd)
    dmap = numpy.sqrt(yyd**2 + xxd**2)

    return dmap
    
def shiftSpot1d(spot, dx, dy, kernels=None, kargs=None):
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
    
def poo(arr, dx, dy, splines=None, binFactor=10, padTo=0, applyPixelResp=False, kargs=None):
    assert dx>=0 and dy>=0

    # Trim raw image to multiples of binFactor pixels.
    maxSize = binFactor*(numpy.array(arr.shape,dtype='i4')/binFactor)
    newSize = (maxSize / binFactor).tolist()
    arr = arr[:maxSize[0],:maxSize[1]].copy()

    # Get our unshifted, binned, reference image.
    if applyPixelResp:
        arrs = applyPixelResponse(arr, binFactor)
        arr00 = pfs_tools.rebin(arrs, *newSize)
    else:
        arr00 = pfs_tools.rebin(arr, *newSize)

    if padTo:
        arr00unpadded = arr00.copy()
        arr00, padSlice = padArray(arr00, padTo)
        
    # Put the rest of this in loop....
    
    # interpolation-shift the binned reference image
    if splines == None:
        splines = 'shift1d'
    splineFuncs = dict(shift2d=shiftSpot2d,
                       shift1d=shiftSpot1d,
                       shiftSpline=shiftSpotSpline,
                       spline=shiftSpotSpline)
    splineFunc = splineFuncs[splines]
    arrShifted, kernels = splineFunc(arr00, float(dx)/binFactor, float(dy)/binFactor, kargs=kargs)

    # And pixel-shift (and optionally apply pixel response) then pixel-bin a comparison image
    arrPlaced = arr * 0
    arrPlaced[dy:,dx:] = arr[slice(None,-dy if dy else None),
                             slice(None,-dx if dx else None)]
    if applyPixelResp:
        arrPlaced = applyPixelResponse(arrPlaced, binFactor)
    arrPlaced = pfs_tools.rebin(arrPlaced, *newSize)

    if padTo:
        arrShifted = unpadArray(arrShifted, padSlice)
        arr00up = unpadArray(arr00, padSlice)
        assert numpy.all(numpy.equal(arr00unpadded, arr00up))
    
    return arr00, arrPlaced.copy(), arrShifted.copy(), kernels

def shiftSpotBy(spot, shiftBy, binTo,
                applyPixelResp=False, shiftFunc=shiftSpot1d,
                oversampleFactor=10, doNorm=False):
    """ shift oversampled spot by 'shiftBy' unbinned pixels, both by shifting and placing and 
    by binning and shifting. 
    """

    binnedOversample = oversampleFactor/binTo
    assert binnedOversample == float(oversampleFactor)/binTo

    try:
        dx, dy = shiftBy
    except:
        dx, dy = shiftBy, shiftBy
        
    binnedShape = (numpy.array(spot.shape,'i2')/(binnedOversample)).tolist()
    binnedSpot = pfs_tools.rebin(spot, *binnedShape)
    if doNorm:
        binnedSpot = binnedSpot / binnedSpot.sum()
    binnedShift = (float(dx)/binnedOversample, 
                   float(dy)/binnedOversample)
    if "debug":
        print ("binTo=%d; oversampling = %d; binnedOversample = %s; binnedShape = %s; binnedShift = %s; dx, dy = %d, %d" % 
               (binTo,
                oversampleFactor,
                binnedOversample,
                binnedShape,
                binnedShift, 
                dx, dy))

    if numpy.abs(binnedShift[0]) >= 1 or numpy.abs(binnedShift[1]) >= 1:
        idx = int(binnedShift[0])
        idy = int(binnedShift[1])
        print "cannot interpolate more than 1 pixel, moving by %d,%d first" % (idx,idy)
        binnedShift = (binnedShift[0]-idx, binnedShift[1]-idy)
        tSpot = numpy.zeros(binnedSpot.shape, dtype='f4')
        tSpot[idy:,idx:] = binnedSpot[slice(None,-idy if idy else None),
                                      slice(None,-idx if idx else None)]
        binnedSpot = tSpot

    # interpolation-shift the binned image
    print "shifting %s by %g,%g" % (binnedSpot.shape, binnedShift[0], binnedShift[1])
    shiftedSpot, kernels = shiftFunc(binnedSpot, *binnedShift)

    # And pixel-shift (and optionally apply pixel response) then pixel-bin a comparison image
    print "placing %s by %d,%d" % (spot.shape, dx,dy)
    placedSpot = numpy.zeros(spot.shape, dtype='f4')
    placedSpot[dy:,dx:] = spot[slice(None,-dy if dy else None),
                               slice(None,-dx if dx else None)]
    if applyPixelResp:
        placedSpot = applyPixelResponse(placedSpot, oversampleFactor)

    placedSpot = pfs_tools.rebin(placedSpot, *binnedShape)
    if doNorm:
        placedSpot = placedSpot / placedSpot.sum()

    return placedSpot, shiftedSpot, kernels

def dispSpot(spotDict, key):
    bin, pad, shift = key
    res = spotDict[key]
    images = res[0]
    freqs = res[1]

    f = plt.figure()

def parts(spotDict, key):
    scale = key[0]**2

    spotShift = spotDict[key][0][2] / scale
    spotPlace = spotDict[key][0][1] / scale
    dspot = spotShift - spotPlace
    dspotNorm = dspot / spotPlace
    freqSpot = spotDict[key][1][2]
    return (spotShift,
            spotPlace,
            dspot,
            dspotNorm,
            freqSpot)
    
def saveForBin(spotDict, bin, filebase):
    filename = "%s_%d.fits" % (filebase, bin)

    hdus = []
    hdus.append(pyfits.PrimaryHDU(spotDict[(bin,0,0,0)][1][2].astype('f4')))
    for s in range(bin):
        skip0 = False
        for dir in (0,s):
            if skip0:
                continue
            if s == 0 and dir == 0:
                skip0 = True
                
            sparts = parts(spotDict, (bin, 0, s, dir))
            hdu = pyfits.ImageHDU(sparts[1].astype('f4'))
            hdu.header['XSHIFT'] = s
            hdu.header['YSHIFT'] = dir
            hdu.header['TYPE'] = 'placed'
            hdus.append(hdu)
            print "added placed hdu for (%d,%d,%d)" % (bin, s, dir)
            hdu = pyfits.ImageHDU(sparts[0].astype('f4'))
            hdu.header['XSHIFT'] = s
            hdu.header['YSHIFT'] = dir
            hdu.header['TYPE'] = 'shifted'
            hdus.append(hdu)
            print "added shifted hdu for (%d,%d,%d)" % (bin, s, dir)
    print "writing %d HDUs" % (len(hdus))
    
    hdulist = pyfits.HDUList(hdus)
    hdulist.writeto(filename, clobber=True)
                
    
def gatherPoo(spot, splines=None, applyPixelResp=False, kargs=None):
    tries = ([1,0],
             [4,0],
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
                      splines=splines,
                      applyPixelResp=applyPixelResp,
                      binFactor=bin, 
                      padTo=pad,
                      kargs=kargs)
            freq = imfreq1d(ret[2])
            all[(bin, pad, shift)] = ret, freq

    return all, kargs
        
def gatherPoo2(spot, bin, splines=None, pad=0, applyPixelResp=False, kargs=None):

    if not kargs:
        kargs = dict(n=3)
        
    all = {}
    for xshift in range(bin):
        for yshift in range(bin):
            print "processing %s with kargs=%s" % ((bin,pad,xshift,yshift), kargs)
            ret = poo(spot, xshift, yshift, 
                      splines=splines,
                      applyPixelResp=applyPixelResp,
                      binFactor=bin, 
                      padTo=pad,
                      kargs=kargs)
            freq = imfreq1d(ret[2])
            all[(bin, pad, xshift,yshift)] = ret, freq

    return all, kargs
        
            
def spotShowImages(im, fftim, iim, imbox, fig, plotids, resClip=0, 
                   pixLims=None,pixscale=0.015,
                   showAliasAt=None, colorbars=(),):

    if im != None:
        plt_im = fig.add_subplot(plotids[0])
        plt_im.xaxis.set_visible(False)
        #plt_im.imshow(plotutils.asinh(im, non_linear=0.0001, scale_max=0.2),extent=imbox)
        plt_im.imshow(im, extent=imbox, interpolation='nearest')
        plt_im.set_ylabel('%g um pixels' % (1000*pixscale))

        l0,l1 = plt_im.get_xlim()
        plt_im.vlines(0,l0,l1,'r', alpha=0.3)
        plt_im.hlines(0,l0,l1,'r', alpha=0.3)
    else:
        plt_im = None

    if fftim != None:
        plt_fft = fig.add_subplot(plotids[1])
        #plt_fft.yaxis.set_visible(False)
        plt_fft.xaxis.set_visible(False)
        plt_fft.imshow(numpy.log10(numpy.abs(fftim)), extent=imbox)

        if showAliasAt:
            plt_fft.vlines(-showAliasAt, -showAliasAt, showAliasAt+1, 'r', alpha=0.6)
            plt_fft.vlines(showAliasAt+1, -showAliasAt, showAliasAt+1, 'r', alpha=0.6)
            plt_fft.hlines(-showAliasAt, -showAliasAt, showAliasAt+1, 'r', alpha=0.6)
            plt_fft.hlines(showAliasAt+1, -showAliasAt, showAliasAt+1, 'r', alpha=0.6)
    else:
        plt_fft = None

    if iim != None:
        plt_iim = fig.add_subplot(plotids[2])
        plt_iim.xaxis.set_visible(False)
        #if resClip:
        #    residIm = numpy.clip(iim,-resClip,resClip)
        #else:
        #    residIm = iim
        residIm = iim

        im3 = plt_iim.imshow(residIm, extent=imbox, vmin=-resClip, vmax=resClip)
        ticks = numpy.sort(numpy.append(numpy.linspace(-resClip,resClip,4),[0]))
        print "setting ticks to %s" % (ticks)
        fig.colorbar(im3, ticks=ticks)
        l0,l1 = plt_iim.get_xlim()
        plt_iim.vlines(0,l0,l1,'r', alpha=0.3)
        plt_iim.hlines(0,l0,l1,'r', alpha=0.3)

    else:
        plt_iim = None

    if pixLims:
        for i,p in enumerate([plt_im,plt_fft,plt_iim]):
            if pixLims[i] and p != None:
                p.set_xlim(-pixLims[i],pixLims[i],auto=False)
                p.set_ylim(-pixLims[i],pixLims[i],auto=False)


    return plt_im, plt_fft, plt_iim

def spotShow(im, scale=10, binning=2):
    fig = plt.figure('spot')
    sgs = gridspec.GridSpec(3,3, hspace=0.1, wspace=0.1)

    im = im/im.sum()
    fftx, ffty, fftim = imfreq1d(im, sampling=1.0/scale)
    iim = freq2im(fftim, doAbs=False)
    imbox = imextent(im, scale=scale, doCenter=True)
    spotShowImages(im, fftim, iim, imbox, fig, [sgs[0,0],sgs[0,1],sgs[0,2]])

    p1 = fig.add_subplot(sgs[6:])
    p1.plot(fftx/scale, numpy.abs(ffty))
    p1.set_yscale('log')

    im = pfs_tools.rebin(im, 
                         im.shape[0]/binning,
                         im.shape[0]/binning)
    binnedScale = float(scale)/binning
    im = im/im.sum()
    fftx, ffty, fftim = imfreq1d(im, sampling=1.0/binnedScale)
    iim = freq2im(fftim, doAbs=False)
    imbox = imextent(im, scale=binnedScale, doCenter=True)
    spotShowImages(im, fftim, iim, imbox, fig, [sgs[1,0],sgs[1,1],sgs[1,2]])
    p1.plot(fftx/binnedScale, numpy.abs(ffty))

    fig.show()

def frdShow(frds,
            band='IR', date='2013-04-18',
            focus=0,
            fiberIdx=0, wavelength=11067, titleExtra='',
            figName='frd', doClear=True,
            doNorm=True, yrange=None):

    
    fig = plt.figure(figName)
    if doClear:
        fig.clf()

    spots = []
    for f in frds:
        d = jegSpots.readSpotFile(pathSpec=dict(band=band,
                                                date=date,
                                                focus=focus,
                                                frd=f))
        spot_w = numpy.where((d['wavelength'] == wavelength) &
                             (d['fiberIdx'] == fiberIdx))[0]
        if len(spot_w) != 1:
            raise RuntimeError("you did not specify a unique spot (%d)" % (len(spot_w)))

        spots.append(d['spot'][spot_w][0])
        print "spot %d max=%g sum=%g" % (f, spots[-1].max(), spots[-1].sum())

    p1 = fig.add_subplot(1,1,1)
    
    refSpot = spots[0]
    normScale = 1.0*refSpot.max()
    midline = refSpot.shape[0]/2
    x = (numpy.arange(2*midline)-midline)/10.0
    colors = ['b','g','r']
    p1.xaxis.grid(True, which='major',
                  markevery=1,
                  color="#a0a0a0", linestyle='-')
    p1.set_title('%s, fiber %d, wavelength %d %s' % (band, 
                                                     fiberIdx,
                                                     wavelength,
                                                     titleExtra))
    p1.set_xlabel('15um pixels')
    p1.set_ylabel('of peak flux in %0.3f sigma spot' % (frds[0]/1000.0))
    p1.set_autoscalex_on(False)
    p1.set_xlim(-5, 5)
    p1.xaxis.set_ticks(numpy.arange(-5,6))

    for i in range(1,len(frds)):
        s = spots[i]
        dline = s[midline,:]/normScale - refSpot[midline,:]/normScale
        p1.plot(x, dline, '+-%s' % (colors[i-1]), label='frd%d - frd%s' % (frds[i], frds[0]))
        dline = s[:,midline]/normScale - refSpot[:,midline]/normScale
        p1.plot(x, dline, '+-%s' % (colors[i-1]))

    p1.set_autoscaley_on(False)
    if yrange != None:
        p1.set_ylim(*yrange)

    p1.plot(x, refSpot[midline,:]/normScale, '-', color='gray', alpha=0.5,
            label='%0.3f sigma FRD' % (frds[0]/1000.0))
    p1.plot(x, refSpot[:,midline]/normScale, '-', color='gray', alpha=0.5)

    p1.legend()

    return spots
                               

def spotShow2(spot0, scale1, scale2, figname='spot', 
              plotWidth=8,
              maxFreq=5, shiftBy=None, 
              applyPixelResp=False,
              unbinnedScale=10, 
              shiftFunc=shiftSpotSpline):
    fig = plt.figure(figname)
    sgs = gridspec.GridSpec(3,3, hspace=0.05, wspace=0.05)

    assert spot0.shape[0] == spot0.shape[1]
    fullSize = spot0.shape[0]
    spot0 = spot0/spot0.sum()

    size1 = scale1 * fullSize/unbinnedScale
    size2 = scale2 * fullSize/unbinnedScale
    print "sizes = %s, %s" % (size1, size2)

    spot1 = pfs_tools.rebin(spot0, size1, size1)
    spot2 = pfs_tools.rebin(spot0, size2, size2)

    spot1 = spot1/spot1.sum()
    spot2 = spot2/spot2.sum()

    if shiftBy != None:
        placedSpot1, shiftedSpot1, kernels = shiftSpotBy(spot0, shiftBy, scale1,
                                                         shiftFunc=shiftFunc, 
                                                         applyPixelResp=applyPixelResp,
                                                         doNorm=True)
    else:
        placedSpot1 = numpy.zeros(spot1.shape)
        shiftedSpot1 = numpy.zeros(spot1.shape)

    fftx, ffty, fftim1 = imfreq1d(spot1, sampling=1.0/scale1)
    imbox1 = imextent(spot1, doCenter=True)
    imlims = [-plotWidth*scale1, plotWidth*scale1]*2

    p1 = fig.add_subplot(sgs[6:])
    if maxFreq:
        p1.set_autoscalex_on(False)
        p1.set_xlim([0, maxFreq])
    p1.set_yscale('log')
    p1.set_xlabel('cycles per 15um CCD pixel')
    p1.plot(fftx, numpy.abs(ffty), label='1.5um pixel')

    fftx, ffty, fftim2 = imfreq1d(spot2, sampling=1.0/scale2)
    imbox2 = imextent(spot2, scale=scale2, doCenter=True)
    imlims = [-plotWidth*scale2, plotWidth*scale2] # imextent(spot1, scale=scale1, doCenter=True)

    if shiftBy != None:
        placedSpot2, shiftedSpot2, kernels = shiftSpotBy(spot0, shiftBy, scale2,
                                                         shiftFunc=shiftFunc,
                                                         applyPixelResp=applyPixelResp,
                                                         doNorm=True)
    else:
        placedSpot2 = numpy.zeros(spot2.shape)
        shiftedSpot2 = numpy.zeros(spot2.shape)
    diffIm2 = placedSpot2 - shiftedSpot2
    print "flux=%g,%g,resid=%g" % (placedSpot2.sum(),shiftedSpot2.sum(),
                                   numpy.abs(diffIm2).sum())

    trimRad = spot2.shape[0]/2
    xr = numpy.arange(spot1.shape[0]) - (spot1.shape[0]+1)/2
    mx,my = numpy.meshgrid(xr,xr)
    trimMask = ((mx < -trimRad) | (mx > trimRad) |
                (my < -trimRad) | (my > trimRad))
    trimMask = distmap(fftim1) > trimRad
    maskedFftim = fftim1 * trimMask
    maskedSpot = freq2im(maskedFftim)
    diffIm1 = maskedSpot

    basePixScale=0.015
    pixLim = 6
    if 'lastminute':
        spotShowImages(spot1, fftim1, None, imbox1, fig, [sgs[0,0],sgs[0,1],None], resClip=0.002,
                       pixLims=[pixLim*scale1]*3, pixscale=basePixScale/scale1, showAliasAt=size2/2.0)
        plots = spotShowImages(spot2, None, diffIm2, imbox2, fig, [sgs[1,0],None,sgs[1,2]], resClip=0.003,
                               pixLims=[pixLim*scale2,None,pixLim*scale2], pixscale=basePixScale/scale2,
                               colorbars=[1])
        plots[2].set_title('placed - shifted\nfraction of total spot flux')
    else:
        spotShowImages(spot1, fftim1, diffIm1, imbox1, fig, [sgs[0,0],sgs[0,1],sgs[0,2]], resClip=0.002,
                       pixLims=[pixLim*scale1]*3, pixscale=basePixScale/scale1, showAliasAt=size2/2.0)
        spotShowImages(spot2, fftim2, diffIm2, imbox2, fig, [sgs[1,0],sgs[1,1],sgs[1,2]], resClip=0.003,
                       pixLims=[pixLim*scale2,None,pixLim*scale2], pixscale=basePixScale/scale2)

    p1.plot(fftx, numpy.abs(ffty), label="15um pixel")
    p1.legend()

    # fig.suptitle('Interpolation at IR detector center')

    return fig

def plot2bins(spot, bin1, bin2, trimmedSize=None, 
              shiftBy=None, plotWidth=6,
              figname='spot', unbinnedScale=10, 
              shiftFunc=shiftSpotSpline, maxFreq=2):

    fullSize = 256
    assert spot.shape == (fullSize, fullSize)

    if trimmedSize:
        inset = (fullSize-trimmedSize)/2
        spot = spot[inset:-inset,inset:-inset].astype('f8')
        
    return spotShow2(spot, bin1, bin2, 
                     plotWidth=plotWidth,
                     shiftBy=shiftBy, figname=figname, 
                     unbinnedScale=unbinnedScale,
                     shiftFunc=shiftFunc,maxFreq=maxFreq)
