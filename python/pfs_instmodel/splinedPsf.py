
#!/usr/bin/env python

try:
    from IPython.core.debugger import Tracer; debug_here = Tracer()
except:
    pass

import os
import pyfits
import cPickle as pickle

import numpy as np
import scipy.interpolate as spinterp

import psf

class SplinedPsf(psf.Psf):
    def __init__(self, band, constructFromSpots=False):
        """ Read in our persisted form. Optionally that from the zemax spots. """

        self.band = band
        self.wave = []
        self.fiber = []
        self.spots = []
        self.imshape = (0,0)
        self.coeffs = []

        # Map self.band to filenames using some config file. For now, hardcode
        dataRoot = os.environ.get('PFS_INSTDATA_DIR', '.')
        psfFilepath = os.path.join(dataRoot, 'data', 'spots', band, 'psf.pck')

        if constructFromSpots:
            spotFilepath = os.path.join(dataRoot, 'data', 'spots', band, 'spots.fits')
            self.makeSplinesFromSpotFile(spotFilepath)
            self._saveToFile(psfFilepath)

        # We don't actually need to do this, but it seems like a good check.
        self.loadFromFile(psfFilepath)

    def __str__(self):
        nSpots = len(self.spots)
        if nSpots > 0:
            return ("%s<%d spots; %d wavelengths (%0.2fAA to %0.2fAA); %d fibers>" %
                    (self.__class__.__name__, nSpots,
                     len(np.unique(self.wave)), np.min(self.wave), np.max(self.wave),
                     len(np.unique(self.fiber))))
        else:
            return ("%s<0 spots>" %
                    (self.__class__.__name__))

    def makeSplinesFromSpotFile(self, filename):
        """ """
        sf = pyfits.open(filename, mode='readonly')
        s = sf[1].data

        self.wave = s['lambda']
        self.fiber = s['fiberIdx']
        self.spots = s['spot'].astype('float64')
        self.imshape = self.spots[0].shape

        self.coeffs = np.zeros(self.imshape,
                               dtype='O')

        # Make a spline for each pixel.
        # XXX - This fiber/wave indexing scheme is not safe in general
        xx = np.unique(self.fiber)
        yy = np.unique(self.wave)

        for ix in range(self.imshape[0]):
            print "splining col %d" % (ix)
            for iy in range(self.imshape[1]):
                self.coeffs[iy, ix] = spinterp.RectBivariateSpline(xx, yy,
                                                                  self.spots[:, iy, ix].reshape(len(xx), len(yy)))
    def psfsAt(self, fibers, waves):
        """ Return a stack of PSFs, instantiated on a rectangular grid.

        Args:
           fibers, waves - two arrays of positions; we instatiate ourselves at the rectangular grid of these.

        Returns:
           - a 3D array of images [len(xs)*len(ys), imshape.x, imshape.y]
           - a list of the (x,y) positions the images are instantiated at.
        """

        newImages = np.zeros(shape=(len(fibers)*len(waves),
                                       self.imshape[0],
                                       self.imshape[1]))

        centers = [(x,y) for x in fibers for y in waves]
        for ix in range(self.imshape[0]):
            print "col %d" % (ix)
            for iy in range(self.imshape[1]):
                newImages[:, iy, ix] = self.coeffs[iy, ix](fibers, waves).flat

        return newImages, centers

    def _saveToFile(self, filename):
        """ Saves our state (esp. the splines), currently to a pickle.

        Very very very unwarranted chumminess with the spline internals.
        Should implement __reduce__, etc. for pickling, FITS/HDF otherwise.

        """

        with open(filename, 'w+') as pfile:
            d = {}

            d['tck'] = [s.tck for s in self.coeffs.flatten()]
            d['fp'] = [s.fp for s in self.coeffs.flatten()]
            d['degrees'] = [s.degrees for s in self.coeffs.flatten()]

            d['wave'] = self.wave
            d['fiber'] = self.fiber
            d['spots'] = self.spots

            print "pickling...."
            pickle.dump(d, pfile, protocol=-1)

    def loadFromFile(self, filename):
        with open(filename, 'r') as pfile:
            d = pickle.load(pfile)

        self.spots = d['spots']
        self.wave = d['wave']
        self.fiber = d['fiber']
        self.imshape = self.spots[0].shape

        self.coeffs = np.zeros(self.imshape,
                               dtype='O')

        fp = np.array(d['fp']).reshape(self.imshape)
        degrees = np.array(d['degrees']).reshape(self.imshape + (-1,))
        tck = np.array(d['tck']).reshape(self.imshape + (-1,))

        for ix in range(self.imshape[0]):
            for iy in range(self.imshape[1]):
                c = spinterp.RectBivariateSpline.__new__(spinterp.RectBivariateSpline)
                c.fp = fp[iy,ix]
                c.degrees = degrees[iy,ix]
                c.tck = tck[iy,ix]
                self.coeffs[iy, ix] = c
