try:
    from IPython.core.debugger import Tracer; debug_here = Tracer()
except:
    pass

import pyfits
import cPickle as pickle

import numpy
import scipy.interpolate as sinterp

import psf

class SplinedPsf(psf.Psf):
    def __init__(self, filename=None):
        self.wave = []
        self.fiber = []
        self.spots = []
        self.imshape = (0,0)
        self.coeffs = []
        
        if filename:
            self.loadFromFile(filename)

    def makeSplinesFromSpotFile(self, filename):
        sf = pyfits.open(filename, mode='readonly')
        s = sf[1].data

        self.wave = s['lambda']
        self.fiber = s['fiberIdx']
        self.spots = s['spot'].astype('float64')
        self.imshape = self.spots[0].shape

        self.coeffs = numpy.zeros(self.imshape,
                                  dtype='O')

        # Make a spline for each pixel
        xx = numpy.unique(self.fiber)
        yy = numpy.unique(self.wave)
        
        for ix in range(self.imshape[0]):
            print "splining col %d" % (ix)
            for iy in range(self.imshape[1]):
                self.coeffs[iy, ix] = sinterp.RectBivariateSpline(xx, yy, 
                                                                  self.spots[:, iy, ix].reshape(len(xx), len(yy)))
    def psfsAt(self, fibers, waves):
        """ Return a stack of PSFs, instantiated on a rectangular grid.

        Args:
           fibers, waves - two arrays of positions; we instatiate ourselves at the rectangular grid of these.
           
        Returns:
           - a 3D array of images [len(xs)*len(ys), imshape.x, imshape.y]
           - a list of the (x,y) positions the images are instantiated at.
        """
        
        newImages = numpy.zeros(shape=(len(fibers)*len(waves),
                                       self.imshape[0],
                                       self.imshape[1]))

        centers = [(x,y) for x in fibers for y in waves]
        for ix in range(self.imshape[0]):
            print "col %d" % (ix)
            for iy in range(self.imshape[1]):
                newImages[:, iy, ix] = self.coeffs[iy, ix](fibers, waves).flat

        return newImages, centers
    
    def saveToFile(self, filename):
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

        self.coeffs = numpy.zeros(self.imshape,
                                  dtype='O')

        fp = numpy.array(d['fp']).reshape(self.imshape)
        degrees = numpy.array(d['degrees']).reshape(self.imshape + (-1,))
        tck = numpy.array(d['tck']).reshape(self.imshape + (-1,))

        for ix in range(self.imshape[0]):
            for iy in range(self.imshape[1]):
                c = sinterp.RectBivariateSpline.__new__(sinterp.RectBivariateSpline)
                c.fp = fp[iy,ix]
                c.degrees = degrees[iy,ix]
                c.tck = tck[iy,ix]
                self.coeffs[iy, ix] = c

                
