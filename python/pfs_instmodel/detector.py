import os
import numpy
import scipy

import pfs_tools.configFile
import exposure

class Detector(object):
    """ Placeholder for _all_ per-camera properties. 'Detector' is abused, since there
        are camera, instrument, and telescope properties folded into this. So this will be 
        split and expanded when requirements become clearer.

         * bad columns
         * vignetting
         * gain, readnoise (have)
         * camera throughput (kinda have).
    """

    def __init__(self, detectorName):
        """ Need to work out how to specify & pass in our properties. Use a "struct" for now. """

        # Map self.band to a filename using some mapper. For now, hardcode
        dataRoot = os.environ.get('PFS_INSTDATA_DIR', '.')
        filepath = os.path.join(dataRoot, 'data', 'detectors', '%s.py' % (detectorName))
        
        self.band = detectorName
        self.config = pfs_tools.configFile.readfile(filepath)

    def makeEmptyExposure(self):
        return exposure.Exposure(self, dtype='u2')

    def simBias(self, shape=None):
        """ Return a bias image based on our properties. """
        
        exp = self.makeEmptyExposure()

        if not shape:
            shape = self.config['ccdSize']

        bias = numpy.random.normal(self.config['bias'],
                                   self.config['readNoise'],
                                   shape).astype('u2')
        ivar = bias*0 + 1/self.config['readNoise']**2
        exp.setImage(bias, ivar)

        return exp
 
    def getResponseSpline(self):
        """ Read in JEG's preliminary detector response.

        This is saved in a text file, whose meat is:

        \  pix   lam(A)     sky(nM) sky(AB)  TPinst  TPsky  cntAB22.5   cntSKY    S/N
        \
        0   9700.40     87.54   17.64   0.060   0.055   1.07e+01   9.42e+02  0.34
        1   9701.21    157.72   17.01   0.060   0.056   1.08e+01   1.71e+03  0.26

        We are interested in the instrument throughput, TPinst
        """

        # Map self.band to a filename using some mapper. For now, hardcode
        dataRoot = os.environ.get('PFS_INSTDATA_DIR', '.')
        filepath = os.path.join(dataRoot, 'data', 'sky', 'sumire%s.dat' % (self.band.upper()))

        a = numpy.genfromtxt(filepath, comments='\\')
        return scipy.interpolate.InterpolatedUnivariateSpline(a[:,1], a[:,4], k=3)
