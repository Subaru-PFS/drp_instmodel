import os
import numpy
import scipy.interpolate

from .utils import configFile
from . import exposure
from .arm import Arm


class Detector(object):
    """ Placeholder for _all_ per-camera properties. 'Detector' is abused, since there
        are camera, instrument, and telescope properties folded into this. So this will be 
        split and expanded when requirements become clearer.

         * bad columns
         * vignetting
         * gain, readnoise (have)
         * camera throughput (kinda have).
    """

    def __init__(self, detectorName, dtype='u2'):
        """ Need to work out how to specify & pass in our properties. Use a "struct" for now. """

        if (len(detectorName) != 2 or
            detectorName[0] not in {'b', 'r', 'n', 'm'} or
            detectorName[1] not in {'1','2','3','4'}):
            raise ValueError('detectorName must be in the form "b2" (%s)' % (detectorName))

        self.detectorName = detectorName

        # Map self.band to a filename using some mapper. For now, hardcode
        dataRoot = os.environ.get('DRP_INSTDATA_DIR', '.')
        filepath = os.path.join(dataRoot, 'data', 'detectors', '%s.py' % (self.armName))
        self.config = configFile.readfile(filepath)
        self.readThroughputSpline()
        self.dtype = dtype

    @property
    def arm(self):
        return Arm(self.detectorName[0])
    @property
    def armName(self):
        armNames = {Arm.BLUE: 'Blue', Arm.RED: 'Red', Arm.NIR: 'NIR', Arm.MID: 'Mid'}
        return armNames[self.arm]
    @property
    def spectrograph(self):
        return int(self.detectorName[1])

    @property
    def gain(self):
        """Electronic gain, in e/ADU"""
        return self.config['gain']

    @property
    def xcPixOffset(self):
        return self.config['ccdSize'][1]/2
    @property
    def ycPixOffset(self):
        return self.config['ccdSize'][0]/2
    @property
    def xcMmOffset(self):
        return self.xcPixOffset*self.config['pixelScale']
    @property
    def ycMmOffset(self):
        return self.ycPixOffset*self.config['pixelScale']

    def makeExposure(self, dtype=None):

        if dtype is None:
            dtype = self.dtype
        exp = exposure.Exposure(self, dtype=dtype)

        return exp

    def getBias(self, exp=None):
        """ Return a bias plane. """

        dtype = 'u2'
        bias = numpy.fix(numpy.random.normal(self.config['bias'],
                                             self.config['readNoise'],
                                             self.config['ccdSize'])).astype(dtype)
        return bias

    def addBias(self, exp, ontoBias=None, useDark=False):
        """ Add our bias to the given exposure. """

        if ontoBias is not None and self.arm != Arm.NIR:
            bias = exp.loadBias(ontoBias, useDark=useDark)
        else:
            bias = self.getBias(exp)

        exp.addPlane('bias', bias)

        return bias

    def saturate(self, exp):
        """ Deal with saturated pixels. Set flags, spread flux, etc. """
        pass
    
    def readout(self, exp, flux, exptime=1.0, ontoBias=None,
                applyFlat=None):
        """ 'Readout' an exposure: add bad columns, hot pixels, etc.

        Input 'flux' image is in units of electrons.
        """

        if applyFlat:
            flat = exp.loadFlat()
            flux0 = flux.copy()
            flux = flux*flat
            exp.addPlane('flat', flux-flux0)
            del flux0
            
        bias = self.addBias(exp, ontoBias=ontoBias, useDark=(exptime > 30))

        rimage = numpy.round(flux/self.gain) + bias  # in ADU
        saturatedPixels = (rimage > 65535)
        lowPixels = (rimage < 0)
        
        rimage = numpy.round(rimage).astype('u2')
        rimage[saturatedPixels] = 65535
        rimage[lowPixels] = 0
        exp.pixelImage = rimage

    def readThroughputSplineJEG(self):
        """ Read in JEG's preliminary detector response.

        This is saved in a text file, whose meat is:

        \  pix   lam(A)     sky(nM) sky(AB)  TPinst  TPsky  cntAB22.5   cntSKY    S/N
        \
        0   9700.40     87.54   17.64   0.060   0.055   1.07e+01   9.42e+02  0.34
        1   9701.21    157.72   17.01   0.060   0.056   1.08e+01   1.71e+03  0.26

        We are interested in the instrument throughput, TPinst
        """

        # Map self.arm to a filename using some mapper. For now, hardcode
        dataRoot = os.environ.get('DRP_INSTDATA_DIR', '.')
        filepath = os.path.join(dataRoot, 'data', 'sky', 'sumire%s.dat' % (self.armName))

        a = numpy.genfromtxt(filepath, comments='\\')
        a[:,1] *= 0.1
        self.throughputSpline = scipy.interpolate.InterpolatedUnivariateSpline(a[:,1], a[:,4],
                                                                               ext='zeros',k=3)
        return self.throughputSpline

    def readThroughputSpline(self):
        """ Read in Yabe's throughput data, as taken from Hirata's ETC.

        630.0   0.0000  1.0000  1.0000  1.0000  1.0000  0.0000
        640.0   0.1008  1.0000  1.0000  1.0000  1.0000  0.1008
        650.0   0.2357  1.0000  1.0000  1.0000  1.0000  0.2357

        We use the first and the second-last columns.
        """

        # Map self.arm to a filename using some mapper. For now, hardcode
        dataRoot = os.environ.get('DRP_INSTDATA_DIR', '.')
        name = {Arm.BLUE: "blu",
                Arm.RED: "red",
                Arm.NIR: "nir",
                Arm.MID: "mid"}[self.arm]
        filepath = os.path.join(dataRoot, 'data', 'throughput', 'pfs_thr_20211220_ext_all_%s.dat' % (name,))

        data = numpy.genfromtxt(filepath, comments='#')
        num = len(data)
        # pchip extrapolation only works if two end points are the same.
        # So pad things out
        wavelength = numpy.empty(num + 4, dtype=data.dtype)
        throughput = numpy.zeros(num + 4, dtype=data.dtype)
        wavelength[2:-2] = data[:, 0]
        throughput[2:-2] = data[:, 9]  # total throughput without the atmosphere
        dw = wavelength[3] - wavelength[2]
        wavelength[1] = wavelength[2] - dw
        wavelength[0] = wavelength[1] - dw
        dw = wavelength[-3] - wavelength[-4]
        wavelength[-2] = wavelength[-3] + dw
        wavelength[-1] = wavelength[-2] + dw
        assert throughput[0] == throughput[1]
        assert throughput[-1] == throughput[-2]

        self.throughputSpline = scipy.interpolate.PchipInterpolator(wavelength, throughput,
                                                                    extrapolate=True)
        return self.throughputSpline

    def throughput(self, waves):
        return self.throughputSpline(waves)
