import os
import numpy

import exposure

class Detector(object):
    def __init__(self, detectorName):
        """ Need to work out how to specify & pass in our properties. Use a "struct" for now. """

        # Map self.band to a filename using some config file. For now, hardcode
        dataRoot = os.environ.get('PFS_INSTDATA_DIR', '.')
        filepath = os.path.join(dataRoot, 'data', 'detectors', '%s.py' % (detectorName))
        
        self.config = self._hackConfig(filepath)

    def _hackConfig(self, configFile):
        gdict = {}
        ldict = {}
        
        execfile(configFile,  gdict, ldict)
        return ldict

    def makeEmptyExposure(self):
        return exposure.Exposure(self, dtype='u2')

    def simBias(self, shape=None):
        exp = self.makeEmptyExposure()

        if not shape:
            shape = self.config['ccdSize']

        bias = numpy.random.normal(self.config['bias'],
                                   self.config['readNoise'],
                                   shape).astype('u2')
        ivar = bias*0 + 1/self.config['readNoise']**2
        exp.setImage(bias, ivar)

        return exp
