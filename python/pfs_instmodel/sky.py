__all__ = ["SkyModel", "StaticSkyModel"]

import os
import numpy
import scipy.interpolate

import pydebug

class SkyModel(object):
    """ Encapsulate a generator for sky spectra. We specify the model and the individual multiplicative and additive terms.

    For starters, we have no terms: there is only one sky.

    model = SkyModel('IR')
    """

    def __init__(self, band, centerAlt=None, centerAz=None, t0=None):
        """ Create a generative model of the sky.


        """
        self.band = band
        self.centerAlt = centerAlt
        self.centerAz = centerAz
        self.t0 = t0

        self._makeSkySpline()
        self._makeExtinctionSpline()

    def __str__(self):
        return ("%s(%s[%0.2f to %0.2f], alt=%s, az=%s t0=%s)" %
                (self.__class__.__name__, self.band,
                 self.minWave, self.maxWave,
                 self.centerAz, self.centerAlt, self.t0))

    def _makeSkySpline(self):
        raise NotImplementedError("_makeSkySpline")
    def _makeExtinctionSpline(self):
        raise NotImplementedError("_makeExtinctionSpline")

    def getSkyAt(self, **argv):
        """ Return a spline for the sky at the given conditions. """
        raise NotImplementedError("getSkyAt")

    def evalSkyAt(self, spacing=1.0, **argv):
        """ evaluate the sky at the given conditions.

        Args:
         - spacing    - the spacing of the wavelengths to evalute at. [1.0 AA]
         - **argv     - passed down to .getSkyAt()

         Returns:
          wave, sky   - arrays 
        """
        x = numpy.arange(self.minWave, self.maxWave+spacing, spacing)
        sky = self.skySpline(x)
        return x, sky

    def evalExtinctionAt(self, spacing=1.0, **argv):
        x = numpy.arange(self.minWave, self.maxWave+spacing, spacing)
        sky = self.extinctionSpline(x)
        return x, sky

    @property
    def minWave(self):
        return self.skySpline.get_knots().min()

    @property
    def maxWave(self):
        return self.skySpline.get_knots().max()

    @property
    def waveRange(self):
        return self.minWave, self.maxWave


class StaticSkyModel(SkyModel):

    def _makeSkySpline(self):
        """ Read in JEG's preliminary sky model. """

        # Map self.band to a filename using some config file. For now, hardcode
        dataRoot = os.environ.get('PFS_INSTDATA_DIR', '.')
        filepath = os.path.join(dataRoot, 'data', 'sky', 'sumire%sskyHR.dat' % (self.band.upper()))

        a = numpy.genfromtxt(filepath, skip_footer=20, comments='\\')
        self.skySpline = scipy.interpolate.InterpolatedUnivariateSpline(a[:,0], a[:,3], k=3)

    def _makeExtinctionSpline(self):
        """ Read in JEG's preliminary extinction model. """

        # Map self.band to a filename using some config file. For now, hardcode
        dataRoot = os.environ.get('PFS_INSTDATA_DIR', '.')
        filepath = os.path.join(dataRoot, 'data', 'sky', 'MKextinction.dat')

        a = numpy.genfromtxt(filepath, comments='\\')
        w = numpy.where((a[:,0] >= self.minWave) & (a[:,0] <= self.maxWave))
        self.extinctionSpline = scipy.interpolate.InterpolatedUnivariateSpline(a[w,0], a[w,1], k=3)

    def getSkyAt(self, **argv):
        """ Return a spline for the sky at the given conditions.

        Just for giggles, scale
        """
        return self.skySpline
