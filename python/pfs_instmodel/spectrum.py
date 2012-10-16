import numpy
import scipy

class Spectrum(object):

    def __init__(self, wave, flux):
        
        # Force a copy of the raw data
        self._wave = wave.astype('f4')
        self._flux = flux.astype('f4')

        self._interp = None

    @property
    def interp(self):
        """ Return (after possibly constructing) a flux = interp(wave) interpolator. 

            We will probably want other interpolators. For now, just spline.
        """
        
        if not self._interp:
            # We have to make sure that the spline is constructed with increasing wavelengths 
            w_i = numpy.argsort(self._wave)
            self._interp = scipy.interpolate.InterpolatedUnivariateSpline(self._wave[w_i], 
                                                                          self._flux[w_i], k=3)
        return self._interp

    def flux(self, wave=None):
        """ Return the flux at the given wavelengths, or all the raw values. 

        Returns:
          wave    - the raw wavelengths if the wave argument was None, else the wave argument
          flux    - the corresponding flux 
        """

        if wave == None:
            wave = self._wave

        return wave, self.interp(wave)
        
