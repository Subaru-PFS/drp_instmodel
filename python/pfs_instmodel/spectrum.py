import numpy
import scipy

import pfs_tools.blackbody

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

    def flux(self, wave=None, waverange=None):
        """ Return the flux at the given wavelengths, or all the raw values. 

        Returns:
          wave    - the raw wavelengths if the wave argument was None, else the wave argument
          flux    - the corresponding flux 
        """

        if wave == None:
            if waverange:
                w = numpy.where((self._wave >= waverange[0]) & (self._wave <= waverange[1]))
                return self._wave[w], self._flux[w]
            else:
                return self._wave, self._flux

        if waverange:
            w = numpy.where((wave >= waverange[0]) & (wave <= waverange[1]))
            wwave = wave[w]
            return wwave, self.interp(wwave)
        else:
            return wave, self.interp(wave)
        
    def __call__(self, wave):
        return self.flux(wave)
    
class FlatSpectrum(Spectrum):
    def __init__(self, detector):
        self.detector = detector

    def flux(self, wave):
        """ return a quartz lamp spectrum, as seen by our detector. 

        Notes
        -----

        We need to apply the instrument response, but that is not available for the new
        optics design. So we just return the full blackbody.

        We require wavelengths to evaluate at, and not just a range or a "full spectrum".
        """
        return wave, pfs_tools.blackbody(wave)
    

