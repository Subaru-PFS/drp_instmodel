import numpy
import scipy

import pfs_tools
from pfs_tools import pydebug

class Spectrum(object):
    """ Track the best information about a spectrum and provide sampling. 

    Note the following use cases:

     * Model spectra
     * Line lists
     * Sampled spectra

    For the sampled spectra (
    """
    
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

    def flux(self, wave=None, waveRange=None):
        """ Return the flux at the given wavelengths, or all the raw values. 

        Parameters
        ----------
        wave : array_like, optional
           The wavelengths to evaluate the flux at. If not passed in, use all the native 
           or sampled wavelengths.
        waveRange : (low, high), optional
           Clip the above wavelength vector to the samples between or including low & high.
           
        Returns
        -------
        wave :  the raw wavelengths if the wave argument was None, else the wave argument
        flux :  the corresponding flux 
        """

        assert waveRange is None or waveRange[0] <= waveRange[1]
        
        if wave == None:
            # No need to interpolate
            if waveRange:
                w = numpy.where((self._wave >= waveRange[0]) & (self._wave <= waveRange[1]))
                return self._wave[w], self._flux[w]
            else:
                return self._wave, self._flux

        if waveRange:
            w = numpy.where((wave >= waveRange[0]) & (wave <= waveRange[1]))
            wwave = wave[w]
            return wwave, self.interp(wwave)
        else:
            return wave, self.interp(wave)
        
    def __call__(self, wave):
        """ Shorthand to fetch the fluxes corresponding to one or more wavelengths. This is the 
        primary method to get fluxes. Subclasses only need to implement .flux(wave) """
        
        return self.flux(wave)[1]
    
class FlatSpectrum(Spectrum):
    def __init__(self, detector, gain=1.0):
        self.detector = detector
        self.scale = 1e11 * gain

    def flux(self, wave):
        """ return a quartz lamp spectrum, as seen by our detector. """

        # Work out the fing scaling, CPL
        return wave, pfs_tools.blackbody(wave, 3800.0) * self.scale

class LineSpectrum(Spectrum):
    """ """

    def flux(self, wave):
        return self.linelist(wave.min(), wave.max())
    
class CombSpectrum(LineSpectrum):
    def __init__(self, spacing=50, gain=10000.0, inset=100):
        """ Create a spectrum which will return a flux of gain at every ~spacing AA, 0 elsewhere. 

        Actually, return a comb spectrum with non-zero values at the full range endpoints and at as
        many points as can be placed between them no closer than the given spacing.
        """
        self.spacing = spacing
        self.gain = gain
        self.inset = inset
        
    def linelist(self, minWave, maxWave):
        """ Return all the lines in the given wave range. """
        
        dw = maxWave-minWave - 2*self.inset

        wave = numpy.linspace(minWave+self.inset, minWave+dw, dw/self.spacing + 1)
        flux = wave*0 + self.gain

        return wave, flux

    def flux0(self, wave):
        minWave = wave.min()
        maxWave = wave.max()
        nWaves = len(wave)

        dw = maxWave-minWave

        idx = numpy.linspace(0, nWaves-1, dw/self.spacing + 1).astype('i4')
        flux = wave*0
        flux[idx] = self.gain

        return wave, flux

