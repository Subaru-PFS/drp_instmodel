import numpy
import os
import scipy

from .utils import pydebug, blackbody

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
        
        if wave is None:
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
    
class SlopeSpectrum(Spectrum):
    def __init__(self, detector, gain=1.0):
        self.detector = detector
        self.scale = 1e2 * gain

    def flux(self, wave):
        """ Return a flat spectrum, with a slight blue-up-to-red tilt. """

        # Make 3800..12700 go to 1-(0.062) to 1+(0.27)
        flux = wave/10000.0 - 1
        flux = flux/10 + 1

        return wave, flux * self.scale

class FlatSpectrum(Spectrum):
    def __init__(self, detector, gain=1.0):
        self.detector = detector
        self.scale = 1e11 * gain

    def flux(self, wave):
        """ return a quartz lamp spectrum, as seen by our detector. """

        # Work out the fing scaling, CPL
        return wave, blackbody(wave*10.0, 3800.0) * self.scale

class LineSpectrum(Spectrum):
    """ """

    def linelist(self, minWave, maxWave):
        """ Return all the lines in the given wave range. """

        raise NotImplementedError()
    
    def flux(self, wave):
        return self.linelist(wave.min(), wave.max())            

class ArcSpectrum(LineSpectrum):
    def __init__(self, lampset=None, gain=2.0):
        """ Create a spectrum which will return a flux of gain at every ~spacing AA, 0 elsewhere. 

        Actually, return a comb spectrum with non-zero values at the full range endpoints and at as
        many points as can be placed between them no closer than the given spacing.
        """

        self.lampset = lampset
        self.loadLines(lampset)
        self.gain = gain

    def __str__(self):
        return("ArcSpectrum(lampset=%s, gain=%s, nlines=%d, waverange=(%g,%g))" %
               (self.lampset, self.gain, len(self.lines),
                self.lines['wave'].min(), self.lines['wave'].max()))
                                                                
    def linelist(self, minWave=0, maxWave=1e6):
        """ Return all the lines in the given wave range. """

        w_w = (self.lines['wave'] >= minWave) & (self.lines['wave'] <= maxWave)
        
        wave = self.lines['wave'][w_w]
        flux = self.lines['flux'][w_w] * self.gain

        return wave, flux

    def loadLines(self, lampset=None):
        """ Hackety hack hack hack. """
    
        # Map self.band to a filename using some config file. For now, hardcode
        dataRoot = os.environ.get('DRP_INSTDATA_DIR', '.')
        filepath = os.path.join(dataRoot, 'data', 'lines', 'nist_all.txt')

        self.lines = numpy.genfromtxt(filepath, usecols=list(range(3)),
                                      dtype=[('wave', 'f4'),
                                             ('name', 'U5',),
                                             ('flux', 'f4')])
        if lampset is not None:
            select = numpy.zeros(len(self.lines), dtype=bool)
            for lamp in lampset:
                select |= self.lines['name'] == lamp
            self.lines = self.lines[select]

        return self.lines
        

class CombSpectrum(LineSpectrum):
    def __init__(self, spacing=50, gain=10000.0, inset=10):
        """ Create a spectrum which will return a flux of gain at every ~spacing AA, 0 elsewhere. 

        Actually, return a comb spectrum with non-zero values at the full range endpoints and at as
        many points as can be placed between them no closer than the given spacing.
        """
        self.spacing = spacing
        self.gain = gain
        self.inset = inset
        
    def linelist(self, minWave, maxWave):
        """ Return all the lines in the given wave range. """
        
        dw = maxWave - minWave - 2*self.inset

        wave = numpy.linspace(minWave + self.inset, minWave + dw, dw/self.spacing + 1)
        flux = wave * self.gain

        return wave, flux


class TextSpectrum(Spectrum):
    """A spectrum read from a text file via ``numpy.genfromtxt``

    If no additional parameters are supplied, we'll guess that the file is
    comprised of wavelength and flux columns, in that order.

    Parameters
    ----------
    filename : `str`
        Path to file to read.
    wavelengthScale : `float`, optional
        Scale by which to multiply wavelengths.
    fluxScale : `float`, optional
        Scale by which to multiply fluxes.
    **kwargs : `dict`, optional
        Keyword arguments for ``numpy.genfromtxt``.
    """
    def __init__(self, filename, wavelengthScale=1.0, fluxScale=1.0, **kwargs):
        if not kwargs:
            kwargs["dtype"] = [('wavelength', 'f4'), ('flux', 'f4')]
        data = numpy.genfromtxt(filename, **kwargs)
        super().__init__(data['wavelength']*wavelengthScale, data['flux']*fluxScale)


class NullSpectrum(Spectrum):
    """A spectrum that is zero everywhere"""
    def __init__(self):
        pass

    @property
    def interp(self):
        def dummyInterpolator(*args, **kwargs):
            return 0.0
        return dummyInterpolator
