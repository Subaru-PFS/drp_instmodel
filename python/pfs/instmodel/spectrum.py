import numpy
import os
import scipy
import scipy.integrate
import astropy.io.fits
import astropy.wcs

from .utils import blackbody

SPEED_OF_LIGHT = 3.0e8*1.0e9  # nm/s


class Spectrum(object):
    """Base class for flux density as a function of wavelength

    All wavelengths are in nanometres (nm).

    All flux densities are F_nu in nanoJanskys (nJy).

    All (integrated) fluxes are in W/m^2.
    """
    def interpolate(self, wavelength):
        """Interpolate the spectrum at the nominated wavelength

        Parameters
        ----------
        wavelength : array_like
            Vector of wavelengths at which to interpolate, nm.

        Returns
        -------
        flux : array_like
            Vector of flux densities at the provided ``wavelength``s.
        """
        raise NotImplementedError("Method must be defined by subclass")

    def interpolateFrequency(self, frequency):
        """Interpolate the spectrum at the nominated freqency

        Parameters
        ----------
        frequency : array_like
            Vector of frequencies at which to interpolate, Hz.

        Returns
        -------
        flux : array_like
            Vector of flux densities at the provided ``frequency``s.
        """
        return self.interpolate(SPEED_OF_LIGHT/frequency)

    def _integrateImpl(self, lower, upper, **kwargs):
        """Integrate the spectrum between a single set of frequency bounds

        Note that this works with frequency rather than wavelength.
        It just does the integration, and doesn't worry about scaling the
        result to the correct units for external users. The spectrum is in nJy,
        and the integration is in frequency, so the output units should be
        nJy.Hz.

        Parameters
        ----------
        lower : `float`
            Lower frequency bound for the integration, Hz.
        upper : `float`
            Upper frequency bound for the integration, Hz.

        Returns
        -------
        flux : `float`
            Integrated flux between the frequency bounds, nJy.Hz.s
        """
        return scipy.integrate.quadrature(self.interpolateFrequency, lower, upper, **kwargs)[0]

    def integrate(self, lower, upper, **kwargs):
        """Integrate the spectrum between multiple wavelength bounds

        Parameters
        ----------
        lower : array_like
            Lower wavelength bounds for the integration, nm.
        upper : array_like
            Upper wavelength bounds for the integration, nm.

        Returns
        -------
        flux : array_like
            Integrated fluxes between the wavelength bounds, W/m^2.
        """
        lowerFreq = SPEED_OF_LIGHT/upper
        upperFreq = SPEED_OF_LIGHT/lower
        scale = 1.0e-9*1.0e-26  # nJy.Hz --> W/m^2
        return scale*numpy.vectorize(self._integrateImpl)(lowerFreq, upperFreq, **kwargs)

    def bounds(self):
        """Return the wavelength bounds of the spectrum.

        The bounds are defined such that outside them, the flux density is zero.

        Returns
        -------
        lower : `float`
            Lower bound of the spectrum, nm.
        upper : `float`
            Upper bound of the spectrum, nm.
        """
        return -numpy.inf, numpy.inf

    def __imul__(self, value):
        raise NotImplementedError("Method must be defined by subclass")

    def normalize(self, filterName, flux):
        """Normalize the spectrum to the nominated magnitude

        Parameters
        ----------
        filterName : `str`
            Name of filter bandpass. Must be listed in our menu of recognised
            filters.
        flux : `float`
            AB flux [nJy] in the filter bandpass.

        Returns
        -------
        norm : `float`
            Normalisation applied.
        """

        menu = {"i": "wHSC-i2.txt"}
        if filterName not in menu:
            raise RuntimeError(f"Unrecognised filter: {filterName}")

        bandpass = TextSpectrum(os.path.join(os.environ["DRP_INSTDATA_DIR"], "data", "hsc", menu[filterName]),
                                wavelengthScale=0.1)
        ab = ConstantSpectrum(flux)  # ABmag reference spectrum, nJy

        options = dict(tol=0.0, rtol=2.0e-3, maxiter=100)
        current = ProductSpectrum(self, bandpass, PhotonCounting()).integrate(*bandpass.bounds(), **options)
        expected = ProductSpectrum(ab, bandpass, PhotonCounting()).integrate(*bandpass.bounds(), **options)
        norm = expected/current
        self *= norm
        return norm

    def __mul__(self, other):
        """Multiply a spectrum by another

        E.g., to apply a transmission function.

        Parameters
        ----------
        other : `Spectrum`
            Spectrum to multiply by.

        Returns
        -------
        result : `Spectrum`
            Multiplied spectrum. May not be the same as ``self``!
        """
        raise NotImplementedError("Subclasses must implement")

    def hasFlux(self):
        """Return whether there is any flux

        This is intended for optimising the case where the fiber can be ignored.
        """
        return True


class TableSpectrum(Spectrum):
    """A Spectrum defined by a lookup table

    Parameters
    ----------
    wavelength : array_like
        Array of wavelengths, nm.
    flux : array_like
        Array of flux density (F_nu) at the corresponding wavelengths, nJy.
    """
    def __init__(self, wavelength, flux):
        # Force a copy of the raw data, and ensure it's sorted
        indices = numpy.argsort(wavelength)
        self.wavelength = wavelength[indices]
        self.flux = flux[indices]
        self.frequency = SPEED_OF_LIGHT/self.wavelength
        self._interp = lambda x: numpy.interp(x, self.wavelength, self.flux, 0.0, 0.0)

    def interpolate(self, wavelength):
        return self._interp(wavelength)

    def _integrateImpl(self, lower, upper):
        select = numpy.where(numpy.logical_and(self.frequency > lower, self.frequency < upper))[0]
        num = len(select)
        lowValue, highValue = self.interpolateFrequency(numpy.array([lower, upper]))
        if num == 0:
            return 0.5*(upper - lower)*(lowValue + highValue)
        frequency = self.frequency[select]
        flux = self.flux[select]
        if num > 1:
            result = -1*numpy.trapz(flux, frequency)  # -1 because frequency is monotonic decreasing
        else:
            result = 0
        # Add in the bit off the ends we haven't integrated
        lowFreq, highFreq = frequency[-1], frequency[0]
        lowFlux, highFlux = flux[-1], flux[0]
        result += 0.5*(lowFreq - lower)*(lowFlux + lowValue)
        result += 0.5*(upper - highFreq)*(highFlux + highValue)
        return result

    def bounds(self):
        return self.wavelength[0], self.wavelength[-1]

    def __imul__(self, value):
        self.flux *= value
        return self

    def resolution(self):
        """Return the average resolution, nm"""
        return (self.wavelength[1:] - self.wavelength[:-1]).mean()

    def __mul__(self, other):
        """Multiply spectrum by another

        We interpolate the lower-resolution spectrum to the higher-resolution
        and multiply.
        """
        if isinstance(other, TableSpectrum) and self.resolution() > other.resolution():
            wavelength = other.wavelength
            flux = self.interpolate(wavelength)*other.flux
        else:
            wavelength = self.wavelength
            flux = self.flux*other.interpolate(self.wavelength)
        return TableSpectrum(wavelength, flux)

    def __rmul__(self, other):
        """Multiply spectrum by another

        Reflected (swapped) arguments version.
        The operation commutes, so this is the same as ``__mul__``.
        """
        return self.__mul__(other)


class FunctionalSpectrum(Spectrum):
    """A spectrum generated from a functional form

    The functional form is provided by the ``interpolate`` method.
    """
    def interpolate(self, wavelength):
        raise NotImplementedError("Subclasses must define")

    def __mul__(self, other):
        """Multiply by another spectrum"""
        if not isinstance(other, FunctionalSpectrum):
            return NotImplemented
        return MultipliedSpectrum(self, other)


class MultipliedSpectrum(FunctionalSpectrum):
    """A spectrum generated from multiplying two `FunctionalSpectrum`s"""
    def __init__(self, left, right):
        self.left = left
        self.right = right

    def interpolate(self, wavelength):
        return self.left.iterpolate(wavelength)*self.right.interpolate(wavelength)


class SlopeSpectrum(FunctionalSpectrum):
    """A spectrum consisting of a slope in F_nu

    Parameters
    ----------
    scale : `float`
        Scale to apply to spectrum.
    """
    def __init__(self, scale=1.0):
        self.scale = scale

    def interpolate(self, wavelength):
        # Make 3800..12700 go to 1-(0.062) to 1+(0.27)
        return ((wavelength/10000.0 - 1)/10 + 1)*self.scale

    def __imul__(self, value):
        self.scale *= value
        return self

    def hasFlux(self):
        return self.scale > 0


class FlatSpectrum(FunctionalSpectrum):
    """A spectrum simulating a flat-field lamp: a black body

    Parameters
    ----------
    scale : `float`
        Scale to apply to spectrum.
    """
    def __init__(self, scale=7.0e-20):  # Scale chosen manually to get reasonable counts in 30 sec
        self.scale = scale

    def interpolate(self, wavelength):
        """return a quartz lamp spectrum, as seen by our detector."""
        return self.scale*blackbody(wavelength, 5000.0)

    def __imul__(self, value):
        self.scale *= value
        return self

    def hasFlux(self):
        return self.scale > 0


class LineSpectrum(Spectrum):
    """Base class for spectra consisting solely of emission lines

    Parameters
    ----------
    wavelength : array_like
        Array of wavelengths of the emission lines, nm.
    flux : array_like
        Array of flux of the corresponding lines, W/m^2.
    """
    def __init__(self, wavelength, flux):
        indices = numpy.argsort(wavelength)
        self.wavelength = wavelength[indices]
        self.frequency = SPEED_OF_LIGHT/self.wavelength
        self.flux = flux[indices]

    def _integrateImpl(self, lower, upper):
        """Integrate the spectrum between a single set of frequency bounds

        Note that, in contrast to `Spectrum`, this works with wavelength
        instead of frequency, because our fluxes are already integrated and
        there's no need to convert to frequency.

        Parameters
        ----------
        lower : `float`
            Lower wavelength bound for the integration, nm.
        upper : `float`
            Upper wavelength bound for the integration, nm.

        Returns
        -------
        flux : `float`
            Integrated flux between the wavelength bounds, W/m^2.
        """
        select = numpy.logical_and(self.wavelength > lower, self.wavelength < upper)
        return self.flux[select].sum()

    def integrate(self, lower, upper):
        return numpy.vectorize(self._integrateImpl)(lower, upper)

    def bounds(self):
        return self.wavelength[0], self.wavelength[-1]

    def __imul__(self, value):
        self.flux *= value
        return self

    def __mul__(self, other):
        return LineSpectrum(self.wavelength, self.flux*other.interpolate(self.wavelength))

    def hasFlux(self):
        return numpy.any(self.flux > 0)


class ArcSpectrum(LineSpectrum):
    """An arc spectrum

    The arc spectrum consists of emission lines from potentially multiple
    ionic species, which are read from a file.

    Parameters
    ----------
    species : `list` of `str`, or `None`
        Ionic species to select from the line list, or `None` for all.
    scale : `float`
        Scaling to provide to the fluxes in the line list.
    """
    def __init__(self, species=None, scale=1.0):
        self.species = species
        self.scale = scale

        dataRoot = os.environ.get('DRP_INSTDATA_DIR', '.')
        filepath = os.path.join(dataRoot, 'data', 'lines', 'nist_all.txt')

        self.lines = numpy.genfromtxt(filepath, usecols=list(range(3)),
                                      dtype=[('wavelength', 'f4'),
                                             ('name', 'U5',),
                                             ('flux', 'f4')])
        if species is not None:
            select = numpy.zeros(len(self.lines), dtype=bool)
            for ss in species:
                select |= self.lines['name'] == ss
            self.lines = self.lines[select]

        # Flux scale in the file is somewhat arbitrary; convert it to something of order 1.
        super().__init__(self.lines["wavelength"], scale*self.lines["flux"]*1.0e-3)

    def __str__(self):
        return("ArcSpectrum(lampset=%s, scale=%s, nlines=%d, waverange=(%g,%g))" %
               (self.species, self.scale, len(self.lines),
                self.wavelength.min(), self.wavelength.max()))


class CombSpectrum(LineSpectrum):
    """Emission line spectrum with a regular spacing in wavelength

    Parameters
    ----------
    spacing : `float`
        Spacing in wavelength between lines.
    scale : `float`
        Scale of the emission lines.
    lower : `float`
        Lower wavelength bound of the emission lines.
    upper : `float`
        Upper wavelength bound of the emission lines.
    offset : `float`
        Offset from the lower bound for the first line.
    """
    def __init__(self, spacing=50, scale=1.0, lower=300.0, upper=1300.0, offset=0.0):
        self.spacing = spacing
        self.scale = scale

        dw = upper - lower
        wavelength = numpy.linspace(lower + offset, lower + dw, dw/self.spacing + 1)
        flux = scale*numpy.ones_like(wavelength)
        super().__init__(wavelength, flux)


class TextSpectrum(TableSpectrum):
    """A spectrum read from a text file via ``numpy.genfromtxt``

    If no additional parameters are supplied, we'll guess that the file is
    comprised of wavelength and flux columns, in that order.

    Parameters
    ----------
    filename : `str`
        Path to file to read.
    wavelengthScale : `float`, optional
        Scale by which to multiply wavelengths to yield nm.
    fluxScale : `float`, optional
        Scale by which to multiply fluxes to yield nJy.
    **kwargs : `dict`, optional
        Keyword arguments for ``numpy.genfromtxt``.
    """
    def __init__(self, filename, wavelengthScale=1.0, fluxScale=1.0, **kwargs):
        if not kwargs:
            kwargs["dtype"] = [('wavelength', 'f4'), ('flux', 'f4')]
        data = numpy.genfromtxt(filename, **kwargs)
        super().__init__(data['wavelength']*wavelengthScale, data['flux']*fluxScale)


class ConstantSpectrum(FunctionalSpectrum):
    """A spectrum that is constant everywhere

    Parameters
    ----------
    value : `float`
        Constant value of spectrum.
    """
    def __init__(self, value=1.0):
        self.value = value

    def interpolate(self, wavelength):
        return self.value*numpy.ones_like(wavelength)

    def _integrateImpl(self, lower, upper):
        return self.value*(upper - lower)

    def __imul__(self, value):
        self.value *= value
        return self

    def hasFlux(self):
        return self.value > 0


class NullSpectrum(FunctionalSpectrum):
    """A spectrum that is zero everywhere"""
    def interpolate(self, wavelength):
        return numpy.zeros_like(wavelength)

    def _integrateImpl(self, lower, upper):
        return 0.0

    def integrate(self, lower, upper):
        """Integrate the spectrum between multiple wavelength bounds

        Parameters
        ----------
        lower : array_like
            Lower wavelength bounds for the integration, nm.
        upper : array_like
            Upper wavelength bounds for the integration, nm.

        Returns
        -------
        flux : array_like
            Integrated fluxes between the wavelength bounds, W/m^2.
        """
        return numpy.zeros_like(lower)

    def __imul__(self, value):
        # No change!
        return self

    def __mul__(self, other):
        return NullSpectrum()

    def __rmul__(self, other):
        return NullSpectrum()

    def hasFlux(self):
        return False


class SumSpectrum(Spectrum):
    """A spectrum that is the sum of one or more other spectra"""
    def __init__(self, *spectra):
        super().__init__()
        self.spectra = spectra

    def interpolate(self, wavelength):
        """Interpolate the spectrum at the nominated wavelength

        Parameters
        ----------
        wavelength : array_like
            Vector of wavelengths at which to interpolate, nm.

        Returns
        -------
        flux : array_like
            Vector of flux densities at the provided ``wavelength``s.
        """
        result = self.spectra[0].interpolate(wavelength)
        for ss in self.spectra[1:]:
            result += ss.interpolate(wavelength)
        return result

    def integrate(self, lower, upper):
        """Integrate the spectrum between multiple wavelength bounds

        Parameters
        ----------
        lower : array_like
            Lower wavelength bounds for the integration, nm.
        upper : array_like
            Upper wavelength bounds for the integration, nm.

        Returns
        -------
        flux : array_like
            Integrated fluxes between the wavelength bounds, W/m^2.
        """
        result = self.spectra[0].integrate(lower, upper)
        for ss in self.spectra[1:]:
            result += ss.integrate(lower, upper)
        return result

    def __imul__(self, value):
        for ss in self.spectra:
            ss *= value
        return self

    def hasFlux(self):
        return any(ss.hasFlux() for ss in self.spectra)


class ProductSpectrum(Spectrum):
    """A spectrum that is the product of one or more spectra

    Possibly useful for convolving a spectrum by a filter curve.
    """
    def __init__(self, *spectra):
        super().__init__()
        self.spectra = spectra

    def interpolate(self, wavelength):
        """Interpolate the spectrum at the nominated wavelength

        Parameters
        ----------
        wavelength : array_like
            Vector of wavelengths at which to interpolate, nm.

        Returns
        -------
        flux : array_like
            Vector of flux densities at the provided ``wavelength``s.
        """
        result = self.spectra[0].interpolate(wavelength)
        for ss in self.spectra[1:]:
            result *= ss.interpolate(wavelength)
        return result

    def hasFlux(self):
        return all(ss.hasFlux() for ss in self.spectra)


class PfsSimSpectrum(TableSpectrum):
    """A spectrum read from a pfsSim file

    Parameters
    ----------
    filename : `str`
        Path to file to read.
    wavelengthScale : `float`, optional
        Scale by which to multiply wavelengths to yield nm.
    fluxScale : `float`, optional
        Scale by which to multiply fluxes to yield nJy.
    """
    def __init__(self, filename, wavelengthScale=1.0, fluxScale=1.0):
        with astropy.io.fits.open(filename) as ff:
            flux = ff[1].data
            wavelength = ff[2].data
        super().__init__(wavelength*wavelengthScale, flux*fluxScale)


class PhotonCounting(Spectrum):
    """A photon-counting spectrum (1/nu)

    Useful for converting a filter curve to a photon-counting filter curve.

    This has no normalisation, because you're using this to set the
    normalisation.
    """
    def interpolate(self, wavelength):
        """Interpolate the spectrum at the nominated wavelength

        Parameters
        ----------
        wavelength : array_like
            Vector of wavelengths at which to interpolate, nm.

        Returns
        -------
        flux : array_like
            Vector of flux densities at the provided ``wavelength``s.
        """
        return self.interpolateFrequency(SPEED_OF_LIGHT/wavelength)

    def interpolateFrequency(self, frequency):
        """Interpolate the spectrum at the nominated freqency

        Parameters
        ----------
        frequency : array_like
            Vector of frequencies at which to interpolate, Hz.

        Returns
        -------
        flux : array_like
            Vector of flux densities at the provided ``frequency``s.
        """
        return 1.0/frequency
