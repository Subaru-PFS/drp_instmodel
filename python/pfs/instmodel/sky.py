__all__ = ["SkyModel", "StaticSkyModel"]

import os
import numpy
import scipy.interpolate

from .spectrum import TableSpectrum


class SkyModel(object):
    """Generator for additive and multiplicative sky spectra

    The additive term is what we usually call "the sky".
    The multiplicative term is what we usually call "extinction".

    Parameters
    ----------
    band : `pfs.instmodel.arm.Arm`
        Arm for sky spectrum.
    zenithDistance : `float`
        Zenith distance for observation (degrees)
    pwv : `float`
        Precipitable water vapour (mm)
    aerosol : `float`
        Angstrom exponent for the aerosols (power-law index of wavelength for
        the optical depth).
    """

    def __init__(self, arm, zenithDistance=45.0, pwv=1.6, aerosol=1.0):
        self.arm = arm
        self.zenithDistance = zenithDistance
        self.pwv = pwv
        self.aerosol = aerosol

    def __str__(self):
        return ("%s(arm=%s, zenithDistance=%f, pwv=%f, aerosol=%f)" %
                (self.__class__.__name__, self.arm, self.zenithDistance, self.pwv, self.aerosol))

    def getSky(self):
        """Return a sky spectrum"""
        raise NotImplementedError("getSky")

    def getExtinction(self):
        """Return an extinction spectrum"""
        raise NotImplementedError("getExtinction")


class StaticSkyModel(SkyModel):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._readSky()
        self._readExtinction()

    def _readSky(self):
        """Read in JEG's preliminary sky model. """

        # Map arm to a filename using some config file. For now, hardcode
        dataRoot = os.environ.get('DRP_INSTDATA_DIR', '.')
        filepath = os.path.join(dataRoot, 'data', 'sky', 'sumire%sskyHR.dat' % (self.arm.toBand().upper()))

        data = numpy.genfromtxt(filepath, skip_footer=20, comments='\\')

        self.skyWavelength = data[:, 0]/10.0  # Convert Angstroms --> nm
        fiberDiameter = 1.1  # arcsec
        fiberArea = numpy.pi*(0.5*fiberDiameter)**2  # arcsec^2
        self.skyFlux = data[:, 2]*3.631e3*fiberArea  # Convert nanoMaggie/arcsec^2 --> nJy

    def _readExtinction(self):
        """Read in the extinction model. """
        dataRoot = os.environ.get('DRP_INSTDATA_DIR', '.')

        filenames = []
        zd = []
        aerosol = []
        pwv = []

        summaryFilename = os.path.join(dataRoot, "data", "atmosphere", "model_atmosphere_metadata.txt")
        summary = numpy.recfromtxt(summaryFilename, comments='#')
        num = len(summary)
        for line in summary:
            filenames.append(os.path.join(dataRoot, "data", "atmosphere", line[0].decode()))
            zd.append(line[1])  # Zenith distance (deg)
            aerosol.append(line[2])  # Aerosol exponent (power law index)
            assert line[3] == 0.0075  # Aerosol tau (optical depth at 1um)
            pwv.append(line[4])  # Precipitable water vapour (mm)
            assert line[5] == 616  # Pressure (hPa)
            assert line[6] == 260   # Ozone (Dobson units)
            assert line[7] == 4.2  # Altitude (km)

        zdUnique = numpy.array(sorted(set(zd)))
        aerosolUnique = numpy.array(sorted(set(aerosol)))
        pwvUnique = numpy.array(sorted(set(pwv)))

        # Simple case: one matches exactly
        if self.zenithDistance in zdUnique and self.pwv in pwvUnique and self.aerosol in aerosolUnique:
            index = numpy.where((numpy.array(zd) == self.zenithDistance) &
                                (numpy.array(aerosol) == self.aerosol) &
                                (numpy.array(pwv) == self.pwv))[0]
            assert len(index) == 1
            self.extWavelength, self.extTransparency = numpy.loadtxt(filenames[index[0]]).T
            return

        # More complicated: we read all the data, and interpolate
        wavelength = []
        transparency = []
        for fn in filenames:
            ww, tt = numpy.loadtxt(fn).T
            wavelength.append(ww)  # Wavelength (nm)
            transparency.append(tt)  # Transparency

        wlUnique = numpy.array(wavelength[0])
        for ww in wavelength[1:]:
            assert numpy.all(numpy.array(ww) == wlUnique)
        self.extWavelength = numpy.array(wavelength[0])
        transparency = numpy.array(transparency)

        zdNum = len(zdUnique)
        aerosolNum = len(aerosolUnique)
        pwvNum = len(pwvUnique)
        assert zdNum*aerosolNum*pwvNum == num

        gridShape = (zdNum, aerosolNum, pwvNum)
        indices = numpy.empty(gridShape, dtype=int)
        for ii in range(num):
            zdIndex = numpy.where(zdUnique == zd[ii])[0]
            aerosolIndex = numpy.where(aerosolUnique == aerosol[ii])[0]
            pwvIndex = numpy.where(pwvUnique == pwv[ii])[0]
            indices[pwvIndex, aerosolIndex, zdIndex] = ii
        indices = numpy.reshape(indices, -1)

        xx = (zdUnique, aerosolUnique, pwvUnique)
        zz = numpy.array([self.zenithDistance, self.aerosol, self.pwv])
        self.extTransparency = numpy.empty_like(self.extWavelength)
        for ii, wl in enumerate(self.extWavelength):
            yy = numpy.reshape(transparency[indices, ii], gridShape)
            self.extTransparency[ii] = scipy.interpolate.interpn(xx, yy, zz)

    def getSky(self):
        """Return a spline for the sky"""
        return TableSpectrum(self.skyWavelength, self.skyFlux)

    def getExtinction(self):
        return TableSpectrum(self.extWavelength, self.extTransparency)
