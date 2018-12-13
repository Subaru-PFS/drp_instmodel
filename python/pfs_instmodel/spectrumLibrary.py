import os
import math

from .makePfsConfig import CatalogId, Lamps
from .spectrum import ArcSpectrum
from .spectrum import FlatSpectrum
from .spectrum import SlopeSpectrum
from .spectrum import CombSpectrum
from .spectrum import TextSpectrum
from .spectrum import ConstantSpectrum
from .spectrum import PfsSimSpectrum

PLANCK = 6.63e-34  # J.s
SPEED_OF_LIGHT = 3.0e8*1.0e9  # nm/s


def fluxForPhotons(photons=1.0, aperture=8.2, wavelength=600):
    """Return the flux in W/m^2 that produces the desired photon flux

    Neglecting the issue of throughput.

    Parameters
    ----------
    photons : `float`
        Photon flux (photons/sec).
    aperture : `float`
        Effective aperture (m).
    wavelength : `float`
        Wavelength (nm).

    Returns
    -------
    flux : `float`
        Flux (W/m^2).
    """
    area = math.pi*(0.5*aperture)**2  # m^2
    return photons/area/wavelength*PLANCK*SPEED_OF_LIGHT


def fluxDensityForPhotons(photons=1.0, aperture=8.2, resolution=8000):
    """Return the flux density in nJy that produces the desired photon flux

    Neglecting the issue of throughput.

    Parameters
    ----------
    photons : `float`
        Photon flux (photons/sec).
    aperture : `float`
        Effective aperture (m).
    resolution : `float`
        Resolving power (lambda/deltaLambda). Needn't be the instrument's
        resolving power; just the resolving power of the element that will
        contain the photon flux (which might be a pixel).

    Returns
    -------
    flux : `float`
        Flux (nJy).
    """
    area = math.pi*(0.5*aperture)**2  # m^2
    scale = 0.015  # photons/s/m^2 per resolution for 1 nJy
    return photons/scale/area*resolution


class SpectrumLibrary:
    def __init__(self, skyModel, skySwindle):
        self.skyModel = skyModel
        self.skySwindle = skySwindle

    def getSpectrum(self, catId, objId):
        if self.skySwindle and catId == CatalogId.SKY:
            # No sky spectrum, just noise
            return self.getNullSpectrum(0)
        return {
            CatalogId.ARC: self.getArcSpectrum,
            CatalogId.QUARTZ: self.getQuartzSpectrum,
            CatalogId.FLUXSTD: self.getFluxStdSpectrum,
            CatalogId.SKY: self.getSkySpectrum,
            CatalogId.NULL: self.getNullSpectrum,
            CatalogId.COMB: self.getCombSpectrum,
            CatalogId.SCIENCE: self.getScienceSpectrum,
        }[catId](objId)

    def getArcSpectrum(self, objId):
        lamps = []
        if objId & Lamps.NE:
            lamps.append("NeI")
        if objId & Lamps.HG:
            lamps.append("HgI")
        if objId & Lamps.XE:
            lamps.append("XeI")
        if objId & Lamps.CD:
            lamps.append("CdI")
        if objId & Lamps.KR:
            lamps.append("KrI")
        return ArcSpectrum(lamps, scale=fluxForPhotons(10000.0))

    def getQuartzSpectrum(self, objId):
        assert objId == 0, "Only one type of quartz spectrum"
        return FlatSpectrum(scale=fluxDensityForPhotons(10000.0))

    def getNullSpectrum(self, objId):
        return ConstantSpectrum(float(objId)*fluxDensityForPhotons(1.0))

    def getFluxStdSpectrum(self, objId):
        assert objId == 0, "Currently only one type of fluxStd spectrum"
        return SlopeSpectrum(scale=fluxDensityForPhotons(1.0))

    def getSkySpectrum(self, objId):
        assert objId == 0, "Currently only one type of sky spectrum"
        return self.skyModel.getSkyAt()

    def getCombSpectrum(self, objId):
        return CombSpectrum(spacing=float(objId), scale=fluxForPhotons(10000.0))

    def getScienceSpectrum(self, objId):
        # XXX ignoring objId for now
        filename = os.path.join(os.environ.get("DRP_INSTDATA_DIR", "."),
                                "data", "objects", "pfsSimObject-00000-0,0-0-00001496.fits")
        return PfsSimSpectrum(filename)
