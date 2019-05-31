import os
import math
from enum import IntFlag
from types import SimpleNamespace

from pfs.datamodel.pfsConfig import TargetType

from .spectrum import ArcSpectrum
from .spectrum import FlatSpectrum
from .spectrum import SlopeSpectrum
from .spectrum import ConstantSpectrum
from .spectrum import NullSpectrum
from .spectrum import PfsSimSpectrum
from .spectrum import SumSpectrum


__all__ = ["fluxForPhotons", "fluxDensityForPhotons", "DomeStatus", "LightSource"]

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


class Lamps(IntFlag):
    """Lamps that are on"""
    NONE = 0x00  # No lamps on
    QUARTZ = 0x01  # Quartz lamp
    NE = 0x02  # Neon lamp
    HG = 0x04  # Mercury lamp
    XE = 0x08  # Xenon lamp
    CD = 0x10  # Cadmium lamp
    KR = 0x20  # Krypton lamp

    @classmethod
    def fromString(cls, string):
        """Construct from a string

        ``QUARTZ`` may only be specified by itself.
        The other lamps may be combined in a whitespace-separated list.

        Parameters
        ----------
        string : `str`
            Whitespace-separated list of flags.

        Returns
        -------
        self : cls
            Interpreted lamps flag.
        """
        lamps = string.upper().split()
        if "QUARTZ" in lamps and len(lamps) > 1:
            raise RuntimeError("Quartz on while arc lamps are on")
        result = cls.NONE
        for ll in lamps:
            result |= getattr(cls, ll)
        return result


class LightSource:
    """Where the light going down the fibers is coming from

    Parameters
    ----------
    domeOpen : `bool`
        Is the dome open?
    lamps : `
    skyModel : `SkyModel`
        Model of the night sky emission.
    pfsDesign : `pfs.datamodel.PfsDesign`
        Top-end design, fiber targetting.
    spectraDir : `str`, optional
        Directory for simulated science spectra. Needed only if you're going
        to get a science spectrum.
    """
    def __init__(self, domeOpen, lamps, skyModel, pfsDesign, spectraDir=None):
        self.domeOpen = domeOpen
        self.lamps = lamps
        if domeOpen and lamps:
            raise RuntimeError("Dome open but lamps on")
        self.skyModel = skyModel
        self.pfsDesign = pfsDesign
        self.spectraDir = spectraDir

    def getSpectrum(self, fiberId):
        """Get the spectrum for a fiber

        Parameters
        ----------
        fiberId : `int`
            Fiber identifier.

        Returns
        -------
        spectrum : `Spectrum`
            Spectrum from fiber.
        """
        target = self.getTargetData(fiberId)
        if target.targetType in (TargetType.BROKEN, TargetType.BLOCKED):
            return NullSpectrum()
        if self.domeOpen:
            if target.targetType == TargetType.SKY:
                # Sky is added separately
                return NullSpectrum()
            if target.targetType == TargetType.FLUXSTD:
                return self.getFluxStdSpectrum(target)
            return self.getScienceSpectrum(target)
        if self.lamps == Lamps.QUARTZ:
            return self.getFlatSpectrum()
        return self.getArcSpectrum()

    def getTargetData(self, fiberId):
        """Get target information for fiber

        Parameters
        ----------
        fiberId : `int`
            Fiber identifier.

        Returns
        -------
        index : `int`
            Index into pfsDesign.
        catId : `int`
            Catalog identifier.
        objId : `int`
            Object identifier.
        tract : `int`
            Tract identifier.
        patch : `str`
            Patch identifier.
        targetType : `pfs.datamodel.TargetType`
            Type of target.
        fiberMags : array of `float`
            Array of magnitudes.
        """
        index = self.pfsDesign.selectFiber(fiberId)
        assert len(index) == 1
        index = index[0]
        catId = self.pfsDesign.catId[index]
        objId = self.pfsDesign.objId[index]
        tract = self.pfsDesign.tract[index]
        patch = self.pfsDesign.patch[index]
        targetType = self.pfsDesign.targetType[index]
        fiberMags = self.pfsDesign.fiberMag[index]
        return SimpleNamespace(index=index, catId=catId, objId=objId, tract=tract, patch=patch,
                               targetType=targetType, fiberMags=fiberMags)

    def getSkySpectrum(self):
        """Return a sky spectrum"""
        return self.skyModel.getSkyAt()

    def getScienceSpectrum(self, target):
        """Return a science spectrum

        A ``catId`` of zero is special-cased to a constant spectrum.
        Otherwise, the simulated science spectrum is read from the appropriate
        file; no scaling is performed.

        Parameters
        ----------
        target : struct
            Target data; output from ``getTargetData``.

        Returns
        -------
        spectrum : `Spectrum`
            Spectrum of object + sky.
        """
        if target.catId == 0:
            return self.getConstantSpectrum(target)
        if self.spectraDir is None:
            raise RuntimeError("No spectraDir specified, so can't load spectrum.")
        filename = ("pfsSimObject-%05d-%s-%03d-%016x.fits" %
                    (target.tract, target.patch, target.catId, target.objId))
        return PfsSimSpectrum(os.path.join(self.spectraDir, str(target.catId), filename))

    def getConstantSpectrum(self, target):
        """Return a constant F_nu spectrum

        This is useful for testing. The spectrum is scaled to match the target
        fiberMag.

        Parameters
        ----------
        target : struct
            Target data; output from ``getTargetData``.

        Returns
        -------
        result : `ConstantSpectrum`
            Flat F_nu spectrum, scaled appropriately.
        """
        # Flat F_nu spectrum, useful for testing
        fiberMag = set(target.fiberMags)
        if len(fiberMag) != 1:
            raise RuntimeError("Insufficient or differing fiberMags for constant spectrum: %s" %
                               (fiberMag,))
        fiberMag = fiberMag.pop()
        return ConstantSpectrum(3631e9*10**(-0.4*fiberMag))

    def getFlatSpectrum(self):
        """Return a flat-field spectrum"""
        return FlatSpectrum()

    def getArcSpectrum(self):
        """Return an arc spectrum"""
        if self.lamps == Lamps.NONE:
            return NullSpectrum()
        lamps = []
        if self.lamps & Lamps.NE:
            lamps.append("NeI")
        if self.lamps & Lamps.HG:
            lamps.append("HgI")
        if self.lamps & Lamps.XE:
            lamps.append("XeI")
        if self.lamps & Lamps.CD:
            lamps.lamps("CdI")
        if self.lamps & Lamps.KR:
            lamps.lamps("KrI")
        return ArcSpectrum(lamps, scale=fluxForPhotons(10000.0))

    def getFluxStdSpectrum(self, target):
        """Return a flux standard spectrum

        XXX We haven't discussed how the pipeline will represent and store
        spectra used for flux standards, so for now let's just make this a
        constant F_nu spectrum.

        Parameters
        ----------
        target : struct
            Target data; output from ``getTargetData``.

        Returns
        -------
        result : `Spectrum`
            Flux standard spectrum
        """
        return self.getConstantSpectrum(target)
