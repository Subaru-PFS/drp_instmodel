import os
import math
from enum import IntFlag
from types import SimpleNamespace
import collections
import numpy as np
import logging

from pfs.datamodel.pfsConfig import TargetType, FiberStatus

from .spectrum import ArcSpectrum
from .spectrum import FlatSpectrum
from .spectrum import ConstantSpectrum
from .spectrum import NullSpectrum
from .spectrum import PfsSimSpectrum
from .spectrum import AmbreSpectrum


__all__ = ["fluxForPhotons", "fluxDensityForPhotons", "LightSource"]

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

    def toFitsCards(self):
        """Generate appropriate FITS headers

        Returns
        -------
        headers : `list` of `tuple`
            FITS cards, each with elements: keyword, value, comment
        """
        mapping = {
            Lamps.QUARTZ: "W_AITQTH",
            Lamps.NE: "W_AITNEO",
            Lamps.HG: "W_AITHGA",
            Lamps.XE: "W_AITXEN",
            Lamps.KR: "W_AITKRY",
        }
        assert (self & Lamps.CD) == 0, "No keyword defined for Cd lamp"
        return [(mapping[key], (self & key) > 0, "Lamp %s is on" % key.name) for key in mapping]


class LightSource:
    """Where the light going down the fibers is coming from

    Parameters
    ----------
    domeOpen : `bool`
        Is the dome open?
    lamps : `Lamps`
        Lamps that are active.
    skyModel : `SkyModel`
        Model of the night sky emission.
    pfsDesign : `pfs.datamodel.PfsDesign`
        Top-end design, fiber targetting.
    spectraDir : `str`, optional
        Directory for simulated science spectra. Needed only if you're going
        to get a science spectrum.
    """
    def __init__(self, domeOpen, lamps, skyModel, pfsDesign, spectraDir=None,
                 logger=None):
        if logger is None:
            logger = logging.getLogger('lightSource')
        self.logger = logger

        self.domeOpen = domeOpen
        self.lamps = lamps
        if domeOpen and lamps != Lamps.NONE:
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
        if target.fiberStatus != FiberStatus.GOOD:
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
        fiberStatus : `pfs.datamodel.FiberStatus`
            Status of fiber.
        fiberFlux : array of `float`
            Array of fluxes in [nJy].
        psfFlux : array of `float`
            List of psf fluxes for each fiber [nJy].
        totalFlux : array of `float`
            List of total fluxes for each fiber [nJy].
        fiberFluxErr : array of `float`
            List of fiber flux errors for each fiber [nJy].
        psfFluxErr : array of `float`
            List of psf flux errors for each fiber [nJy].
        totalFluxErr : array of `float`
            List of total flux errors for each fiber [nJy].
        """
        index = self.pfsDesign.selectFiber(fiberId)
        if isinstance(index, (collections.Sequence, np.ndarray)):
            assert len(index) == 1
            index = index[0]
        catId = self.pfsDesign.catId[index]
        objId = self.pfsDesign.objId[index]
        tract = self.pfsDesign.tract[index]
        patch = self.pfsDesign.patch[index]
        targetType = self.pfsDesign.targetType[index]
        fiberStatus = self.pfsDesign.fiberStatus[index]
        fiberFlux = dict(zip(self.pfsDesign.filterNames[index],
                         self.pfsDesign.fiberFlux[index]))
        psfFlux = dict(zip(self.pfsDesign.filterNames[index],
                       self.pfsDesign.psfFlux[index]))
        totalFlux = dict(zip(self.pfsDesign.filterNames[index],
                         self.pfsDesign.totalFlux[index]))
        fiberFluxErr = dict(zip(self.pfsDesign.filterNames[index],
                            self.pfsDesign.fiberFluxErr[index]))
        psfFluxErr = dict(zip(self.pfsDesign.filterNames[index],
                          self.pfsDesign.psfFluxErr[index]))
        totalFluxErr = dict(zip(self.pfsDesign.filterNames[index],
                            self.pfsDesign.totalFluxErr[index]))
        return SimpleNamespace(index=index, catId=catId, objId=objId,
                               tract=tract, patch=patch,
                               targetType=targetType, fiberStatus=fiberStatus,
                               fiberFlux=fiberFlux,
                               psfFlux=psfFlux,
                               totalFlux=totalFlux,
                               fiberFluxErr=fiberFluxErr,
                               psfFluxErr=psfFluxErr,
                               totalFluxErr=totalFluxErr
                               )

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
        self.logger.info('Getting science spectrum for target '
                         f'catId={target.catId}, objId={target.objId:#0x} ..')
        if target.catId == 0:
            return self.getConstantSpectrum(target)
        if self.spectraDir is None:
            raise RuntimeError("No spectraDir specified, so can't load spectrum.")
        catMenu = {1: "lowz_COSMOS",
                   }
        catDir = catMenu.get(target.catId, str(target.catId))
        filename = ("pfsSimObject-%05d-%05d-%s-%016x.fits" %
                    (target.catId, target.tract, target.patch, target.objId))
        return PfsSimSpectrum(os.path.join(self.spectraDir, catDir, filename))

    def getConstantSpectrum(self, target):
        """Return a constant F_nu spectrum

        This is useful for testing. The spectrum is scaled to match the target
        fiberFlux.

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
        fiberFlux = set(target.fiberFlux)
        if len(fiberFlux) != 1:
            raise RuntimeError("Insufficient or differing fiberFlux for constant spectrum: %s" %
                               (fiberFlux,))
        fflux = fiberFlux.pop()  # Flux is in nJy
        return ConstantSpectrum(fflux)

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
            lamps.append("ArI")
        if self.lamps & Lamps.XE:
            lamps.append("XeI")
        if self.lamps & Lamps.CD:
            lamps.append("CdI")
        if self.lamps & Lamps.KR:
            lamps.append("KrI")
        return ArcSpectrum(lamps, scale=fluxForPhotons(10000.0))

    def getFluxStdSpectrum(self, target):
        """Return a flux standard spectrum

        Parameters
        ----------
        target : struct
            Target data; output from ``getTargetData``.

        Returns
        -------
        result : `Spectrum`
            Flux standard spectrum
        """
        menu = ["p6500_g+4.0_m0.0_t01_z+0.00_a+0.00.AMBRE_Extp.fits",
                "p6500_g+4.0_m0.0_t01_z-1.00_a+0.00.AMBRE_Extp.fits",
                "p7000_g+4.0_m0.0_t01_z+0.00_a+0.00.AMBRE_Extp.fits",
                "p7000_g+4.0_m0.0_t01_z-1.00_a+0.00.AMBRE_Extp.fits",
                "p7500_g+4.0_m0.0_t01_z+0.00_a+0.00.AMBRE_Extp.fits",
                "p7500_g+4.0_m0.0_t01_z-1.00_a+0.00.AMBRE_Extp.fits",
                ]
        filename = os.path.join(os.environ["DRP_INSTDATA_DIR"], "data", "objects", "fluxCal",
                                menu[target.objId % len(menu)])
        spectrum = AmbreSpectrum(filename)

        filterName = "i"
        spectrum.normalize(filterName, target.fiberFlux[filterName])
        return spectrum
