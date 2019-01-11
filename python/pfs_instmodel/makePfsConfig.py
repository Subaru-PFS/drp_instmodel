import enum
import numpy as np

import lsst.afw.geom

from pfs.datamodel.pfsConfig import PfsConfig, TargetType


@enum.unique
class CatalogId(enum.IntEnum):
    """Catalog identifier, which controls what pool the spectra come from"""
    ARC = 1  # Arc spectrum
    QUARTZ = 2  # Quartz spectrum
    FLUXSTD = 3  # Flux standard spectrum
    SKY = 4  # Sky spectrum
    NULL = 5  # Blocked and broken
    COMB = 6  # Comb of emission lines
    SCIENCE = 7  # Science spectrum


class Lamps(enum.IntFlag):
    NONE = 0x00  # No lamps
    NE = 0x01  # Neon
    HG = 0x02  # Mercury
    XE = 0x04  # Xenon
    CD = 0x08  # Cadmium
    KR = 0x10  # Krypton


def makePfsConfig(pfiDesignId, expId, fiberIds, catIds, objIds, targetTypes,
                  fiberMags=None, filterNames=None, raBoresight=0.0*lsst.afw.geom.degrees,
                  decBoresight=0.0*lsst.afw.geom.degrees, rng=None):
    """Build a ``PfsConfig``

    The top-end settings (``ra``, ``dec``, ``pfiNominal``, ``pfiCenter``) are
    random, but everything else is set up as normal.

    Parameters
    ----------
    pfiDesignId : `int`
        Identifier for the top-end design. For our purposes, this is just a
        unique integer.
    expId : `int`
        Identifier for the exposure. For our purposes, this is just a unique
        integer.
    fiberIds : `numpy.ndarray` of `int`
        Array of identifiers for fibers that will be lit.
    catIds : `numpy.ndarray` of `int`
        Array of identifiers for catalogs from which the spectra will be taken.
    objIds : `numpy.ndarray` of `int`
        Array of identifiers for objects, which link to the particular spectra.
    targetTypes : `numpy.ndarray` of `int`
        Array of `pfs.datamodel.TargetType` enumerated values.
    fiberMags : `list` of `numpy.ndarray` of `float`
        List of magnitudes for each fiber.
    filterNames : `list` of `list` of `str`
        List of filter names for each fiber.
    raBoresight : `lsst.afw.geom.Angle`
        Right Ascension of the boresight.
    decBoresight : `lsst.afw.geom.Angle`
        Declination of the boresight.
    rng : `numpy.random.RandomState`, optional
        Random number generator. If not specified, we use the default from
        ``numpy``, which has a non-deterministic seed.

    Returns
    -------
    config : `pfs.datamodel.PfsConfig`
        Configuration of the top-end.
    """
    FIELD_OF_VIEW = 1.5*lsst.afw.geom.degrees
    PFI_SCALE = 800000.0/FIELD_OF_VIEW.asDegrees()  # microns/degree; guess, but not currently important
    PFI_ERRORS = 10  # microns
    if rng is None:
        rng = np.random
    tract = np.zeros_like(fiberIds, dtype=int)
    patch = ["0,0" for _ in fiberIds]

    num = len(fiberIds)
    boresight = lsst.afw.geom.SpherePoint(raBoresight, decBoresight)
    radius = np.sqrt(rng.uniform(size=num))*0.5*FIELD_OF_VIEW.asDegrees()  # degrees
    theta = rng.uniform(size=num)*2*np.pi  # radians
    coords = [boresight.offset(tt*lsst.afw.geom.radians, rr*lsst.afw.geom.degrees) for
              rr, tt in zip(radius, theta)]
    ra = np.array([cc.getRa().asDegrees() for cc in coords])
    dec = np.array([cc.getDec().asDegrees() for cc in coords])
    pfiNominal = (PFI_SCALE*np.array([(rr*np.cos(tt), rr*np.sin(tt)) for
                                      rr, tt in zip(radius, theta)])).astype(np.float32)
    pfiCenter = (pfiNominal + rng.normal(scale=PFI_ERRORS, size=(len(fiberIds), 2))).astype(np.float32)

    if fiberMags is None:
        fiberMags = [[] for _ in fiberIds]
    if filterNames is None:
        filterNames = [[] for _ in fiberIds]

    return PfsConfig(pfiDesignId, expId, raBoresight.asDegrees(), decBoresight.asDegrees(),
                     fiberIds, tract, patch, ra, dec, catIds, objIds, targetTypes,
                     fiberMags, filterNames, pfiCenter, pfiNominal)


def makeArcConfig(pfiDesignId, expId, arcLamps, fiberIds, rng=None):
    """Build a ``PfsConfig`` for an arc

    Parameters
    ----------
    pfiDesignId : `int`
        Identifier for the top-end design. For our purposes, this is just a
        unique integer.
    expId : `int`
        Identifier for the exposure. For our purposes, this is just a unique
        integer.
    arcLamps : `Lamps` or `int`
        Flag indicating which lamps are on, e.g., `Lamps.NE | Lamps.HG`.
    fiberIds : `numpy.ndarray` of `int`
        Array of identifiers for fibers that will be lit.
    rng : `numpy.random.RandomState`, optional
        Random number generator. If not specified, we use the default from
        ``numpy``, which has a non-deterministic seed.

    Returns
    -------
    config : `pfs.datamodel.PfsConfig`
        Configuration of the top-end.
    """
    catIds = int(CatalogId.ARC)*np.ones_like(fiberIds, dtype=int)
    objIds = arcLamps*np.ones_like(fiberIds, dtype=int)
    targetTypes = int(TargetType.SCIENCE)*np.ones_like(fiberIds, dtype=int)
    config = makePfsConfig(pfiDesignId, expId, fiberIds, catIds, objIds, targetTypes)
    return config


def makeFlatConfig(pfiDesignId, expId, fiberIds, rng=None):
    """Build a ``PfsConfig`` for a flat

    Parameters
    ----------
    pfiDesignId : `int`
        Identifier for the top-end design. For our purposes, this is just a
        unique integer.
    expId : `int`
        Identifier for the exposure. For our purposes, this is just a unique
        integer.
    fiberIds : `numpy.ndarray` of `int`
        Array of identifiers for fibers that will be lit.
    rng : `numpy.random.RandomState`, optional
        Random number generator. If not specified, we use the default from
        ``numpy``, which has a non-deterministic seed.

    Returns
    -------
    config : `pfs.datamodel.PfsConfig`
        Configuration of the top-end.
    """
    catIds = int(CatalogId.QUARTZ)*np.ones_like(fiberIds, dtype=int)
    objIds = np.zeros_like(fiberIds, dtype=int)
    targetTypes = int(TargetType.SCIENCE)*np.ones_like(fiberIds, dtype=int)
    config = makePfsConfig(pfiDesignId, expId, fiberIds, catIds, objIds, targetTypes)
    return config


def makeCombConfig(pfiDesignId, expId, fiberIds, spacing=50, rng=None):
    """Build a ``PfsConfig`` for a comb

    Parameters
    ----------
    pfiDesignId : `int`
        Identifier for the top-end design. For our purposes, this is just a
        unique integer.
    expId : `int`
        Identifier for the exposure. For our purposes, this is just a unique
        integer.
    fiberIds : `numpy.ndarray` of `int`
        Array of identifiers for fibers that will be lit.
    spacing : `int`
        Spacing between teeth in the comb. (An integer so it can be stored in
        the ``objId``.)
    rng : `numpy.random.RandomState`, optional
        Random number generator. If not specified, we use the default from
        ``numpy``, which has a non-deterministic seed.

    Returns
    -------
    config : `pfs.datamodel.PfsConfig`
        Configuration of the top-end.
    """
    catIds = int(CatalogId.COMB)*np.ones_like(fiberIds, dtype=int)
    objIds = int(spacing)*np.ones_like(fiberIds, dtype=int)
    targetTypes = int(TargetType.SCIENCE)*np.ones_like(fiberIds, dtype=int)
    config = makePfsConfig(pfiDesignId, expId, fiberIds, catIds, objIds, targetTypes)
    return config


def makeConstantConfig(pfiDesignId, expId, fiberIds, value=1, rng=None):
    """Build a ``PfsConfig`` with constant spectra

    Parameters
    ----------
    pfiDesignId : `int`
        Identifier for the top-end design. For our purposes, this is just a
        unique integer.
    expId : `int`
        Identifier for the exposure. For our purposes, this is just a unique
        integer.
    fiberIds : `numpy.ndarray` of `int`
        Array of identifiers for fibers that will be lit.
    value : `int`
        Value of spectrum. (An integer so it can be stored in
        the ``objId``.)
    rng : `numpy.random.RandomState`, optional
        Random number generator. If not specified, we use the default from
        ``numpy``, which has a non-deterministic seed.

    Returns
    -------
    config : `pfs.datamodel.PfsConfig`
        Configuration of the top-end.
    """
    catIds = int(CatalogId.NULL)*np.ones_like(fiberIds, dtype=int)
    objIds = int(value)*np.ones_like(fiberIds, dtype=int)
    targetTypes = int(TargetType.SCIENCE)*np.ones_like(fiberIds, dtype=int)
    config = makePfsConfig(pfiDesignId, expId, fiberIds, catIds, objIds, targetTypes)
    return config


def makeScienceConfig(pfiDesignId, expId, fiberIds,
                      fracSky=0.2, fracFluxStd=0.1,
                      minScienceMag=18.0, maxScienceMag=24.0, fluxStdMag=18.0,
                      raBoresight=0.0*lsst.afw.geom.degrees,
                      decBoresight=0.0*lsst.afw.geom.degrees,
                      rng=None):
    """Build a ``PfsConfig`` for a science exposure

    Fibers are randomly assigned to sky, flux standards or science targets
    based on the nominated fractions.

    Science targets have a random distribution of magnitudes in the nominated
    range.

    Parameters
    ----------
    pfiDesignId : `int`
        Identifier for the top-end design. For our purposes, this is just a
        unique integer.
    expId : `int`
        Identifier for the exposure. For our purposes, this is just a unique
        integer.
    fiberIds : `numpy.ndarray` of `int`
        Array of identifiers for fibers that will be lit.
    fracSky : `float`
        Fraction of fibers that will be devoted to sky.
    fracFluxStd : `float`
        Fraction of fibers that will be devoted to flux standards.
    minScienceMag : `float`
        Minimum magnitude of science targets.
    maxScienceMag : `float`
        Maximum magnitude of science targets.
    fluxStdMag : `float`
        Magnitude of (all) flux standards.
    raBoresight : `lsst.afw.geom.Angle`
        Right Ascension of the boresight.
    decBoresight : `lsst.afw.geom.Angle`
        Declination of the boresight.
    rng : `numpy.random.RandomState`, optional
        Random number generator. If not specified, we use the default from
        ``numpy``, which has a non-deterministic seed.

    Returns
    -------
    config : `pfs.datamodel.PfsConfig`
        Configuration of the top-end.
    """
    if rng is None:
        rng = np.random

    objIds = np.ones_like(fiberIds, dtype=int)
    numFibers = len(fiberIds)
    numSky = int(fracSky*numFibers)
    numFluxStd = int(fracFluxStd*numFibers)
    numScience = numFibers - (numSky + numFluxStd)

    targetTypes = np.array([int(TargetType.SKY)]*numSky +
                           [int(TargetType.FLUXSTD)]*numFluxStd +
                           [int(TargetType.SCIENCE)]*numScience)
    rng.shuffle(targetTypes)
    assert len(targetTypes) == len(fiberIds)

    # Currently only have one type of spectrum (objId=0) for each kind of target, so makes this easy
    objIds = np.zeros_like(fiberIds, dtype=int)
    if numScience > 0:
        scienceIds = np.arange(numScience, dtype=int)
        rng.shuffle(scienceIds)
        objIds[targetTypes == TargetType.SCIENCE] = scienceIds
    catTranslation = {TargetType.SKY: CatalogId.SKY,
                      TargetType.FLUXSTD: CatalogId.FLUXSTD,
                      TargetType.BROKEN: CatalogId.NULL,
                      TargetType.BLOCKED: CatalogId.NULL,
                      TargetType.SCIENCE: CatalogId.SCIENCE}
    catIds = np.array([int(catTranslation[tt]) for tt in targetTypes])

    noMagTypes = (TargetType.SKY, TargetType.BROKEN, TargetType.BLOCKED)
    filterNames = [["i"] if tt not in noMagTypes else [] for tt in targetTypes]
    scienceMags = 22.0*np.ones(numScience, dtype=float)
    mags = np.zeros_like(fiberIds, dtype=float)
    mags[targetTypes == TargetType.SCIENCE] = scienceMags
    mags[targetTypes == TargetType.FLUXSTD] = fluxStdMag
    fiberMags = [np.array([mm]) if tt not in noMagTypes else [] for tt, mm in zip(targetTypes, mags)]

    config = makePfsConfig(pfiDesignId, expId, fiberIds, catIds, objIds, targetTypes,
                           fiberMags=fiberMags, filterNames=filterNames, raBoresight=raBoresight,
                           decBoresight=decBoresight, rng=rng)
    return config
