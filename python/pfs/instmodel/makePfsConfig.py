import enum
import numpy as np

import lsst.afw.geom

from pfs.datamodel.pfsConfig import PfiDesign, TargetType, PfsConfig

FLUXSTD_MAG = 18.0  # ABmag


def makePfsConfig(pfiDesign, expId, rng=None, pfiErrors=10.0):
    """Build a ``PfsConfig`` from a ``PfiDesign``

    The ``PfiDesign`` supplies just about everything we need, except for the
    ``pfiCenter``, which will be a random modification of the ``pfiNominal``.

    Parameters
    ----------
    pfiDesign : `pfs.datamodel.PfiDesign`
        Design for the top-end.
    expId : `int`
        Identifier for the exposure. For our purposes, this is just a unique
        integer.
    rng : `numpy.random.RandomState`, optional
        Random number generator. If not specified, we use the default from
        ``numpy``, which has a non-deterministic seed.
    pfiErrors : `float`
        Standard deviation of cobra positioning errors, microns.
    """
    if rng is None:
        rng = np.random
    centering = rng.normal(scale=pfiErrors, size=pfiDesign.pfiNominal.shape)
    pfiCenter = pfiDesign.pfiNominal + centering
    return PfsConfig.fromPfiDesign(pfiDesign, expId, pfiCenter)


def makePfiDesign(pfiDesignId, fiberIds, catIds, objIds, targetTypes,
                  fiberMags=None, filterNames=None, raBoresight=0.0*lsst.afw.geom.degrees,
                  decBoresight=0.0*lsst.afw.geom.degrees, rng=None):
    """Build a ``PfiDesign``

    The top-end settings (``ra``, ``dec``, ``pfiNominal``, ``pfiCenter``) are
    random, but everything else is set up as normal.

    Parameters
    ----------
    pfiDesignId : `int`
        Identifier for the top-end design. For our purposes, this is just a
        unique integer.
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
    design : `pfs.datamodel.PfiDesign`
        Design of the top-end.
    """
    FIELD_OF_VIEW = 1.5*lsst.afw.geom.degrees
    PFI_SCALE = 800000.0/FIELD_OF_VIEW.asDegrees()  # microns/degree; guess, but not currently important
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

    if fiberMags is None:
        fiberMags = [[] for _ in fiberIds]
    if filterNames is None:
        filterNames = [[] for _ in fiberIds]

    return PfiDesign(pfiDesignId, raBoresight.asDegrees(), decBoresight.asDegrees(),
                     fiberIds, tract, patch, ra, dec, catIds, objIds, targetTypes,
                     fiberMags, filterNames, pfiNominal)


def makeScienceDesign(pfiDesignId, fiberIds,
                      fracSky=0.2, fracFluxStd=0.1,
                      minScienceMag=18.0, maxScienceMag=24.0,
                      fluxStdMag=18.0,
                      scienceCatId=0, scienceObjId=None,
                      raBoresight=0.0*lsst.afw.geom.degrees,
                      decBoresight=0.0*lsst.afw.geom.degrees,
                      rng=None):
    """Build a ``PfiDesign`` for a science exposure

    Fibers are randomly assigned to sky, flux standards or science targets
    based on the nominated fractions.

    Science targets have a random distribution of magnitudes in the nominated
    range.

    Parameters
    ----------
    pfiDesignId : `int`
        Identifier for the top-end design. For our purposes, this is just a
        unique integer.
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
    scienceCatId : `int`, optional
        Catalog identifier for science targets.
    scienceObjId : iterable of `int`, optional
        Object identifiers for science targets.
    fluxStdMag : `float`
        Magnitude of flux standard targets.
    raBoresight : `lsst.afw.geom.Angle`
        Right Ascension of the boresight.
    decBoresight : `lsst.afw.geom.Angle`
        Declination of the boresight.
    rng : `numpy.random.RandomState`, optional
        Random number generator. If not specified, we use the default from
        ``numpy``, which has a non-deterministic seed.

    Returns
    -------
    design : `pfs.datamodel.PfiDesign`
        Design of the top-end.
    """
    if rng is None:
        rng = np.random

    numFibers = len(fiberIds)
    numSky = int(fracSky*numFibers)
    numFluxStd = int(fracFluxStd*numFibers)
    numScience = numFibers - (numSky + numFluxStd)

    targetTypes = np.array([int(TargetType.SKY)]*numSky +
                           [int(TargetType.FLUXSTD)]*numFluxStd +
                           [int(TargetType.SCIENCE)]*numScience)
    rng.shuffle(targetTypes)
    assert len(targetTypes) == len(fiberIds)

    objId = np.zeros_like(fiberIds, dtype=int)  # Object ID for all fibers
    catId = np.zeros_like(objId)
    objIdStart = 1
    if numFluxStd > 0:
        catId[targetTypes == TargetType.FLUXSTD] = 0
        fluxStdObjId = np.arange(numFluxStd, dtype=int) + objIdStart
        objIdStart += numFluxStd
        rng.shuffle(fluxStdObjId)
        objId[targetTypes == TargetType.FLUXSTD] = fluxStdObjId

    if numScience > 0:
        if scienceObjId is None:
            scienceObjId = np.arange(numScience, dtype=int)
            if scienceCatId == 0:
                scienceObjId += objIdStart
                objIdStart += numScience
        rng.shuffle(scienceObjId)
        if len(scienceObjId) > numScience:
            scienceObjId = scienceObjId[:numScience]
        objId[targetTypes == TargetType.SCIENCE] = scienceObjId
        catId[targetTypes == TargetType.SCIENCE] = scienceCatId

    noMagTypes = (TargetType.SKY, TargetType.BROKEN, TargetType.BLOCKED)
    filterNames = [["i"] if tt not in noMagTypes else [] for tt in targetTypes]
    scienceMags = rng.uniform(minScienceMag, maxScienceMag, numScience)
    mags = np.zeros_like(fiberIds, dtype=float)
    mags[targetTypes == TargetType.SCIENCE] = scienceMags
    mags[targetTypes == TargetType.FLUXSTD] = fluxStdMag
    fiberMags = [np.array([mm]) if tt not in noMagTypes else [] for tt, mm in zip(targetTypes, mags)]

    return makePfiDesign(pfiDesignId, fiberIds, catId, objId, targetTypes,
                         fiberMags=fiberMags, filterNames=filterNames, raBoresight=raBoresight,
                         decBoresight=decBoresight, rng=rng)
