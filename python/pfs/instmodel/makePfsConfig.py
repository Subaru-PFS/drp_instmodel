import numpy as np
import astropy.units as u

import lsst.afw.geom
import lsst.log

from pfs.datamodel.pfsConfig import PfsDesign, TargetType, FiberStatus, PfsConfig

FLUXSTD_MAG = 18.0  # ABmag

logger = lsst.log.Log.getLogger("pfs.instmodel.makePfsConfig")


def makePfsConfig(pfsDesign, expId, rng=None, pfiErrors=10.0):
    """Build a ``PfsConfig`` from a ``PfsDesign``

    The ``PfsDesign`` supplies just about everything we need, except for the
    ``pfiCenter``, which will be a random modification of the ``pfiNominal``.

    Parameters
    ----------
    pfsDesign : `pfs.datamodel.PfsDesign`
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
    centering = rng.normal(scale=pfiErrors, size=pfsDesign.pfiNominal.shape)
    pfiCenter = pfsDesign.pfiNominal + centering
    return PfsConfig.fromPfsDesign(pfsDesign, expId, pfiCenter)


def makePfsDesign(pfsDesignId, fiberIds, catIds, objIds, targetTypes,
                  fiberStatus,
                  fiberFlux=None,
                  psfFlux=None,
                  totalFlux=None,
                  fiberFluxErr=None, psfFluxErr=None, totalFluxErr=None,
                  filterNames=None, raBoresight=0.0*lsst.afw.geom.degrees,
                  decBoresight=0.0*lsst.afw.geom.degrees,
                  posAng=0.0*lsst.afw.geom.degrees,
                  rng=None):
    """Build a ``PfsDesign``

    The top-end settings (``ra``, ``dec``, ``pfiNominal``, ``pfiCenter``) are
    random, but everything else is set up as normal.

    Parameters
    ----------
    pfsDesignId : `int`
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
    fiberStatus : `numpy.ndarray` of `int`
        Array of `pfs.datamodel.FiberStatus` enumerated values.
    fiberFlux : `list` of `numpy.ndarray` of `float`
        List of fiber fluxes for each fiber [nJy].
    psfFlux : `list` of `numpy.ndarray` of `float`
        List of psf fluxes for each fiber [nJy].
    totalFlux : `list` of `numpy.ndarray` of `float`
        List of total fluxes for each fiber [nJy].
    fiberFluxErr : `list` of `numpy.ndarray` of `float`
        List of fiber flux errors for each fiber [nJy].
    psfFluxErr : `list` of `numpy.ndarray` of `float`
        List of psf flux errors for each fiber [nJy].
    totalFluxErr : `list` of `numpy.ndarray` of `float`
        List of total flux errors for each fiber [nJy].
    filterNames : `list` of `list` of `str`
        List of filter names for each fiber.
    raBoresight : `lsst.afw.geom.Angle`
        Right Ascension of the boresight.
    decBoresight : `lsst.afw.geom.Angle`
        Declination of the boresight.
    posAng : `lsst.afw.geom.Angle`
        position angle of the PFI.
    rng : `numpy.random.RandomState`, optional
        Random number generator. If not specified, we use the default from
        ``numpy``, which has a non-deterministic seed.

    Returns
    -------
    design : `pfs.datamodel.PfsDesign`
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

    if fiberFlux is None:
        fiberFlux = [[] for _ in fiberIds]
    if psfFlux is None:
        psfFlux = [[] for _ in fiberIds]
    if totalFlux is None:
        totalFlux = [[] for _ in fiberIds]
    if fiberFluxErr is None:
        fiberFluxErr = [[] for _ in fiberIds]
    if psfFluxErr is None:
        psfFluxErr = [[] for _ in fiberIds]
    if totalFluxErr is None:
        totalFluxErr = [[] for _ in fiberIds]
    if filterNames is None:
        filterNames = [[] for _ in fiberIds]

    return PfsDesign(pfsDesignId, raBoresight.asDegrees(),
                     decBoresight.asDegrees(),
                     posAng.asDegrees(),
                     fiberIds, tract, patch, ra, dec, catIds, objIds,
                     targetTypes, fiberStatus,
                     fiberFlux, psfFlux, totalFlux,
                     fiberFluxErr, psfFluxErr, totalFluxErr,
                     filterNames, pfiNominal)


def makeScienceDesign(pfsDesignId, fiberIds,
                      fracSky=0.2, fracFluxStd=0.1,
                      minScienceMag=18.0, maxScienceMag=24.0,
                      fluxStdMag=18.0,
                      scienceCatId=0, scienceObjId=None,
                      raBoresight=0.0*lsst.afw.geom.degrees,
                      decBoresight=0.0*lsst.afw.geom.degrees,
                      posAng=0.0*lsst.afw.geom.degrees,
                      rng=None, unlitFiberIds=None):
    """Build a ``PfsDesign`` for a science exposure

    Fibers are randomly assigned to sky, flux standards or science targets
    based on the nominated fractions.

    Science targets have a random distribution of magnitudes in the nominated
    range.

    Parameters
    ----------
    pfsDesignId : `int`
        Identifier for the top-end design. For our purposes, this is just a
        unique integer.
    fiberIds : `numpy.ndarray` of `int`
        Array of identifiers for fibers that will be lit.
    fracSky : `float`
        Fraction of fibers that will be devoted to sky.
    fracFluxStd : `float`
        Fraction of fibers that will be devoted to flux standards.
    minScienceMag : `float`
        Minimum AB magnitude of science targets.
    maxScienceMag : `float`
        Maximum AB magnitude of science targets.
    scienceCatId : `int`, optional
        Catalog identifier for science targets.
    scienceObjId : iterable of `int`, optional
        Object identifiers for science targets.
    fluxStdMag : `float`
        AB Magnitude of flux standard targets.
    raBoresight : `lsst.afw.geom.Angle`
        Right Ascension of the boresight.
    decBoresight : `lsst.afw.geom.Angle`
        Declination of the boresight.
    posAng : `lsst.afw.geom.Angle`
        position angle of the PFI.
    rng : `numpy.random.RandomState`, optional
        Random number generator. If not specified, we use the default from
        ``numpy``, which has a non-deterministic seed.
    unlitFiberIds : `numpy.ndarray` of `int`
        Array of identifiers for fibers that will not be lit. These will be
        flagged as blocked.

    Returns
    -------
    design : `pfs.datamodel.PfsDesign`
        Design of the top-end.
    """
    if rng is None:
        rng = np.random
    if unlitFiberIds is None:
        unlitFiberIds = np.array([], dtype=int)

    allFiberIds = np.concatenate((fiberIds, unlitFiberIds))

    numLit = len(fiberIds)
    numUnlit = len(unlitFiberIds)
    numFibers = numLit + numUnlit
    numSky = int(fracSky*numLit)
    numFluxStd = int(fracFluxStd*numLit)
    numScience = numLit - (numSky + numFluxStd)

    targetTypes = np.array([int(TargetType.SKY)]*numSky +
                           [int(TargetType.FLUXSTD)]*numFluxStd +
                           [int(TargetType.SCIENCE)]*numScience +
                           [int(TargetType.UNASSIGNED)]*numUnlit)
    fiberStatus = np.array([int(FiberStatus.GOOD)]*numLit +
                           [int(FiberStatus.UNILLUMINATED)]*numUnlit)
    rng.shuffle(targetTypes[:numLit])
    assert len(targetTypes) == numFibers
    assert len(fiberStatus) == numFibers

    objId = np.zeros(numFibers, dtype=int)  # Object ID for all fibers
    catId = np.zeros(numFibers, dtype=int)  # Catalog ID for all fibers
    objIdStart = 1
    if numFluxStd > 0:
        catId[targetTypes == TargetType.FLUXSTD] = 0
        fluxStdObjId = np.arange(numFluxStd, dtype=int) + objIdStart
        objIdStart += numFluxStd
        rng.shuffle(fluxStdObjId)
        objId[targetTypes == TargetType.FLUXSTD] = fluxStdObjId

    if numSky > 0:
        catId[targetTypes == TargetType.SKY] = 0
        skyObjId = np.arange(numSky, dtype=int) + objIdStart
        objIdStart += numSky
        rng.shuffle(skyObjId)
        objId[targetTypes == TargetType.SKY] = skyObjId

    if numScience > 0:
        if scienceObjId is None:
            scienceObjId = np.arange(numScience, dtype=int)
            if scienceCatId == 0:
                scienceObjId += objIdStart
                objIdStart += numScience
        scienceObjId = rng.choice(scienceObjId, numScience, len(scienceObjId) < numScience)
        objId[targetTypes == TargetType.SCIENCE] = scienceObjId
        catId[targetTypes == TargetType.SCIENCE] = scienceCatId

    if numUnlit > 0:
        catId[fiberStatus == FiberStatus.UNILLUMINATED] = -1
        objId[fiberStatus == FiberStatus.UNILLUMINATED] = -1

    createMags = (((targetTypes == TargetType.SCIENCE) ^ (targetTypes == TargetType.FLUXSTD)) &
                  (fiberStatus == FiberStatus.GOOD))

    filterNames = [["i"] if xx else [] for xx in createMags]
    scienceMags = rng.uniform(minScienceMag, maxScienceMag, numScience)
    mags = np.zeros(numFibers, dtype=float)
    mags[targetTypes == TargetType.SCIENCE] = scienceMags
    mags[targetTypes == TargetType.FLUXSTD] = fluxStdMag

    fluxes = [(mm*u.ABmag).to_value(u.nJy) for mm in mags]
    fiberFlux = [np.array([fx]) if xx else[]
                 for xx, fx in zip(createMags, fluxes)]
    psfFlux = fiberFlux.copy()
    totalFlux = fiberFlux.copy()

    # Assigning 1% of flux as error
    fiberFluxErr = [fFlux * 0.01 if fFlux else [] for fFlux in fiberFlux]
    psfFluxErr = fiberFluxErr.copy()
    totalFluxErr = fiberFluxErr.copy()

    indices = np.argsort(allFiberIds)
    return makePfsDesign(pfsDesignId, allFiberIds[indices], catId[indices], objId[indices],
                         targetTypes[indices], fiberStatus[indices],
                         fiberFlux=[fiberFlux[ii] for ii in indices],
                         psfFlux=[psfFlux[ii] for ii in indices],
                         totalFlux=[totalFlux[ii] for ii in indices],
                         fiberFluxErr=[fiberFluxErr[ii] for ii in indices],
                         psfFluxErr=[psfFluxErr[ii] for ii in indices],
                         totalFluxErr=[totalFluxErr[ii] for ii in indices],
                         filterNames=[filterNames[ii] for ii in indices],
                         raBoresight=raBoresight, decBoresight=decBoresight,
                         posAng=posAng,
                         rng=rng)
