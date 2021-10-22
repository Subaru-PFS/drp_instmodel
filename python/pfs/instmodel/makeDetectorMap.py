from argparse import ArgumentParser

import numpy as np

from .detector import Detector
from .splinedPsf import SplinedPsf
from .slit import Slit

__all__ = ("makeDetectorMap",)


def interpolate(x, x1, x2, y1, y2):
    """Interpolate between two points

    Parameters
    ----------
    x : `float`
        X coordinate at which to interpolate.
    x1, x2 : `float`
        X coordinate of points.
    y1, y2 : array_like
        Y coordinate(s) of points.

    Returns
    -------
    y : array_like
        Interpolated y coordinate(s).
    """
    weight2 = (x - x1)/(x2 - x1)
    weight1 = 1.0 - weight2
    return y1*weight1 + y2*weight2


def makeDetectorMap(detectorName, date="2020-01-01T00:00:00"):
    """Make and write a detectorMap

    Parameters
    ----------
    detectorName : `str`
        Name of detector (arm and spectrograph, e.g., ``r1``).

    Returns
    -------
    detectorMap : `pfs.drp.stella.DetectorMap`
        Mapping of fiberId,wavelength to x,y.
    """
    detector = Detector(detectorName)
    detMap = SplinedPsf(detector, spotID={}).makeDetectorMap(obsdate=date)

    # The engineering fibers don't have spots that would allow us to include them in the detectorMap
    # directly, so we will attempt to guess where they should be based on where the other fibers are.
    numFibers = len(detMap)
    engineeringFibers = Slit(detector.spectrograph, allHoles=True).engineeringFibers
    xCenterKnots = {ff: detMap.getXCenterSpline(ff).getX() for ff in detMap.fiberId}
    xCenterValues = {ff: detMap.getXCenterSpline(ff).getY() for ff in detMap.fiberId}
    wavelengthKnots = {ff: detMap.getWavelengthSpline(ff).getX() for ff in detMap.fiberId}
    wavelengthValues = {ff: detMap.getWavelengthSpline(ff).getY() for ff in detMap.fiberId}
    for fiberId in engineeringFibers:
        index = min(np.searchsorted(detMap.fiberId, fiberId), numFibers - 2)
        low = detMap.fiberId[index]
        high = detMap.fiberId[index + 1]
        xCenterKnots[fiberId] = interpolate(fiberId, low, high, detMap.getXCenterSpline(low).getX(),
                                            detMap.getXCenterSpline(high).getX())
        xCenterValues[fiberId] = interpolate(fiberId, low, high, detMap.getXCenterSpline(low).getY(),
                                             detMap.getXCenterSpline(high).getY())
        wavelengthKnots[fiberId] = interpolate(fiberId, low, high, detMap.getWavelengthSpline(low).getX(),
                                               detMap.getWavelengthSpline(high).getX())
        wavelengthValues[fiberId] = interpolate(fiberId, low, high, detMap.getWavelengthSpline(low).getY(),
                                                detMap.getWavelengthSpline(high).getY())

    fiberId = np.concatenate((detMap.fiberId, + engineeringFibers))
    fiberId.sort()
    return type(detMap)(
        detMap.bbox, fiberId,
        [xCenterKnots[ff] for ff in fiberId],
        [xCenterValues[ff] for ff in fiberId],
        [wavelengthKnots[ff] for ff in fiberId],
        [wavelengthValues[ff] for ff in fiberId],
        None, None,
        detMap.visitInfo,
        detMap.metadata,
    )


def main():
    """CLI to make and write a detectorMap"""
    parser = ArgumentParser()
    parser.add_argument("detector", help="Detector name (arm and spectrograph, e.g., r1)")
    parser.add_argument("filename", help="Filename for output detectorMap")
    parser.add_argument("--date", default="2020-01-01T00:00:00", help="calibTime for detectorMap")
    args = parser.parse_args()
    detMap = makeDetectorMap(args.detector, args.date)
    detMap.writeFits(args.filename)
    return detMap


if __name__ == "__main__":
    main()
