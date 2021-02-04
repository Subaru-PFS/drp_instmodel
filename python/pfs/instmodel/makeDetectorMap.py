from argparse import ArgumentParser

from .detector import Detector
from .splinedPsf import SplinedPsf

__all__ = ("makeDetectorMap",)


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
    return SplinedPsf(detector, spotID={}).makeDetectorMap(obsdate=date)


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
