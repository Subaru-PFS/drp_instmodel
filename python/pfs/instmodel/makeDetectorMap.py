from argparse import ArgumentParser

from .detector import Detector
from .splinedPsf import SplinedPsf

__all__ = ("makeDetectorMap",)


def makeDetectorMap(detectorName):
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
    return SplinedPsf(detector, spotID={}).makeDetectorMap()


def main():
    """CLI to make and write a detectorMap"""
    parser = ArgumentParser()
    parser.add_argument("detector", help="Detector name (arm and spectrograph, e.g., r1)")
    parser.add_argument("filename", help="Filename for output detectorMap")
    args = parser.parse_args()
    detMap = makeDetectorMap(args.detector)
    detMap.writeFits(args.filename)
    return detMap


if __name__ == "__main__":
    main()
