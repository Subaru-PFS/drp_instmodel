import numpy as np
from pfs.instmodel.slit import Slit, NUM_FIBERS
from pfs.instmodel.utils.generators import keepOdd, keepEven
from pfs.instmodel.makePfsConfig import makeScienceDesign

"""Generate a PfsDesign"""

SPECTROGRAPHS = [1, 2, 3, 4]  # Spectrograph numbers


def holeToFiberId(holes, spectrograph):
    """Convert hole numbers to fiberIds

    Parameters
    ----------
    holes : array_like of `int`
        Hole numbers (1..651).
    spectrograph : `int`
        Spectrograph number (1..4).

    Returns
    -------
    fiberId : `numpy.ndarray` of `int`
        Fiber identifiers.
    """
    return NUM_FIBERS*(spectrograph - 1) + np.array(holes)


def parseFibers(fibers, spectrograph):
    """Parse the list of fibers

    Parameters
    ----------
    fibers : `str`
        Whitespace-separated list of fiber IDs, or symbolic name representing a
        list of fiber IDs (single, double, lam, all, odd, even, fifteen).
    spectrograph : `int`
        Spectrograph identifier (1-4).

    Returns
    -------
    fiberIds : `list` of `int`
        Fiber identifiers.
    """
    allFibers = Slit(spectrograph).scienceFibers
    menu = {"single": holeToFiberId([315], spectrograph),
            "double": holeToFiberId([311, 315], spectrograph),
            "lam": holeToFiberId([2, 65, 191, 254, 315, 337, 400, 463, 589, 650], spectrograph),
            "all": allFibers,
            "odd": [ii for ii in keepOdd(allFibers)],
            "even": [ii for ii in keepEven(allFibers)],
            "fifteen": allFibers[np.linspace(0, len(allFibers), 15, False, dtype=int)],
    }

    if fibers.lower() in menu:
        return menu[fibers.lower()]
    fibers = (int(ff) for ff in fibers.split())
    return [ff for ff in fibers if ff in allFibers]


def parseObjId(text):
    """Parse the list of object IDs

    Parameters
    ----------
    text : `str`
        Whitespace-separated list of object IDs.

    Returns
    -------
    objIds : `list` of `int`
        Object identifiers.
    """
    return [int(ii) for ii in text.split()]


def main():
    """Main entrypoint when running as a script"""
    from argparse import ArgumentParser
    parser = ArgumentParser(description="Generate a PfsDesign")
    parser.add_argument("--fibers", required=True,
                        help="Fibers to light (single, double, lam, all, odd, even, fifteen; or "
                             "space-delimited integers)")
    parser.add_argument("--spectrograph", type=int, nargs="+", default=[1, 2, 3, 4],
                        help="Spectrographs to be illuminated")
    parser.add_argument("--pfsDesignId", type=int, required=True, help="Top-end design identifier")
    parser.add_argument("--fracSky", type=float, default=0.2, help="Fraction of fibers to use for sky")
    parser.add_argument("--fracFluxStd", type=float, default=0.1,
                        help="Fraction of fibers to use for flux standards")
    parser.add_argument("--minScienceMag", type=float, default=18.0,
                        help="Minimum AB magnitude of science objects")
    parser.add_argument("--maxScienceMag", type=float, default=24.0,
                        help="Maximum AB magnitude of science objects")
    parser.add_argument("--fluxStdMag", type=float, default=18.0, help="AB Magnitude of flux standards")
    parser.add_argument("--scienceCatId", type=int, default=0, help="Catalog ID for science targets")
    parser.add_argument("--scienceObjId", type=parseObjId,
                        help="Object IDs for science targets (space-delimited integers)")
    parser.add_argument("--seed", type=int, help="RNG seed")
    parser.add_argument("--dirName", default=".", help="Output directory")
    parser.add_argument("--arms", default="brn", help="Arms to be exposed")
    args = parser.parse_args()

    rng = np.random.RandomState(args.seed) if args.seed is not None else None
    spectrographs = args.spectrograph or SPECTROGRAPHS
    litFibers = np.concatenate([parseFibers(args.fibers, sp) for sp in spectrographs])
    allFibers = np.concatenate([Slit(ss).scienceFibers for ss in SPECTROGRAPHS])
    unlitFibers = np.array(sorted(set(allFibers) - set(litFibers)))

    pfsDesign = makeScienceDesign(args.pfsDesignId, litFibers, args.arms,
                                  args.fracSky, args.fracFluxStd,
                                  args.minScienceMag, args.maxScienceMag,
                                  args.fluxStdMag,
                                  args.scienceCatId, args.scienceObjId,
                                  rng=rng,
                                  unlitFiberIds=unlitFibers)
    pfsDesign.write(args.dirName)


if __name__ == "__main__":
    main()
