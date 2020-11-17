from functools import partial
import numpy as np
from pfs.datamodel.pfsConfig import PfsDesign, FiberStatus
from pfs.instmodel.utils.generators import keepOdd, keepEven


def keepFibers(design, generator):
    """Keep a selection of fibers, blocking the remainder

    Parameters
    ----------
    design : `pfs.datamodel.PfsDesign`
        PfsDesign to modify.
    generator : callable
        Generator object that operates on an array of ``fiberId`` to provide
        the ``fiberId`` values that will be kept.
    """
    selection = set(generator(design.fiberId))
    design.fiberStatus = np.array([tt if ff in selection else FiberStatus.BLACKSPOT
                                   for ff, tt in zip(design.fiberId, design.fiberStatus)])


def shuffleFibers(design):
    """Shuffle the fibers

    Parameters
    ----------
    design : `pfs.datamodel.PfsDesign`
        PfsDesign to modify.
    """
    indices = np.arange(len(design))
    np.random.shuffle(indices)

    design.tract = design.tract[indices]
    design.patch = [design.patch[ii] for ii in indices]
    design.ra = design.ra[indices]
    design.dec = design.dec[indices]
    design.catId = design.catId[indices]
    design.objId = design.objId[indices]
    design.targetType = design.targetType[indices]
    design.fiberStatus = design.fiberStatus[indices]
    design.fiberFlux = [design.fiberFlux[ii] for ii in indices]
    design.filterNames = [design.filterNames[ii] for ii in indices]
    design.pfiNominal = design.pfiNominal[indices]


# Menu of available transmutations (mapping of name --> callable that will modify a PfsDesign)
TRANSMUTATIONS = {"odd": partial(keepFibers, generator=keepOdd),
                  "even": partial(keepFibers, generator=keepEven),
                  "shuffle": shuffleFibers
                  }


def transmutePfsDesign(inPfsDesignId, transmutation, outPfsDesignId, dirName="."):
    """Transumte a PfsDesign

    Parameters
    ----------
    inPfsDesignId : `int`
        Input pfsDesignId.
    transmutation : `str`
        Transmutation to apply. Must be one of:
        - ``odd``: keep odd fibers only.
        - ``even``: keep even fibers only.
        - ``shuffle``: shuffle fibers.
    outPfsDesignId : `int`
        Output pfsDesignId.
    dirName : `str`, optional
        Directory containing PfsDesign files.

    Returns
    -------
    design : `pfs.datamodel.PfsDesign`
        Transmuted PfsDesign.
    """
    design = PfsDesign.read(inPfsDesignId)
    if transmutation not in TRANSMUTATIONS:
        raise RuntimeError(f"Unrecognised transmutation: {transmutation} not in {TRANSMUTATIONS.keys()}")

    TRANSMUTATIONS[transmutation](design)
    design.pfsDesignId = outPfsDesignId
    design.write(dirName=dirName)
    return design


def main():
    """Main entrypoint when running as a script"""
    from argparse import ArgumentParser
    parser = ArgumentParser(description="Transmute a PfsDesign")
    parser.add_argument("input", type=int, help="Input pfsDesignId")
    parser.add_argument("transmutation", choices=TRANSMUTATIONS.keys(), help="Transmutation to apply")
    parser.add_argument("output", type=int, help="Output pfsDesignId")
    parser.add_argument("--dir", default=".", help="Directory with PfsDesign files")
    args = parser.parse_args()
    transmutePfsDesign(args.input, args.transmutation, args.output, dirName=args.dir)


if __name__ == "__main__":
    main()
