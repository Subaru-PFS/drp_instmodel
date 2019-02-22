import numpy as np
import pfs.instmodel.fieldDefinitions

"""Generate a PfiDesign"""

__all__ = ["run", "main"]


def run(field, pfiDesignId, seed=None, dirName="."):
    """Make a PfiDesign

    Parameters
    ----------
    field : `str`
        Field type; must be defined in ``pfs.instmodel.fieldDefinitions``.
    pfiDesignId : `int`
        Top-end design identifier.
    seed : `int`, optional
        RNG seed.
    dirName : `str`
        Output directory.

    Returns
    -------
    pfiDesign : `pfs.datamodel.PfiDesign`
        Top-end design.
    """
    if field not in pfs.instmodel.fieldDefinitions.__all__:
        raise RuntimeError("Unrecognised field: %s" % (field,))
    rng = np.random.RandomState(seed) if seed is not None else None
    pfiDesign = getattr(pfs.instmodel.fieldDefinitions, field)(rng=rng)
    pfiDesign.pfiDesignId = pfiDesignId
    pfiDesign.write(dirName)
    return pfiDesign


def main():
    """Main entrypoint when running as a script"""
    from argparse import ArgumentParser
    parser = ArgumentParser(description="Generate a PfiDesign")
    parser.add_argument("--field", choices=pfs.instmodel.fieldDefinitions.__all__, required=True,
                        help="Field type")
    parser.add_argument("--pfiDesignId", type=int, required=True, help="Top-end design identifier")
    parser.add_argument("--seed", type=int, help="RNG seed")
    parser.add_argument("--dirName", default=".", help="Output directory")
    args = parser.parse_args()
    run(args.field, args.pfiDesignId, args.seed, args.dirName)


if __name__ == "__main__":
    main()
