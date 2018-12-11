#!/usr/bin/env python

from importlib import reload
from contextlib import contextmanager

import argparse
import logging
import os
import re
from types import SimpleNamespace

from .utils import schema
import pfs_instmodel.simImage as simImage
import pfs_instmodel.sky as pfsSky
from .spectrumLibrary import SpectrumLibrary
import pfs_instmodel.spectrum as pfsSpectrum
reload(pfsSpectrum)


@contextmanager
def pdbOnException(enable=True):
    """Context manager to drop into pdb if there's an exception"""
    try:
        yield
    except Exception:
        if enable:
            import traceback
            traceback.print_exc()
            import pdb
            pdb.post_mortem()
        raise


def makeSim(detector, fieldName, pfiDesignId=0, expId=0, fiberFilter=None,
            frd=None, focus=0, date=None, psf=None, dtype='u2',
            everyNth=20,
            addNoise=True, combSpacing=50, shiftPsfs=True,
            constantPsf=False, constantX=False,
            xOffset=0.0, yOffset=0.0,
            realBias=None,
            logger=None):
    """ Construct a simulated image.

    Parameters
    ----------

    detector : str
       The name of the detector ("r1", "n3", "b4", etc)
    fieldName : str
       The name of the field definition with the fiber locations and targeting.
    fiberFilter : list of integers, optional
       Only process the given fiber IDs.

    Returns
    -------

    sim : a SimImage object. Notable member is .image

    Notes
    -----

    We don't know how to generate anything other than sky spectra yet.

    The fieldName is currently just an entry in the fixed file
    :download:`examples/sampleField.py <../../examples/sampleField.py>`
    """

    if logger is None:
        logger = logging.getLogger()

    logger.info("args everyNth: %s", everyNth)

    simID = dict(detector=detector, frd=frd, focus=focus, date=date)

    sim = simImage.SimImage(detector, simID=simID, psf=psf, dtype=dtype,
                            everyNth=everyNth,
                            addNoise=addNoise,
                            constantPsf=constantPsf, constantX=constantX,
                            slitOffset=(xOffset/1000.0, yOffset/1000.0),
                            logger=logger)
    skyModel = pfsSky.StaticSkyModel(sim.detector.armName)  # plus field info....
    config = loadConfig(fieldName)
    config.pfiDesignId = pfiDesignId
    config.expId = expId

    logger.info("addNoise=%s" % (addNoise))

    fibers = config.fiberId
    library = SpectrumLibrary(skyModel)
    spectra = [library.getSpectrum(catId, objId) for catId, objId in zip(config.catId, config.objId)]
    sim.addFibers(fibers,
                  spectra=spectra,
                  shiftPsfs=shiftPsfs)
    return SimpleNamespace(image=sim, config=config)


def displayImage(img):
    import ds9
    disp = ds9.ds9()
    disp.set_np2arr(img)


def loadConfig(fieldName):
    """ Load the given field definition.

    Currently just looks in a static file (examples/sampleField.py).
    """
    # Decide on where to save field definitions, and add the usual path crap
    configs = schema.loadParFile(os.path.join(os.environ["DRP_INSTMODEL_DIR"],
                                              "examples", "sampleField.py"))
    return configs[fieldName]


def expandRangeArg(arg):
    """ Generate an array from a range expression. No error checking.

    Expression: R,R,R
    R: INT or INT-INT

    So 0,1,100-103,300 yields [0,1,100,101,102,103,300]
    """

    fullRange = []
    if not arg:
        return fullRange

    rangeList = re.split('[, ]', arg)
    for r in rangeList:
        parts = r.split('-')
        if len(parts) == 1:
            fullRange.append(int(parts[0]))
        else:
            fullRange.extend(range(int(parts[0]), int(parts[1])+1))
        # print "after %s, r=%s" % (r, fullRange)
    return fullRange


def str2bool(v):
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')


def str2file(v):
    if v.lower() in ('yes', 'true', 't'):
        return True
    elif v.lower() in ('no', 'false', 'f'):
        return False

    try:
        return int(v, base=10)
    except ValueError:
        return v


def main(args=None):
    """ called by __main__, or pass in a string as if it were a command line.

    Get arg help from the command line.

    Returns
    -------
    sim : the SimImage object
    """
    if isinstance(args, str):
        import shlex
        args = shlex.split(args)
    elif args is None:
        import sys
        args = sys.argv[1:]

    helpDoc = """
Examples
--------

Generate an image file of two 3-fiber groups of sky spectra,
currently as defined in :download:`examples/sampleField/py <../../examples/sampleField.py>`

   --detector=r1 --output=sim.fits --field=field1
"""
    # Configure the default formatter and logger.
    logging.basicConfig(
        datefmt="%Y-%m-%d %H:%M:%S",
        format="%(asctime)s.%(msecs)03dZ %(name)-12s %(levelno)s %(filename)s:%(lineno)d %(message)s"
    )
    logger = logging.getLogger()

    parser = argparse.ArgumentParser(description="generate a simulated image", 
                                     formatter_class=argparse.RawDescriptionHelpFormatter,
                                     epilog=helpDoc)
    parser.add_argument('-d', '--detector', action='store', required=True)
    parser.add_argument('-F', '--field', action='store', required=True)
    parser.add_argument('-e', '--expId', type=int, required=True, help="Exposure identifier")
    parser.add_argument('-p', '--pfiDesignId', type=int, default=0, help="pfiDesignId")
    parser.add_argument('-f', '--fibers', action='store', default=None)
    parser.add_argument('-v', '--verbose', action='store_true')
    parser.add_argument('--exptime', action='store', default=1, type=float)
    parser.add_argument('--focus', action='store', default=0, type=int)
    parser.add_argument('--frd', action='store', default=23, type=int)
    parser.add_argument('--date', action='store')
    parser.add_argument('--dtype', action='store', default='u2')
    parser.add_argument('--everyNth', action='store', default=20, type=int,
                        help='how many (oversampled) pixels to linearly interpolate over.')
    parser.add_argument('--xoffset', action='store', type=float, default=0.0,
                        help='shift in slit position along slit, in microns')
    parser.add_argument('--yoffset', action='store', type=float, default=0.0,
                        help='shift in slit position along dispersion, in microns')
    parser.add_argument('--noNoise', action='store_true')
    parser.add_argument('--allOutput', action='store_true',
                        help='whether to add (many) additional HDUs abut the simulation')
    parser.add_argument('--realBias', action='store', type=str2file, default='True')
    parser.add_argument('--realFlat', action='store', type=str2bool, default='False',
                        help='Apply an imaging flat. Use False/None to avoid.')
    parser.add_argument('--shiftPsfs', action='store_false')
    parser.add_argument('--imagetyp', action='store', default=None,
                        help='IMAGETYP,EXPTIME pair')
    parser.add_argument('--combSpacing', action='store', type=float, default=50)
    parser.add_argument('--constantPsf', action='store', type=float, default=0,
                        help='Use a single PSF for the entire field.')
    parser.add_argument('--constantX', action='store_true',
                        help='Use the middle X-coordinate for all of each fiber.')
    parser.add_argument('--compress', action='store', default=None,
                        help='fitsio FITS compression type. e.g. RICE')
    parser.add_argument('--pdb', default=False, action='store_true', help="Launch pdb on exception?")

    parser.add_argument('--ds9', action='store_true', default=False)

    args = parser.parse_args(args)
    logger.setLevel(logging.INFO if not args.verbose else logging.DEBUG)
    logger.debug('starting simImage logging')

    if args.realBias is False:
        args.realBias = None

    fibers = expandRangeArg(args.fibers)

    with pdbOnException(args.pdb):
        sim = makeSim(args.detector, fieldName=args.field,
                      pfiDesignId=args.pfiDesignId,
                      expId=args.expId,
                      fiberFilter=fibers,
                      frd=args.frd, focus=args.focus, date=args.date,
                      dtype=args.dtype,
                      everyNth=args.everyNth,
                      addNoise=not args.noNoise,
                      combSpacing=args.combSpacing,
                      shiftPsfs=args.shiftPsfs,
                      constantPsf=args.constantPsf,
                      constantX=args.constantX,
                      xOffset=args.xoffset,
                      yOffset=args.yoffset,
                      realBias=args.realBias,
                      logger=logger)

    site = 'F'  # Fake
    category = 'A'  # Science (as opposed to metrology camera, etc.)
    visit = args.expId
    spectrograph = int(args.detector[1])
    armNum = {'b': 1, 'r': 2, 'n': 3, 'm': 4}[args.detector[0]]
    imageName = "PF%1s%1s%06d%1d%1d.fits" % (site, category, visit, spectrograph, armNum)
    with pdbOnException(args.pdb):
        sim.image.writeTo(imageName, addNoise=not args.noNoise,
                          exptime=args.exptime, pfiDesignId=args.pfiDesignId,
                          compress=args.compress, allOutput=args.allOutput,
                          realBias=args.realBias, realFlat=args.realFlat,
                          imagetyp=args.imagetyp)
        sim.config.write()

    if args.ds9:
        displayImage(sim.image)
    return sim


if __name__ == "__main__":
    main()
