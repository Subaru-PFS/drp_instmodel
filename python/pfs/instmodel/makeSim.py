from contextlib import contextmanager

import argparse
import logging
import os
import re
from types import SimpleNamespace
import numpy as np

from pfs.datamodel.pfsConfig import TargetType
import pfs.instmodel.simImage as simImage
import pfs.instmodel.sky as pfsSky
from pfs.datamodel import PfsDesign
from pfs.instmodel.makePfsConfig import makePfsConfig
from .lightSource import LightSource, Lamps


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


def makeSim(detector, pfsDesignId=0, expId=0, fiberFilter=None,
            lamps=Lamps.NONE, spectraDir=None,
            frd=None, focus=0, date=None, psf=None, dtype='u2',
            everyNth=20,
            addNoise=True, domeOpen=True, combSpacing=50, shiftPsfs=True,
            constantPsf=False, constantX=False,
            xOffset=0.0, yOffset=0.0,
            realBias=None,
            dirName=".",
            logger=None):
    """ Construct a simulated image.

    Parameters
    ----------

    detector : str
       The name of the detector ("r1", "n3", "b4", etc)
    fiberFilter : list of integers, optional
       Only process the given fiber IDs.

    Returns
    -------
s
    sim : a SimImage object. Notable member is .image
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
    design = PfsDesign.read(pfsDesignId, dirName=dirName)
    config = makePfsConfig(design, expId)

    logger.info("addNoise=%s" % (addNoise))

    fibers = config.fiberId
    if domeOpen:
        doSkyForFiber = [tt in set([TargetType.SCIENCE, TargetType.SKY, TargetType.FLUXSTD]) for
                         tt in config.targetType]
    else:
        doSkyForFiber = np.zeros_like(fibers, dtype=bool)
    source = LightSource(domeOpen, lamps, skyModel, design, spectraDir)
    spectra = [source.getSpectrum(fiberId) for fiberId in fibers]
    sim.addFibers(fibers,
                  spectra=spectra,
                  doSkyForFiber=doSkyForFiber,
                  skySpectrum=skyModel.getSkyAt(),
                  shiftPsfs=shiftPsfs)
    return SimpleNamespace(image=sim, config=config)


def displayImage(img):
    import ds9
    disp = ds9.ds9()
    disp.set_np2arr(img)


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
    parser.add_argument('-e', '--expId', type=int, required=True, help="Exposure identifier")
    parser.add_argument('-p', '--pfsDesignId', type=int, required=True, help="pfsDesignId")
    parser.add_argument('--dirName', default=".", help="Directory in which to write")
    parser.add_argument('--spectraDir', default=".", help="Directory from which to read spectra")
    parser.add_argument('--pfsConfig', default=False, action="store_true", help="Generate pfsConfig?")
    parser.add_argument('--lamps', default="", help="List of lamps that are on (QUARTZ,NE,HG,XE,CD,KR)")
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
    parser.add_argument('--domeClosed', action='store_true', default=False, help="Close the dome (so no sky)")
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
    parser.add_argument('--detectorMap', help="Name for detectorMap file")

    parser.add_argument('--ds9', action='store_true', default=False)

    args = parser.parse_args(args)
    logger.setLevel(logging.INFO if not args.verbose else logging.DEBUG)
    logger.debug('starting simImage logging')

    if args.realBias is False:
        args.realBias = None

    fibers = expandRangeArg(args.fibers)
    lamps = Lamps.fromString(args.lamps)

    with pdbOnException(args.pdb):
        sim = makeSim(args.detector,
                      pfsDesignId=args.pfsDesignId,
                      expId=args.expId,
                      fiberFilter=fibers,
                      lamps=lamps,
                      spectraDir=args.spectraDir,
                      frd=args.frd, focus=args.focus, date=args.date,
                      dtype=args.dtype,
                      everyNth=args.everyNth,
                      addNoise=not args.noNoise,
                      domeOpen=not args.domeClosed,
                      combSpacing=args.combSpacing,
                      shiftPsfs=args.shiftPsfs,
                      constantPsf=args.constantPsf,
                      constantX=args.constantX,
                      xOffset=args.xoffset,
                      yOffset=args.yoffset,
                      realBias=args.realBias,
                      dirName=args.dirName,
                      logger=logger)

    site = 'F'  # Fake
    category = 'A'  # Science (as opposed to metrology camera, etc.)
    visit = args.expId
    spectrograph = int(args.detector[1])
    armNum = {'b': 1, 'r': 2, 'n': 3, 'm': 4}[args.detector[0]]
    imageName = "PF%1s%1s%06d%1d%1d.fits" % (site, category, visit, spectrograph, armNum)
    with pdbOnException(args.pdb):
        sim.image.writeTo(os.path.join(args.dirName, imageName), addNoise=not args.noNoise,
                          exptime=args.exptime, pfsDesignId=args.pfsDesignId,
                          compress=args.compress, allOutput=args.allOutput,
                          realBias=args.realBias, realFlat=args.realFlat,
                          imagetyp=args.imagetyp, addCards=lamps.toFitsCards())
        if args.pfsConfig:
            sim.config.write()
        if args.detectorMap:
            sim.image.psf.makeDetectorMap(args.detectorMap)

    if args.ds9:
        displayImage(sim.image)
    return sim


if __name__ == "__main__":
    if True:
        main()
    else:
        import cProfile
        cProfile.run("main()", "profile.pstats")