from contextlib import contextmanager

import argparse
import logging
import os
import re
from types import SimpleNamespace
import numpy as np
from datetime import datetime, timedelta
import yaml

from pfs.datamodel.pfsConfig import TargetType, FiberStatus
from pfs.instmodel.arm import Arm
import pfs.instmodel.simImage as simImage
import pfs.instmodel.sky as pfsSky
from pfs.datamodel import PfsDesign
from pfs.instmodel.makePfsConfig import makePfsConfig
from .lightSource import LightSource, Lamps
from .slit import Slit
from .detector import Detector


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


def makeSim(detector, pfsDesignId=0, fiberFilter=None,
            lamps=Lamps.NONE, objSpectraDir=None,
            catConfig=None,
            frd=None, focus=0, date=None, psf=None, dtype='u2',
            everyNth=20,
            domeOpen=True, combSpacing=50, shiftPsfs=True,
            constantPsf=False, constantX=False,
            zenithDistance=45.0, aerosol=1.0, pwv=1.6, extinctSky=False,
            xOffset=0.0, yOffset=0.0,
            realBias=None,
            dirName=".",
            logger=None):
    """ Construct a simulated image.

    Parameters
    ----------

    detector : str
       The name of the detector ("r1", "n3", "b4", etc)
    pfsDesignId: int
        The ID of the corresponding pfsDesign file.
        The location of that file will be read from the ``dirName`` argument.
    fiberFilter : list of integers, optional
       Only process the given fiber IDs.

    Returns
    -------
    sim : a SimImage object. Notable member is .image

    Raise:
    -------
    yaml.YAMLError: if the ``catConfig`` file cannot be loaded.
    """

    if logger is None:
        logger = logging.getLogger()

    logger.info("args everyNth: %s", everyNth)

    simID = dict(detector=detector, frd=frd, focus=focus, date=date)

    arm = Arm.fromDetectorName(detector)
    skyModel = pfsSky.StaticSkyModel(arm, zenithDistance, aerosol, pwv, extinctSky)
    sim = simImage.SimImage(detector, skyModel, simID=simID, psf=psf, dtype=dtype,
                            everyNth=everyNth,
                            constantPsf=constantPsf, constantX=constantX,
                            slitOffset=(xOffset/1000.0, yOffset/1000.0),
                            logger=logger)
    design = PfsDesign.read(pfsDesignId, dirName=dirName)

    specFibers = set(Slit(Detector(detector).spectrograph).scienceFibers)
    fibers = np.array(sorted(set(design.fiberId).intersection(specFibers)))
    if domeOpen:
        skyTargetType = set([TargetType.SCIENCE, TargetType.SKY, TargetType.FLUXSTD])
        skyFiberStatus = set([FiberStatus.GOOD])
        doSkyForFiber = [design.targetType[ii] in skyTargetType and design.fiberStatus[ii] in skyFiberStatus
                         for ii in range(len(design))]
    else:
        doSkyForFiber = np.zeros_like(fibers, dtype=bool)

    if objSpectraDir:
        catConfigPath = os.path.join(objSpectraDir, catConfig)

        # Read in catalog configuration file.
        # This is a YAML file that contains a sequence of one or more mapping nodes,
        # where each mapping node consists of three key-value pairs:
        # - the `catId`,
        # - the catalog `name`, and
        # - the `rel_location` (relative location).
        with open(catConfigPath, 'r') as stream:
            catConfigData = yaml.safe_load(stream)

        # Convert yaml data to a dict indexed by catId,
        # with values corresponding to the relative location
        # of the catalog data with respect to a root directory
        # Note: could have designed YAML file such that a load would
        # return the required dict,
        # but that would make the YAML file less readable.
        catDict = dict()
        for catalogInfo in catConfigData:
            catDict[catalogInfo['catId']] = catalogInfo['rel_location']
    else:
        catDict = None

    source = LightSource(domeOpen, lamps, skyModel, design, objSpectraDir,
                         catDict)
    spectra = [source.getSpectrum(fiberId) for fiberId in fibers]
    sim.addFibers(fibers,
                  spectra=spectra,
                  doSkyForFiber=doSkyForFiber,
                  shiftPsfs=shiftPsfs)
    return SimpleNamespace(image=sim, design=design)


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
    parser.add_argument('-d', '--detector', action='store', required=True, help="Detector name, e.g., r1")
    parser.add_argument('--visit', type=int, action="append", required=True, help="Visit number")
    parser.add_argument('-p', '--pfsDesignId', type=int, required=True, help="pfsDesignId")
    parser.add_argument('--dirName', default=".", help="Directory in which to write")
    parser.add_argument('--objSpectraDir', nargs='?',
                        help="Directory from which to read object spectra")
    parser.add_argument('--catConfig', default='catalog_config.yaml',
                        help=("Name of catalog config file for object spectra."
                              "Location of file defined by --objSpectraDir"))
    parser.add_argument('--pfsConfig', default=False, action="store_true", help="Generate pfsConfig?")
    parser.add_argument('--type', required=True, choices=("bias", "dark", "flat", "arc", "object"),
                        help="Type of image")
    parser.add_argument('--imagetyp', help="Value for IMAGETYP header")
    parser.add_argument('--lamps', default="", help="List of lamps that are on (QUARTZ,NE,HG,XE,CD,KR)")
    parser.add_argument('-f', '--fibers', action='store', default=None)
    parser.add_argument('-v', '--verbose', action='store_true')
    parser.add_argument('--exptime', action='store', default=0.0, type=float)
    parser.add_argument('--focus', action='store', default=0, type=int)
    parser.add_argument('--frd', action='store', default=23, type=int)
    parser.add_argument('--date', default="2020-01-01T00:00:00", help="Base datetime to use for DATE-OBS")
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
    parser.add_argument('--combSpacing', action='store', type=float, default=50)
    parser.add_argument('--constantPsf', action='store', type=float, default=0,
                        help='Use a single PSF for the entire field.')
    parser.add_argument('--constantX', action='store_true',
                        help='Use the middle X-coordinate for all of each fiber.')
    parser.add_argument('--compress', action='store', default=None,
                        help='fitsio FITS compression type. e.g. RICE')
    parser.add_argument('--pdb', default=False, action='store_true', help="Launch pdb on exception?")
    parser.add_argument('--zenithDistance', type=float, default=45.0, help="Zenith distance (degrees)")
    parser.add_argument('--aerosol', type=float, default=1.0, help="Aerosol power-law index")
    parser.add_argument('--pwv', type=float, default=1.6, help="Precipitable water vapour (mm)")
    parser.add_argument('--extinctSky', default=False, action="store_true", help="Apply extinction to sky?")

    parser.add_argument('--ds9', action='store_true', default=False)

    args = parser.parse_args(args)
    logger.setLevel(logging.INFO if not args.verbose else logging.DEBUG)
    logger.debug('starting simImage logging')

    if args.realBias is False:
        args.realBias = None

    fibers = expandRangeArg(args.fibers)
    lamps = Lamps.fromString(args.lamps)

    lamps = Lamps.fromString(args.lamps)
    if args.type == "arc" and lamps in (Lamps.NONE, Lamps.QUARTZ):
        raise RuntimeError("Arc requested, but no lamps specified")
    if args.type == "flat":
        if lamps == Lamps.NONE:
            lamps = Lamps.QUARTZ
        if lamps != Lamps.QUARTZ:
            raise RuntimeError("Flat requested, but lamps specified (%s)" % (args.lamps,))
    if args.type == "bias" and args.exptime != 0:
        raise RuntimeError("Bias requested, but non-zero exposure time specified (%f)" % (args.exptime))
    if args.type == "dark" and lamps != Lamps.NONE:
        raise RuntimeError("Dark requested, but lamps specified (%s)" % (args.lamps,))
    if args.type == "object" and not args.objSpectraDir:
        raise RuntimeError("Object type requested, "
                           "but objSpectraDir is not specified.")

    with pdbOnException(args.pdb):
        sim = makeSim(args.detector,
                      pfsDesignId=args.pfsDesignId,
                      fiberFilter=fibers,
                      lamps=lamps,
                      objSpectraDir=args.objSpectraDir,
                      frd=args.frd, focus=args.focus,
                      dtype=args.dtype,
                      everyNth=args.everyNth,
                      domeOpen=args.type == "object",
                      combSpacing=args.combSpacing,
                      shiftPsfs=args.shiftPsfs,
                      constantPsf=args.constantPsf,
                      constantX=args.constantX,
                      zenithDistance=args.zenithDistance,
                      aerosol=args.aerosol,
                      pwv=args.pwv,
                      extinctSky=args.extinctSky,
                      xOffset=args.xoffset,
                      yOffset=args.yoffset,
                      realBias=args.realBias,
                      dirName=args.dirName,
                      catConfig=args.catConfig,
                      logger=logger)

    site = 'F'  # Fake
    category = 'A'  # Science (as opposed to metrology camera, etc.)
    spectrograph = int(args.detector[1])
    armNum = {'b': 1, 'r': 2, 'n': 3, 'm': 4}[args.detector[0]]
    date = datetime.fromisoformat(args.date)

    for visit in args.visit:
        imageName = "PF%1s%1s%06d%1d%1d.fits" % (site, category, visit, spectrograph, armNum)
        with pdbOnException(args.pdb):
            image = sim.image.clone()
            config = makePfsConfig(sim.design, visit)
            timestamp = (date + timedelta(minutes=visit)).isoformat()

            header = lamps.toFitsCards()
            header.append(("W_VISIT", visit, "Visit identifier"))
            header.append(("IMAGETYP", (args.imagetyp or args.type).upper(), "Image type"))
            image.writeTo(os.path.join(args.dirName, imageName), addNoise=not args.noNoise,
                          exptime=args.exptime, pfsDesignId=args.pfsDesignId, timestamp=timestamp,
                          compress=args.compress, allOutput=args.allOutput,
                          realBias=args.realBias, realFlat=args.realFlat,
                          imagetyp=args.type.upper(), addCards=header)
            if args.pfsConfig:
                config.write()

    if args.ds9:
        displayImage(sim.image)
    return sim


if __name__ == "__main__":
    if True:
        main()
    else:
        import cProfile
        cProfile.run("main()", "profile.pstats")
