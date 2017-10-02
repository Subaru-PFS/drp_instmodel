#!/usr/bin/env python

import argparse
import logging
import os
import re

import pfs_tools.par
import pfs_instmodel.simImage as simImage
import pfs_instmodel.sky as pfsSky
import pfs_instmodel.spectrum as pfsSpectrum
reload(pfsSpectrum)

def makeSim(band, fieldName, fiberFilter=None,
            frd=23, focus=0, date='2016-10-26', psf=None, dtype='u2',
            addNoise=True, combSpacing=50, shiftPsfs=True,
            constantPsf=False, constantX=False,
            xOffset=0.0, yOffset=0.0,
            realBias=None,
            logger=None):
    """ Construct a simulated image. 

    Parameters
    ----------

    band : str
       The name of the wavelength band (for PFS: IR, Red, Blue)
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
        
    simID = dict(band=band, frd=frd, focus=focus, date=date)

    sim = simImage.SimImage(band, simID=simID, psf=psf, dtype=dtype,
                            addNoise=addNoise,
                            constantPsf=constantPsf, constantX=constantX,
                            slitOffset=(xOffset/1000.0, yOffset/1000.0),
                            logger=logger)
    skyModel = pfsSky.StaticSkyModel(band) # plus field info....
    flatSpectrum = pfsSpectrum.FlatSpectrum(sim.detector, gain=20.0)
    slopeSpectrum = pfsSpectrum.SlopeSpectrum(sim.detector, gain=20.0)
    combSpectrum = pfsSpectrum.CombSpectrum(spacing=combSpacing, 
                                            gain=200000.0)
    
    field = loadField(fieldName)

    logger.info("addNoise=%s" % (addNoise))

    fibers = []
    spectra = []
    for f in field:
        if fiberFilter and f.fiberId not in fiberFilter:
            continue
        if f.type == 'UNPLUGGED':
            continue
        if f.type == 'SKY':
            fibers.append(f.fiberId)
            spectra.append(skyModel.getSkyAt(ra=f.ra, dec=f.dec))
        elif f.type == 'SIMFLAT':
            fibers.append(f.fiberId)
            spectra.append(flatSpectrum)
        elif f.type == 'SIMCOMB':
            fibers.append(f.fiberId)
            spectra.append(combSpectrum)
        elif f.type == 'SIMSLOPE':
            fibers.append(f.fiberId)
            spectra.append(slopeSpectrum)
        elif f.type == 'SIMARC':
            fibers.append(f.fiberId)
            arcSpectrum = pfsSpectrum.ArcSpectrum(*f.args)
            spectra.append(arcSpectrum)
        elif f.type == 'OBJECT':
            raise RuntimeError("sorry, we don't do %s spectra yet" % f.type)

            # Per JEG, we expect to insert object spectra differently
            # from object spectra, so add them in one at a time.
            fibers.append(f.fiberId)
            spectra.append(skyModel.getSkyAt(ra=f.ra, dec=f.dec))
            
            fibers.append(f.fiberId)
            spectra.append(fetchSpectrumSomehow(f.object))   # Boom for now.
        else:
            raise RuntimeError("sorry, we don't do %s spectra yet" % f.type)

    sim.addFibers(fibers,
                  spectra=spectra,
                  shiftPsfs=shiftPsfs)
    return sim

def displayImage(img):
    import ds9
    disp = ds9.ds9()
    disp.set_np2arr(img)
    
def loadField(fieldName):
    """ Load the given field definition. 

    Currently just looks in a static file (examples/sampleField.py).  

    """
    
    # Decide on where to save field definitions, and add the usual path crap
    fields = pfs_tools.par.loadParFile(os.path.join(os.environ["DRP_INSTMODEL_DIR"], "examples", "sampleField.py"))
    return fields[fieldName]
    
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
            fullRange.extend(range(int(parts[0]),int(parts[1])+1))
        # print "after %s, r=%s" % (r, fullRange)
    return fullRange

def str2bool(v):
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

def main(args=None):
    """ called by __main__, or pass in a string as if it were a command line. 

    Get arg help from the command line.
    
    Returns
    -------

    sim : the SimInage object
    
    """
    if isinstance(args, basestring):
        import shlex
        args = shlex.split(args)

    helpDoc = \
""" 
Examples
--------

Generate an image file of two 3-fiber groups of sky spectra,
currently as defined in :download:`examples/sampleField/py <../../examples/sampleField.py>`
    
   --band=IR --output=irsim.fits --field=field1
"""

    # Configure the default formatter and logger.
    logging.basicConfig(datefmt = "%Y-%m-%d %H:%M:%S",
                        format = "%(asctime)s.%(msecs)03dZ %(name)-12s %(levelno)s %(filename)s:%(lineno)d %(message)s")
    logger = logging.getLogger()
    
    parser = argparse.ArgumentParser(description="generate a simulated image", 
                                     formatter_class=argparse.RawDescriptionHelpFormatter,
                                     epilog=helpDoc)
    parser.add_argument('-b', '--band', action='store', required=True)
    parser.add_argument('-F', '--field', action='store', required=True)
    parser.add_argument('-o', '--output', action='store', default=None)
    parser.add_argument('-f', '--fibers', action='store', default=None)
    parser.add_argument('-v', '--verbose', action='store_true')
    parser.add_argument('--focus', action='store', default=0, type=int)
    parser.add_argument('--frd', action='store', default=23, type=int)
    parser.add_argument('--date', action='store', default='2016-10-26')
    parser.add_argument('--dtype', action='store', default='u2')
    parser.add_argument('--xoffset', action='store', type=float, default=0.0,
                        help='shift in slit position along slit, in microns')
    parser.add_argument('--yoffset', action='store', type=float, default=0.0,
                        help='shift in slit position along dispersion, in microns')
    parser.add_argument('--noNoise', action='store_true')
    parser.add_argument('--allOutput', action='store_true',
                        help='whether to add (many) additional HDUs abut the simulation')
    parser.add_argument('--realBias', action='store', type=str2bool, default=True)
    parser.add_argument('--realFlat', action='store', type=str2bool, default=True,
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
    
    parser.add_argument('-d', '--ds9', action='store_true', default=False)

    args = parser.parse_args(args)
    logger.setLevel(logging.INFO if not args.verbose else logging.DEBUG)
    logger.debug('starting simImage logging')
    
    fibers = expandRangeArg(args.fibers)

    sim = makeSim(args.band, fieldName=args.field, 
                  fiberFilter=fibers,
                  frd=args.frd, focus=args.focus, date=args.date,
                  dtype=args.dtype,
                  addNoise=not args.noNoise,
                  combSpacing=args.combSpacing,
                  shiftPsfs=args.shiftPsfs,
                  constantPsf=args.constantPsf,
                  constantX=args.constantX,
                  xOffset=args.xoffset,
                  yOffset=args.yoffset,
                  realBias=args.realBias,
                  logger=logger)
    if args.output != 'no':
        sim.writeTo(args.output, addNoise=not args.noNoise,
                    compress=args.compress, allOutput=args.allOutput,
                    realBias=args.realBias, realFlat=args.realFlat,
                    imagetyp=args.imagetyp)
    if args.ds9:
        displayImage(sim.image)
    return sim

if __name__ == "__main__":
    main()
