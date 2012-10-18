
import argparse
import numpy
import os
import re

import pfs_tools.par
import pfs_instmodel.simImage as simImage
import pfs_instmodel.sky as pfsSky
import pfs_instmodel.spectrum as pfsSpectrum

def makeSim(band, fieldName=None, fibers=None, everyNthPsf=50):
    """ Construct a simulated image. 

    Parameters
    ----------

    band : str
       The name of the wavelength band (for PFS: IR, Red, Blue)
    fieldName : str, optional
       The name of the field definition with the fiber locations and targeting.
    fibers : list of integers, optional
       [REMOVE now that we have field defs. ]
    everyNthPsf : int, optional
       How many (sub-)pixels we can use the same PSF on.

    Returns
    -------

    sim : a SimImage object. Notable member is .image

    Notes
    -----

    We don't know how to generate anything other than sky spectra yet.

    The fieldName is currently just an entry in the fixed file 
    :download:`examples/sampleField.py <../../examples/sampleField.py>`
    """
    
    sim = simImage.SimImage(band)
    skyModel = pfsSky.StaticSkyModel(band) # plus field info....
    flatSpectrum = pfsSpectrum.FlatSpectrum(sim.detector)
    combSpectrum = pfsSpectrum.CombSpectrum(spacing=50)
    
    if fieldName:
        field = loadField(fieldName)

        fibers = []
        spectra = []
        for f in field:
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
            else:
                raise RuntimeError("sorry, we don't do %s spectra yet" % f.type)
    else:
        if not fibers:
            # Interesting: At this point I don't know how many fibers there are.
            fibers = numpy.concatenate([numpy.arange(3),
                                        numpy.arange(3) + 100,
                                        300 - numpy.arange(3)])
        spectra = [skyModel.getSkyAt() for f in fibers]
        
    sim.addFibers(fibers,
                  spectra=spectra,
                  everyNthPsf=everyNthPsf)
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
    fields = pfs_tools.par.loadParFile(os.path.join(os.environ["PFS_INSTMODEL_DIR"], "examples", "sampleField.py"))
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

def saveSim(sim, outputFile):
    import pyfits

    pyfits.writeto(outputFile, sim.image, checksum=True, clobber=True)
    
def main(args=None):
    """ called by __main__, or pass in a string as if it were a command line. 

    Get arg help from the command line.
    
    Returns
    -------

    sim : the SimInage object
    
    """
    if isinstance(args, basestring):
        args = args.split()

    helpDoc = \
""" 
Examples
--------

Generate an image file of two 3-fiber groups of sky spectra,
currently as defined in :download:`examples/sampleField/py <../../examples/sampleField.py>`
    
   --band=IR --output=irsim.fits --field=field1
"""
    
    parser = argparse.ArgumentParser(description="generate a simulated image", 
                                     formatter_class=argparse.RawDescriptionHelpFormatter,
                                     epilog=helpDoc)
    parser.add_argument('-b', '--band', action='store', required=True)
    parser.add_argument('-o', '--output', action='store', default=None)
    parser.add_argument('-f', '--fibers', action='store', default=None)
    parser.add_argument('-F', '--field', action='store', default=None)
    parser.add_argument('--everyNth', action='store', type=int, default=50)
    parser.add_argument('-d', '--ds9', action='store_true', default=False)

    args = parser.parse_args(args)

    fibers = expandRangeArg(args.fibers)

    sim = makeSim(args.band, fieldName=args.field, 
                  fibers=fibers, everyNthPsf=args.everyNth)
    if args.output:
        saveSim(sim, args.output)
    if args.ds9:
        displayImage(sim.image)
    return sim

if __name__ == "__main__":
    main()
