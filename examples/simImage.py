
import argparse
import numpy
import os
import re

import pfs_tools.par
import pfs_instmodel.simImage as simImage
import pfs_instmodel.sky as pfsSky

def displayImage(img):
    import ds9
    disp = ds9.ds9()
    disp.set_np2arr(img)
    
def loadField(fieldName):
    """ load the given field definition. Currently just looks in a static file (examples/sampleField.py).  """
    
    # Decide on where to save field definitions, and add the usual path crap
    fields = pfs_tools.par.loadParFile(os.path.join(os.environ["PFS_INSTMODEL_DIR"], "examples", "sampleField.py"))
    return fields[fieldName]
    
def makeSim(band, fieldName=None, fibers=None, everyNthPsf=50):
    sim = simImage.SimImage(band)
    sky = pfsSky.StaticSkyModel(band) # plus field info....

    if fieldName:
        field = loadField(fieldName)

        fibers = []
        spectra = []
        for f in field:
            if f.type == 'UNPLUGGED':
                continue
            if f.type == 'SKY':
                fibers.append(f.fiberId)
                spectra.append(sky)
            else:
                raise RuntimeError("sorry, we don't do %s spectra yet" % f.type)
    else:
        if not fibers:
            # Interesting: At this point I don't know how many fibers there are.
            fibers = numpy.concatenate([numpy.arange(3),
                                        numpy.arange(3) + 100,
                                        300 - numpy.arange(3)])
        spectra = [sky]*len(fibers)
        
    sim.addFibers(fibers,
                  spectra=spectra,
                  everyNthPsf=everyNthPsf)
    return sim

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
    if isinstance(args, basestring):
        args = args.split()

    helpDoc = \
"""
Examples
--------

Generate an image of two 3-fiber groups of sky specta:
   -b IR -o irsim.fits -f 0-2,290-292
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
