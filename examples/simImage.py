

import argparse
import numpy
import re
import sys

import pfs_instmodel.simImage as simImage
import pfs_instmodel.sky as pfsSky

def displayImage(img):
    import ds9
    disp = ds9.ds9()
    disp.set_np2arr(img)
    
def makeSim(band, fibers=None, everyNthPsf=1):
    sim = simImage.SimImage(band)
    sky = pfsSky.StaticSkyModel(band)

    if not fibers:
        # Interesting: At this point I don't know how many fibers there are.
        fibers = numpy.concatenate([numpy.arange(3),
                                    numpy.arange(3) + 100,
                                    300 - numpy.arange(3)])
    spectra = [sky]*len(fibers)
        
    simage = sim.addFibers(fibers,
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

def getField(fieldName):
    pass

def main(args):
    if isinstance(args, basestring):
        args = args.split()
        
    parser = argparse.ArgumentParser(description="generate a simulated image")
    parser.add_argument('-b', '--band', action='store', required=True)
    parser.add_argument('-o', '--output', action='store', default=None)
    parser.add_argument('-f', '--fibers', action='store', default=None)
    parser.add_argument('-F', '--field', action='store', default=None)
    parser.add_argument('--everyNth', action='store', type=int, default=50)
    parser.add_argument('-d', '--ds9', action='store_true', default=False)

    res = parser.parse_args(args)

    fibers = expandRangeArg(res.fibers)

    sim = makeSim(res.band, fibers=fibers, everyNthPsf=res.everyNth)
    if res.ds9:
        displayImage(sim.image)
    return sim

if __name__ == "__main__":
    main(sys.argv)
    


