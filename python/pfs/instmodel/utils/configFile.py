
from past.builtins import execfile
def readfile(configFile, outDict=None):
    """ execute a python file, and return the final environment. """
    
    gdict = {}
    ldict = {}
        
    execfile(configFile, gdict) # , ldict)
    if outDict:
        outDict.update(gdict)
    return gdict
