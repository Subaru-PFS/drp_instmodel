class SpecMap(object):
    """ Accumulate and provide the various fiberNum, wavelength, x, y mappings. 

    I wonder if we could just have a single mapper?
      1. start with tracing the fibers on a flat. Get fiberNum->[(pixelx,pixely)]
      2. add arc/wavelength solution. Add fiberNum->[(pixelx,wavelength)]
      
    """
    
    def __init__(self, img, ivar=None):
        """ Construct ourselves from the given flat image. """
        pass

    def addWavelengthSolution(self, waveSolution):
        pass
    
    def fibersToXy(self, fibers):
        """ Return the (x, y) vectors for the given fibers. """
        pass

    def fibersToWavelength(self, fibers):
        """ Return the (x, y, wavelength) vectors for the given fibers. """
        pass

    def fibersToXySplines(self, fibers):
        """ Return the (x,y) generating splines for the given fibers. """
        pass

    def fibersToWavelengthSplines(self, fibers):
        """ Return the (x,wavelength) generating splines for the given fibers. """
        pass

    
