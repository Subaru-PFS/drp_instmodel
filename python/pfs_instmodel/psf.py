class Psf(object):

    def psfAt(self, fiber, wave):
        """ Instantiate a single PSF at the given position.
        """
        
        ims, ctrs =  self.psfsAt([fiber], [wave])
        return ims[0]
    
