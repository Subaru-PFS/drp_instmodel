class Psf(object):

    def psfAt(self, fiber, wave):
        ims, ctrs =  self.psfsAt([fiber], [wave])
        return ims[0]
    
