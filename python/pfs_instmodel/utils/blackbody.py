import numpy as np

def blackbody(lam, temp):
    """ Crude blackbody(Angstroms, Kelvin) -> ergs/(s cm**2 nm) """

    c_nm = 2.99792458e17                 # nm/s
    hc = 1.988e-6

    fnu = hc/lam**3 * 1/(np.exp(1.43868e8/(lam*temp)) - 1.0)

    flam = fnu*c_nm/lam**2

    return flam.astype('f4')

