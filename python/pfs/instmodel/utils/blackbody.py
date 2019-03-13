import numpy as np

def blackbody(lam, temp):
    """Blackbody(Angstroms, Kelvin) -> nJy/sr"""

    speedOfLight = 3.0e8  # m/s
    planck = 6.63e-34  # J.s
    boltzman = 1.38e-23  # J/K
    wavelength = lam*1.0e-9  # m
    scale = 1.0e26*1.0e9   # W/m^2/Hz/sr --> nJy/sr
    flux = 2*planck*speedOfLight/wavelength**3/(np.exp(planck*speedOfLight/wavelength/boltzman/temp) - 1)
    return scale*flux
