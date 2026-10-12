"""
Classes that store important physical constants in specified units.

Author: Angela Xue
Date: Oct 2026
"""

from numpy import pi

class CGS:
    """ Physical constants in CGS units. """

    c = 2.99792458e10           # Light speed
    h = 6.626e-27               # Planck's constant
    m_Pl = 2.176e-5             # Planck mass
    l_Pl = 1.616e-33            # Planck length

    # Subatomic units.
    mp = 1.67262171e-24         # Proton mass
    mn = 1.6749e-24             # Neutron mass
    me = 9.1093826e-28          # Electron mass
    amu = 1.66054e-24            # Atomic mass unit
    a0 = 5.29177210544e-9       # Bohr radius
    sigmaT = 6.65               # Thomson cross section for an electron
    e = 4.80320425e-10          # Elementary proton charge
    eV = 1.602176463158e-12     # Electron Volt in erg.

    # Electromagnetic constants.
    epsilon0 = 1./(4*pi)        # permitivity = 1; % of vaccum space in cgs Taken from KS script
    mu0 = 4*pi/c**2             # permeability= 1; % of vaccum space in cgs Taken from KS script

    # Statistical units.
    kB = 1.3806505e-16          # Boltzmann constant
    R = 8.31446261815324        # Gas constant
    c2 = h*c/kB                 # Radiation constant
    avog = 6.02214076e23        # Avogadro's number

    # Astrophysical units.
    AU = 1.49597892e13          # Astronomical unit
    pc = 3.0856775807e18        # Parcec
    yr = 365.25 * 24 * 60**2    # Year
    G = 6.67430e-8              # Newton's constant
    Msun = 1.989e+33            # Solar mass