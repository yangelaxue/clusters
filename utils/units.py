"""
Classes that store important physical constants in specified units.
"""

class CGS:
    """
    Physical constants in CGS units.
    """

    AU = 1.49597892e13          # Astronomical unit
    pc = 3.0856775807e18        # Parcec
    mp = 1.67262171e-24         # Proton mass
    me = 9.1093826e-28          # Electron mass
    c = 2.99792458e10           # Light speed
    yr = 365.25 * 24 * 60**2    # Year
    e = 4.80320425e-10          # Elementary proton charge
    eV = 1.602176463158e-12     # Electron Volt in erg.
    kB = 1.3806505e-16          # Boltzmann constant
    sigmaT = 6.65               # Thomson cross section for an electron
    G = 6.67430e-8              # Newton's constant
    
    Msun = 1.989e+33            # Solar mass


    CONST_eVtoK = 11606

    L_unit = 'cm'
    v_unit = 'cm/s'
    t_unit = 's'
    rho_unit = 'g/cm^3'
    force_unit = 'dyne'
    p_unit = 'dyne/cm^2'
    B_unit = 'G'