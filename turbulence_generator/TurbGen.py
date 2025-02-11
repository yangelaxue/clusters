"""
TurbGen written in Python, based on C. Federrath's work:

--------------------------

    Please see Federrath et al. (2010, A&A 512, A81) for details and cite :)

    DESCRIPTION
    Contains functions to compute the time-dependent physical turbulent vector
    field used to drive or initialise turbulence in hydro codes such as AREPO,
    FLASH, GADGET, PHANTOM, PLUTO, QUOKKA, etc.
    The driving sequence follows an Ornstein-Uhlenbeck (OU) process.

    For example applications see Federrath et al. (2008, ApJ 688, L79);
    Federrath et al. (2010, A&A 512, A81); Federrath (2013, MNRAS 436, 1245);
    Federrath et al. (2021, Nature Astronomy 5, 365)

    AUTHOR: Christoph Federrath, 2008-2024

--------------------------

Angela Xue, 2025
"""

sol_weight = 1.
ncmp = 1.

def set_solenoidal_weight_normalisation():
    sol_weight_norm = (3.0/ncmp)**.5 * 3.0**.5 * 1.0/(1.0 - 2.0*sol_weight + ncmp*sol_weight**2)**.5

