"""Approximate ice-fluid water equilibrium functions.

This module provides approximation functions (low-order polynomial and
power-law fits) for properties of ice and fluid water, e.g. the melting
and sublimation curves of ice.

:Functions:

* :func:`liqpressure`: Approximate the pressure at which ice at the
  given temperature will melt.
* :func:`liqtemperature`: Approximate the temperature at which ice at
  the given pressure will melt.
* :func:`vappressure`: Approximate the vapour pressure over ice at the
  given temperature.

"""

__all__ = ['liqpressure','liqtemperature','vappressure']

import numpy
import constants0

_TTP = constants0.TTP
_PTPE = constants0.PTPE
_C_PMELT = ((1195393.37,3), (80818.3159,25.75), (3338.2686,103.75))
_C_TMELT = (-1.66356104484551e-7,-2.13519241979406e-13,
    3.52967405341877e-20,-2.73184525236281e-26)
_C_PSUBL = ((-21.2144006,0.00333333333), (27.3203819,1.20666667),
    (-6.1059813,1.70333333))


## Functions
def liqpressure(temp):
    """Approximate ice melting pressure.
    
    Approximate the pressure at which ice at the given temperature will
    melt.
    
    :arg float temp: Temperature in K.
    :returns: Melting pressure in Pa.
    """
    tau = temp/_TTP
    pres = 1.
    for (a,b) in _C_PMELT:
        pres += a * (1 - tau**b)
    pres *= _PTPE
    return pres

def liqtemperature(pres):
    """Approximate ice melting temperature.
    
    Approximate the temperature at which ice at the given pressure will
    melt.
    
    :arg float pres: Pressure in Pa.
    :arg float temp: Melting temperature in K.
    """
    psi = pres/_PTPE - 1
    temp = 0.
    for coeff in _C_TMELT[::-1]:
        temp = temp*psi + coeff
    temp = (1 + temp*psi)*_TTP
    return temp

def vappressure(temp):
    """Approximate vapour pressure over ice.
    
    Approximate the vapour pressure over ice at the given temperature.
    
    :arg float temp: Temperature in K.
    :returns: Vapour pressure in Pa.
    """
    tau = temp/_TTP
    earg = 0.
    for (a,b) in _C_PSUBL:
        earg += a * tau**(b-1)
    pres = _PTPE * numpy.exp(earg)
    return pres

