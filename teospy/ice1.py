"""Ice Gibbs free energy.

This module implements the Gibbs free energy of ice and its derivatives with respect to temperature and pressure. This function applies to hexagonal ice 1 (denoted Ih), the relevant phase for most geophysical applications.

:Functions:

* ice_g: Ice Gibbs free energy with derivatives.

"""

__all__ = ['ice_g']

import cmath
import constants0

# Single constants
_TTP = constants0.TTP
_PTP = constants0.PTP
_PATM = constants0.PATM
_PI0 = _PATM/_PTP

# Constants used in empirical equations
_GCOEFFS = (
    (-6.32020233335886e5, 0.655022213658955, -1.89369929326131e-8,
        3.39746123271053e-15, -5.56464869058991e-22),
    complex(44.7050716285388,65.6876847463481),
    (complex(-72.597457432922,-78.100842711287),
        complex(-5.57107698030123e-5,4.64578634580806e-5),
        complex(2.34801409215913e-11,-2.85651142904972e-11)),
    (complex(3.68017112855051e-2,5.10878114959572e-2),
        complex(0.337315741065416,0.335449415919309)),
    -3327.33756492168
)


## Auxiliary functions
def _ice_g0(temp,pres):
    """Calculate ice Gibbs free energy.
    
    Calculate the specific Gibbs free energy of ice.
    
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :returns: Gibbs free energy in J/kg.
    """
    # Reduced variables
    tn = temp/_TTP
    pn = pres/_PTP
    _PI0 = _PATM/_PTP
    g = 0.
    
    # Power series and entropy components
    G0, TK, RK1, RK2, S0 = _GCOEFFS
    for (l,n) in enumerate(_GCOEFFS[0]):
        g += n * (pn-_PI0)**l
    g -= _GCOEFFS[4] * _TTP * tn
    
    # Residual terms including complex numbers
    sr = [_GCOEFFS[1], complex(0.0,0.0)]
    for (k,rk) in enumerate(_GCOEFFS[2]):
        sr[1] += rk * (pn-_PI0)**k
    for (tk,s) in zip(_GCOEFFS[3],sr):
        term = (tk-tn)*cmath.log(tk-tn) + (tk+tn)*cmath.log(tk+tn)
        term -= 2*tk*cmath.log(tk) + tn**2 / tk
        g += _TTP * (s*term).real
    return g

def _ice_dgdt(temp,pres):
    """Calculate ice Gibbs free energy T-derivative.
    
    Calculate the derivative of the specific Gibbs free energy of ice
    with respect to temperature.
    
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :returns: Gibbs free energy derivative in J/kg/K.
    """
    # Reduced variables
    tn = temp/_TTP
    pn = pres/_PTP
    _PI0 = _PATM/_PTP
    g_t = 0.
    
    # Power series and entropy components
    g_t += -_GCOEFFS[4]
    
    # Residual terms including complex numbers
    sr = [_GCOEFFS[1], complex(0.0,0.0)]
    for (k,rk) in enumerate(_GCOEFFS[2]):
        sr[1] += rk * (pn-_PI0)**k
    for (tk,s) in zip(_GCOEFFS[3],sr):
        term = -cmath.log(tk-tn) + cmath.log(tk+tn) - 2*tn/tk
        g_t += (s*term).real
    return g_t

def _ice_dgdp(temp,pres):
    """Calculate ice Gibbs free energy P-derivative.
    
    Calculate the derivative of the specific Gibbs free energy of ice
    with respect to pressure.
    
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :returns: Gibbs free energy derivative in J/kg/Pa.
    """
    # Reduced variables
    tn = temp/_TTP
    pn = pres/_PTP
    _PI0 = _PATM/_PTP
    g_p = 0.
    
    # Power series and entropy components
    for (l,n) in enumerate(_GCOEFFS[0]):
        if l > 0: g_p += n * l*(pn-_PI0)**(l-1) / _PTP
    
    # Residual terms including complex numbers
    s = complex(0.0,0.0)
    for (k,rk) in enumerate(_GCOEFFS[2]):
        if k > 0:
            s += rk * k*(pn-_PI0)**(k-1) / _PTP
    tk = _GCOEFFS[3][1]
    term = (tk-tn)*cmath.log(tk-tn) + (tk+tn)*cmath.log(tk+tn)
    term -= 2*tk*cmath.log(tk) + tn**2/tk
    g_p += _TTP * (s*term).real
    return g_p

def _ice_d2gdt2(temp,pres):
    """Calculate ice Gibbs free energy TT-derivative.
    
    Calculate the second derivative of the specific Gibbs free energy of
    ice with respect to temperature.
    
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :returns: Gibbs free energy derivative in J/kg/K^2.
    """
    # Reduced variables
    tn = temp/_TTP
    pn = pres/_PTP
    _PI0 = _PATM/_PTP
    g_tt = 0.
    
    # Residual terms including complex numbers
    sr = [_GCOEFFS[1], complex(0.0,0.0)]
    for (k,rk) in enumerate(_GCOEFFS[2]):
        sr[1] += rk * (pn-_PI0)**k
    for (tk,s) in zip(_GCOEFFS[3],sr):
        term = 1./(tk-tn) + 1./(tk+tn) - 2./tk
        g_tt += (s*term).real / _TTP
    return g_tt

def _ice_d2gdtdp(temp,pres):
    """Calculate ice Gibbs free energy TP-derivative.
    
    Calculate the mixed derivative of the specific Gibbs free energy of
    ice with respect to temperature and pressure.
    
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :returns: Gibbs free energy derivative in J/kg/K/Pa.
    """
    # Reduced variables
    tn = temp/_TTP
    pn = pres/_PTP
    _PI0 = _PATM/_PTP
    
    # Residual terms including complex numbers
    s = complex(0.0,0.0)
    for (k,rk) in enumerate(_GCOEFFS[2]):
        if k > 0:
            s += rk * k*(pn-_PI0)**(k-1) / _PTP
    tk = _GCOEFFS[3][1]
    term = -cmath.log(tk-tn) + cmath.log(tk+tn) - 2*tn/tk
    g_tp = (s*term).real
    return g_tp

def _ice_d2gdp2(temp,pres):
    """Calculate ice Gibbs free energy PP-derivative.
    
    Calculate the second derivative of the specific Gibbs free energy of ice with respect to pressure.
    
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :returns: Gibbs free energy derivative in J/kg/Pa^2.
    """
    # Reduced variables
    tn = temp/_TTP
    pn = pres/_PTP
    _PI0 = _PATM/_PTP
    g_pp = 0.
    
    # Power series and entropy components
    for (l,n) in enumerate(_GCOEFFS[0]):
        if l > 1:
            g_pp += n * l*(l-1) * (pn-_PI0)**(l-2) / _PTP**2
    
    # Residual terms including complex numbers
    s = _GCOEFFS[2][2] * 2. / _PTP**2
    tk = _GCOEFFS[3][1]
    term = (tk-tn)*cmath.log(tk-tn) + (tk+tn)*cmath.log(tk+tn)
    term -= 2*tk*cmath.log(tk) + tn**2/tk
    g_pp += _TTP * (s*term).real
    return g_pp

def _ice_d3gdt3(temp,pres):
    """Calculate ice Gibbs free energy TTT-derivative.
    
    Calculate the third derivative of the specific Gibbs free energy of
    ice with respect to temperature.
    
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :returns: Gibbs free energy derivative in J/kg/K^3.
    """
    # Reduced variables
    tn = temp/_TTP
    pn = pres/_PTP
    _PI0 = _PATM/_PTP
    g_ttt = 0.
    
    # Residual terms including complex numbers
    sr = [_GCOEFFS[1], complex(0.0,0.0)]
    for (k,rk) in enumerate(_GCOEFFS[2]):
        sr[1] += rk * (pn-_PI0)**k
    for (tk,s) in zip(_GCOEFFS[3],sr):
        term = 1./(tk-tn)**2 - 1./(tk+tn)**2
        g_ttt += (s * term).real / _TTP**2
    return g_ttt

def _ice_d3gdt2dp(temp,pres):
    """Calculate ice Gibbs free energy TTP-derivative.
    
    Calculate the mixed temperature-temperature-pressure derivative of
    the specific Gibbs free energy of ice.
    
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :returns: Gibbs free energy derivative in J/kg/K^2/Pa.
    """
    # Reduced variables
    tn = temp/_TTP
    pn = pres/_PTP
    _PI0 = _PATM/_PTP
    
    # Residual terms including complex numbers
    s = complex(0.0,0.0)
    for (k,rk) in enumerate(_GCOEFFS[2]):
        if k > 0:
            s += rk * k*(pn-_PI0)**(k-1) / _PTP
    tk = _GCOEFFS[3][1]
    term = 1./(tk-tn) + 1./(tk+tn) - 2./tk
    g_ttp = (s*term).real / _TTP
    return g_ttp

def _ice_d3gdtdp2(temp,pres):
    """Calculate ice Gibbs free energy TPP-derivative.
    
    Calculate the temperature-pressure-pressure derivative of the
    specific Gibbs free energy of ice.
    
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :returns: Gibbs free energy derivative in J/kg/K/Pa^2.
    """
    # Reduced variables
    tn = temp/_TTP
    pn = pres/_PTP
    _PI0 = _PATM/_PTP
    
    # Residual terms including complex numbers
    s = _GCOEFFS[2][2] * 2.0/_PTP**2
    tk = _GCOEFFS[3][1]
    term = -cmath.log(tk-tn) + cmath.log(tk+tn) - 2.0*tn/tk
    g_tpp = (s*term).real
    return g_tpp

def _ice_d3gdp3(temp,pres):
    """Calculate ice Gibbs free energy PPP-derivative.
    
    Calculate the third derivative of the specific Gibbs free energy of
    ice with respect to pressure.
    
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :returns: Gibbs free energy derivative in J/kg/Pa^3.
    """
    
    # Reduced variables
    tn = temp/_TTP
    pn = pres/_PTP
    _PI0 = _PATM/_PTP
    g_ppp = 0.
    
    # Power series terms
    for (l,n) in enumerate(_GCOEFFS[0]):
        if l > 2:
            g_ppp += n * l*(l-1)*(l-2) * (pn-_PI0)**(l-3)
    g_ppp /= _PTP**3
    return g_ppp


## Public functions
def ice_g(drvt,drvp,temp,pres,chkbnd=True,stacklevel=2):
    """Calculate ice Gibbs free energy.
    
    Calculate the specific Gibbs free energy of ice or its derivatives
    with respect to temperature and pressure. Derivatives up to third
    order are available.
    
    :arg int drvt: Number of temperature derivatives.
    :arg int drvp: Number of pressure derivatives.
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg bool chkbnd: If True (default) then warnings are raised when
        the given values are valid but outside the recommended bounds.
    :arg int stacklevel: Controls how many levels deep to raise the
        warning from (default 2).
    :returns: Gibbs free energy, in units of
        (J/kg) / K^drvt / Pa^drvp.
    :raises ValueError: If temp or pres are nonpositive.
    :raises RuntimeWarning: If temp or pres are outside the recommended
        bounds and chkbnd is True.
    :raises ValueError: if drvt<0, drvp<0, or (drvt+drvp)>3.
    
    :Examples:
    
    >>> ice_g(0,0,270.,1e5)
    -3786.74963128
    >>> ice_g(1,0,270.,1e5)
    1244.97335506
    >>> ice_g(0,1,270.,1e5)
    1.09029713624e-03
    >>> ice_g(2,0,270.,1e5)
    -7.67955356337
    >>> ice_g(1,1,270.,1e5)
    1.72604208697e-07
    >>> ice_g(0,2,270.,1e5)
    -1.27811223643e-13
    """
    constants0.chkicebnds(temp,pres,chkbnd=chkbnd,stacklevel=stacklevel)
    
    # Run through various derivative cases
    if (drvt,drvp) == (0,0):
        g = _ice_g0(temp,pres)
    elif (drvt,drvp) == (1,0):
        g = _ice_dgdt(temp,pres)
    elif (drvt,drvp) == (0,1):
        g = _ice_dgdp(temp,pres)
    elif (drvt,drvp) == (2,0):
        g = _ice_d2gdt2(temp,pres)
    elif (drvt,drvp) == (1,1):
        g = _ice_d2gdtdp(temp,pres)
    elif (drvt,drvp) == (0,2):
        g = _ice_d2gdp2(temp,pres)
    elif (drvt,drvp) == (3,0):
        g = _ice_d3gdt3(temp,pres)
    elif (drvt,drvp) == (2,1):
        g = _ice_d3gdt2dp(temp,pres)
    elif (drvt,drvp) == (1,2):
        g = _ice_d3gdtdp2(temp,pres)
    elif (drvt,drvp) == (0,3):
        g = _ice_d3gdp3(temp,pres)
    else:
        errmsg = 'Derivatives {0} not recognized'.format((drvt,drvp))
        raise ValueError(errmsg)
    
    return g


