"""Icy air Gibbs energy and related properties.

This module provides the Gibbs function for ice-saturated (icy) air and
related thermodynamic properties. The primary variables are the total
dry air fraction, temperature, and pressure. The 'total' fraction here
is the mass fraction of dry air in the total parcel (including ice), and
uses the variable ``wair``. The dry air mass fraction in humid air uses
the variable ``airf``.

:Examples:

>>> iceair_g(0,0,0,.5,270.,1e5)
-2595.57166634
>>> iceair_g(1,0,0,.5,270.,1e5)
2382.35592988
>>> iceair_g(0,0,1,.5,270.,1e5)
0.389645501224
>>> iceair_g(1,1,0,.5,270.,1e5)
-1269.41767949
>>> cp(0.5,270.,1e5)
1893.230554
>>> enthalpy(0.5,270.,1e5)
-167366.990802
>>> lapserate(0.5,270.,1e5)
2.283832444e-04
>>> solidfraction(0.5,270.,1e5)
0.498525089434

:Functions:

* :func:`iceair_g`: Icy air Gibbs free energy with derivatives.
* :func:`cp`: Icy air isobaric heat capacity.
* :func:`density`: Icy air density.
* :func:`enthalpy`: Icy air enthalpy.
* :func:`entropy`: Icy air entropy.
* :func:`expansion`: Icy air thermal expansion coefficient.
* :func:`kappa_t`: Icy air isothermal compressibility.
* :func:`lapserate`: Icy air adiabatic lapse rate.
* :func:`solidfraction`: Total mass fraction of ice in icy air.
* :func:`vapourfraction`: Total mass fraction of water vapour in icy
  air.

"""

__all__ = ['iceair_g','cp','density','enthalpy','entropy','expansion','kappa_t',
    'lapserate','solidfraction','vapourfraction']

import numpy
from teospy import constants0
from teospy import ice1
from teospy import air2
from teospy import ice2
from teospy import maths3
from teospy import air3a
from teospy import iceair4a

_CHKTOL = constants0.CHKTOL
_chkhumbnds = constants0.chkhumbnds
_chkicebnds = constants0.chkicebnds
_ice_g = ice1.ice_g
_air_f = air2.air_f
_eq_pressure = air2.eq_pressure
_eq_vappot = air2.eq_vappot
_newton = maths3.newton
_eq_atpe = iceair4a.eq_atpe


## Gibbs function
def iceair_g(drvw,drvt,drvp,wair,temp,pres,airf=None,dhum=None,
    chkvals=False,chktol=_CHKTOL,airf0=None,dhum0=None,chkbnd=False,
    mathargs=None):
    """Calculate icy air Gibbs free energy with derivatives.
    
    Calculate the specific Gibbs free energy of icy air or its
    derivatives with respect to total dry air fraction, temperature,
    and pressure.
    
    :arg int drvw: Number of dry fraction derivatives.
    :arg int drvt: Number of temperature derivatives.
    :arg int drvp: Number of pressure derivatives.
    :arg float wair: Total dry air fraction in kg/kg.
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg airf: Dry air fraction in humid air in kg/kg.
    :type airf: float or None
    :arg dhum: Humid air density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dhum: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg airf0: Initial guess for the dry fraction in kg/kg. If None
        (default) then `iceair4a._approx_tp` is used.
    :type airf0: float or None
    :arg dhum0: Initial guess for the humid air density in kg/m3. If
        None (default) then `iceair4a._approx_tp` is used.
    :type dhum0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Gibbs free energy in units of
        (J/kg) / (kg/kg)^drvw / K^drvt / Pa^drvp.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    :raises RuntimeWarning: If air with the given parameters would be
        unsaturated.
    
    :Examples:
    
    >>> iceair_g(0,0,0,.5,270.,1e5)
    -2595.57166634
    >>> iceair_g(1,0,0,.5,270.,1e5)
    2382.35592988
    >>> iceair_g(0,1,0,.5,270.,1e5)
    610.264515318
    >>> iceair_g(0,0,1,.5,270.,1e5)
    0.389645501224
    >>> iceair_g(2,0,0,.5,270.,1e5)
    0.0
    >>> iceair_g(1,1,0,.5,270.,1e5)
    -1269.41767949
    >>> iceair_g(1,0,1,.5,270.,1e5)
    0.777110408175
    >>> iceair_g(0,2,0,.5,270.,1e5)
    -7.011965016
    >>> iceair_g(0,1,1,.5,270.,1e5)
    1.601415320e-03
    >>> iceair_g(0,0,2,.5,270.,1e5)
    -3.911839890e-06
    """
    airf, __, __, dhum = _eq_atpe(temp=temp,pres=pres,airf=airf,dhum=dhum,
        chkvals=chkvals,chktol=chktol,airf0=airf0,dhum0=dhum0,chkbnd=chkbnd,
        mathargs=mathargs)
    if airf <= wair:
        warnmsg = 'Air with the given parameters is unsaturated'
        warnings.warn(warnmsg,RuntimeWarning)
        g = air3a.air_g(drvw,drvt,drvp,wair,temp,pres,dhum=dhum)
        return g
    w = wair / airf
    
    # Simple derivative cases
    if (drvw,drvt,drvp) == (0,0,0):
        fh = _air_f(0,0,0,airf,temp,dhum)
        fh_d = _air_f(0,0,1,airf,temp,dhum)
        gi = _ice_g(0,0,temp,pres)
        g = w*(fh + dhum*fh_d) + (1-w)*gi
        return g
    elif (drvw,drvt,drvp) == (1,0,0):
        fh_a = _air_f(1,0,0,airf,temp,dhum)
        g_w = fh_a
        return g_w
    elif (drvw,drvt,drvp) == (0,1,0):
        gh_t = _air_f(0,1,0,airf,temp,dhum)
        gi_t = _ice_g(1,0,temp,pres)
        g_t = w*gh_t + (1-w)*gi_t
        return g_t
    elif (drvw,drvt,drvp) == (0,0,1):
        gh_p = dhum**(-1)
        gi_p = _ice_g(0,1,temp,pres,chkbnd=chkbnd)
        g_p = w*gh_p + (1-w)*gi_p
        return g_p
    elif (drvw,drvt,drvp) == (2,0,0):
        g_ww = 0.
        return g_ww
    elif (drvw,drvt,drvp) == (1,1,0):
        gh_t = _air_f(0,1,0,airf,temp,dhum)
        gi_t = _ice_g(1,0,temp,pres)
        g_wt = (gh_t - gi_t) / airf
        return g_wt
    elif (drvw,drvt,drvp) == (1,0,1):
        gh_p = dhum**(-1)
        gi_p = _ice_g(0,1,temp,pres)
        g_wp = (gh_p - gi_p) / airf
        return g_wp
    
    # Derivative cases requiring inversion
    __, __, __, pg_ad = iceair4a._diff_tp(airf,dhum,temp,pres)
    if (drvw,drvt,drvp) == (0,2,0):
        ph_t = _eq_pressure(0,1,0,airf,temp,dhum)
        gv_t = _eq_vappot(0,1,0,airf,temp,dhum)
        gi_t = _ice_g(1,0,temp,pres)
        pg_t = numpy.array([ph_t,gv_t-gi_t])
        ad_t = numpy.linalg.solve(pg_ad,-pg_t)
        
        fh_t = _air_f(0,1,0,airf,temp,dhum)
        fh_at = _air_f(1,1,0,airf,temp,dhum)
        fh_tt = _air_f(0,2,0,airf,temp,dhum)
        fh_td = _air_f(0,1,1,airf,temp,dhum)
        gi_t = _ice_g(1,0,temp,pres)
        gi_tt = _ice_g(2,0,temp,pres)
        g_ta = w/airf*(airf*fh_at - fh_t + gi_t)
        g_td = w*fh_td
        g_tx = numpy.array([g_ta,g_td])
        g_tt = w*fh_tt + (1-w)*gi_tt + g_tx.dot(ad_t)
        return g_tt
    elif (drvw,drvt,drvp) == (0,1,1):
        ph_t = _eq_pressure(0,1,0,airf,temp,dhum)
        gv_t = _eq_vappot(0,1,0,airf,temp,dhum)
        gi_t = _ice_g(1,0,temp,pres)
        pg_t = numpy.array([ph_t,gv_t-gi_t])
        ad_t = numpy.linalg.solve(pg_ad,-pg_t)
        
        gi_p = _ice_g(0,1,temp,pres)
        gi_tp = _ice_g(1,1,temp,pres)
        g_pa = -w/airf * (dhum**(-1) - gi_p)
        g_pd = -w/dhum**2
        g_px = numpy.array([g_pa,g_pd])
        g_tp = (1-w)*gi_tp + g_px.dot(ad_t)
        return g_tp
    elif (drvw,drvt,drvp) == (0,0,2):
        gi_p = _ice_g(0,1,temp,pres)
        pg_p = numpy.array([1.,gi_p])
        ad_p = numpy.linalg.solve(pg_ad,pg_p)
        
        gi_pp = _ice_g(0,2,temp,pres)
        g_pa = -w/airf * (dhum**(-1) - gi_p)
        g_pd = -w/dhum**2
        g_px = numpy.array([g_pa,g_pd])
        g_pp = (1-w)*gi_pp + g_px.dot(ad_p)
        return g_pp
    
    # Should not have made it this far!
    errmsg = 'Derivatives {0} not recognized'.format((drvw,drvt,drvp))
    raise ValueError(errmsg)


## Thermodynamic properties
def cp(wair,temp,pres,airf=None,dhum=None,chkvals=False,chktol=_CHKTOL,
    airf0=None,dhum0=None,chkbnd=False,mathargs=None):
    """Calculate icy air isobaric heat capacity.
    
    Calculate the isobaric (constant pressure) heat capacity of icy air.
    
    :arg float wair: Total dry air fraction in kg/kg.
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg airf: Dry air fraction in humid air in kg/kg.
    :type airf: float or None
    :arg dhum: Humid air density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dhum: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg airf0: Initial guess for the dry fraction in kg/kg. If None
        (default) then `iceair4a._approx_tp` is used.
    :type airf0: float or None
    :arg dhum0: Initial guess for the humid air density in kg/m3. If
        None (default) then `iceair4a._approx_tp` is used.
    :type dhum0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Heat capacity in J/kg/K.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    :raises RuntimeWarning: If air with the given parameters would be
        unsaturated.
    
    :Examples:
    
    >>> cp(0.5,270.,1e5)
    1893.230554
    """
    g_tt = iceair_g(0,2,0,wair,temp,pres,airf=airf,dhum=dhum,chkvals=chkvals,
        chktol=chktol,airf0=airf0,dhum0=dhum0,chkbnd=chkbnd,mathargs=mathargs)
    cp = -temp * g_tt
    return cp

def density(wair,temp,pres,airf=None,dhum=None,chkvals=False,
    chktol=_CHKTOL,airf0=None,dhum0=None,chkbnd=False,mathargs=None):
    """Calculate icy air density.
    
    Calculate the density of icy air, the total density of the parcel
    including ice mass.
    
    :arg float wair: Total dry air fraction in kg/kg.
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg airf: Dry air fraction in humid air in kg/kg.
    :type airf: float or None
    :arg dhum: Humid air density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dhum: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg airf0: Initial guess for the dry fraction in kg/kg. If None
        (default) then `iceair4a._approx_tp` is used.
    :type airf0: float or None
    :arg dhum0: Initial guess for the humid air density in kg/m3. If
        None (default) then `iceair4a._approx_tp` is used.
    :type dhum0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Density in kg/m3.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    :raises RuntimeWarning: If air with the given parameters would be
        unsaturated.
    
    :Examples:
    
    >>> density(0.5,270.,1e5)
    2.56643538000
    """
    g_p = iceair_g(0,0,1,wair,temp,pres,airf=airf,dhum=dhum,chkvals=chkvals,
        chktol=chktol,airf0=airf0,dhum0=dhum0,chkbnd=chkbnd,mathargs=mathargs)
    rho = g_p**(-1)
    return rho

def enthalpy(wair,temp,pres,airf=None,dhum=None,chkvals=False,
    chktol=_CHKTOL,airf0=None,dhum0=None,chkbnd=False,mathargs=None):
    """Calculate icy air enthalpy.
    
    Calculate the specific enthalpy of icy air.
    
    :arg float wair: Total dry air fraction in kg/kg.
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg airf: Dry air fraction in humid air in kg/kg.
    :type airf: float or None
    :arg dhum: Humid air density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dhum: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg airf0: Initial guess for the dry fraction in kg/kg. If None
        (default) then `iceair4a._approx_tp` is used.
    :type airf0: float or None
    :arg dhum0: Initial guess for the humid air density in kg/m3. If
        None (default) then `iceair4a._approx_tp` is used.
    :type dhum0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Enthalpy in J/kg.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    :raises RuntimeWarning: If air with the given parameters would be
        unsaturated.
    
    :Examples:
    
    >>> enthalpy(0.5,270.,1e5)
    -167366.990802
    """
    airf, __, __, dhum = _eq_atpe(temp=temp,pres=pres,airf=airf,dhum=dhum,
        chkvals=chkvals,chktol=chktol,airf0=airf0,dhum0=dhum0,chkbnd=chkbnd,
        mathargs=mathargs)
    if airf <= wair:
        warnmsg = 'Air with the given parameters is unsaturated'
        warnings.warn(warnmsg,RuntimeWarning)
        h = air3b.enthalpy(wair,temp,pres,dhum0=dhum0,mathargs=mathargs)
        return h
    g = iceair_g(0,0,0,wair,temp,pres,airf=airf,dhum=dhum)
    g_t = iceair_g(0,1,0,wair,temp,pres,airf=airf,dhum=dhum)
    h = g - temp*g_t
    return h

def entropy(wair,temp,pres,airf=None,dhum=None,chkvals=False,
    chktol=_CHKTOL,airf0=None,dhum0=None,chkbnd=False,mathargs=None):
    """Calculate icy air entropy.
    
    Calculate the specific entropy of icy air.
    
    :arg float wair: Total dry air fraction in kg/kg.
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg airf: Dry air fraction in humid air in kg/kg.
    :type airf: float or None
    :arg dhum: Humid air density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dhum: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg airf0: Initial guess for the dry fraction in kg/kg. If None
        (default) then `iceair4a._approx_tp` is used.
    :type airf0: float or None
    :arg dhum0: Initial guess for the humid air density in kg/m3. If
        None (default) then `iceair4a._approx_tp` is used.
    :type dhum0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Entropy in J/kg/K.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    :raises RuntimeWarning: If air with the given parameters would be
        unsaturated.
    
    :Examples:
    
    >>> entropy(0.5,270.,1e5)
    -610.264515318
    """
    g_t = iceair_g(0,1,0,wair,temp,pres,airf=airf,dhum=dhum,chkvals=chkvals,
        chktol=chktol,airf0=airf0,dhum0=dhum0,chkbnd=chkbnd,mathargs=mathargs)
    s = -g_t
    return s

def expansion(wair,temp,pres,airf=None,dhum=None,chkvals=False,
    chktol=_CHKTOL,airf0=None,dhum0=None,chkbnd=False,mathargs=None):
    """Calculate icy air expansion coefficient.
    
    Calculate the thermal expansion coefficient of icy air.
    
    :arg float wair: Total dry air fraction in kg/kg.
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg airf: Dry air fraction in humid air in kg/kg.
    :type airf: float or None
    :arg dhum: Humid air density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dhum: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg airf0: Initial guess for the dry fraction in kg/kg. If None
        (default) then `iceair4a._approx_tp` is used.
    :type airf0: float or None
    :arg dhum0: Initial guess for the humid air density in kg/m3. If
        None (default) then `iceair4a._approx_tp` is used.
    :type dhum0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Expansion coefficient in J/kg/K.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    :raises RuntimeWarning: If air with the given parameters would be
        unsaturated.
    
    :Examples:
    
    >>> expansion(0.5,270.,1e5)
    4.109928935e-03
    """
    airf, __, __, dhum = _eq_atpe(temp=temp,pres=pres,airf=airf,dhum=dhum,
        chkvals=chkvals,chktol=chktol,airf0=airf0,dhum0=dhum0,chkbnd=chkbnd,
        mathargs=mathargs)
    if airf <= wair:
        warnmsg = 'Air with the given parameters is unsaturated'
        warnings.warn(warnmsg,RuntimeWarning)
        alpha = air3b.expansion(wair,temp,pres,dhum0=dhum0,mathargs=mathargs)
        return alpha
    g_p = iceair_g(0,0,1,wair,temp,pres,airf=airf,dhum=dhum)
    g_tp = iceair_g(0,1,1,wair,temp,pres,airf=airf,dhum=dhum)
    alpha = g_tp / g_p
    return alpha

def kappa_t(wair,temp,pres,airf=None,dhum=None,chkvals=False,
    chktol=_CHKTOL,airf0=None,dhum0=None,chkbnd=False,mathargs=None):
    """Calculate icy air isothermal compressibility.
    
    Calculate the isothermal compressibility of icy air.
    
    :arg float wair: Total dry air fraction in kg/kg.
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg airf: Dry air fraction in humid air in kg/kg.
    :type airf: float or None
    :arg dhum: Humid air density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dhum: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg airf0: Initial guess for the dry fraction in kg/kg. If None
        (default) then `iceair4a._approx_tp` is used.
    :type airf0: float or None
    :arg dhum0: Initial guess for the humid air density in kg/m3. If
        None (default) then `iceair4a._approx_tp` is used.
    :type dhum0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Compressibility in 1/Pa.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    :raises RuntimeWarning: If air with the given parameters would be
        unsaturated.
    
    :Examples:
    
    >>> kappa_t(0.5,270.,1e5)
    1.003948429e-05
    """
    airf, __, __, dhum = _eq_atpe(temp=temp,pres=pres,airf=airf,dhum=dhum,
        chkvals=chkvals,chktol=chktol,airf0=airf0,dhum0=dhum0,chkbnd=chkbnd,
        mathargs=mathargs)
    if airf <= wair:
        warnmsg = 'Air with the given parameters is unsaturated'
        warnings.warn(warnmsg,RuntimeWarning)
        kappa = air3b.kappa_t(wair,temp,pres,dhum0=dhum0,mathargs=mathargs)
        return kappa
    g_p = iceair_g(0,0,1,wair,temp,pres,airf=airf,dhum=dhum)
    g_pp = iceair_g(0,0,2,wair,temp,pres,airf=airf,dhum=dhum)
    kappa = -g_pp / g_p
    return kappa

def lapserate(wair,temp,pres,airf=None,dhum=None,chkvals=False,
    chktol=_CHKTOL,airf0=None,dhum0=None,chkbnd=False,mathargs=None):
    """Calculate icy air adiabatic lapse rate.
    
    Calculate the adiabatic lapse rate of icy air.
    
    :arg float wair: Total dry air fraction in kg/kg.
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg airf: Dry air fraction in humid air in kg/kg.
    :type airf: float or None
    :arg dhum: Humid air density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dhum: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg airf0: Initial guess for the dry fraction in kg/kg. If None
        (default) then `iceair4a._approx_tp` is used.
    :type airf0: float or None
    :arg dhum0: Initial guess for the humid air density in kg/m3. If
        None (default) then `iceair4a._approx_tp` is used.
    :type dhum0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Lapse rate in K/Pa.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    :raises RuntimeWarning: If air with the given parameters would be
        unsaturated.
    
    :Examples:
    
    >>> lapserate(0.5,270.,1e5)
    2.283832444e-04
    """
    airf, __, __, dhum = _eq_atpe(temp=temp,pres=pres,airf=airf,dhum=dhum,
        chkvals=chkvals,chktol=chktol,airf0=airf0,dhum0=dhum0,chkbnd=chkbnd,
        mathargs=mathargs)
    if airf <= wair:
        warnmsg = 'Air with the given parameters is unsaturated'
        warnings.warn(warnmsg,RuntimeWarning)
        gamma = air3b.lapserate(wair,temp,pres,dhum0=dhum0,mathargs=mathargs)
        return gamma
    g_tt = iceair_g(0,2,0,wair,temp,pres,airf=airf,dhum=dhum)
    g_tp = iceair_g(0,1,1,wair,temp,pres,airf=airf,dhum=dhum)
    gamma = -g_tp / g_tt
    return gamma

def solidfraction(wair,temp,pres,airf=None,dhum=None,chkvals=False,
    chktol=_CHKTOL,airf0=None,dhum0=None,chkbnd=False,mathargs=None):
    """Calculate icy air ice fraction.
    
    Calculate the mass fraction of ice in icy air.
    
    :arg float wair: Total dry air fraction in kg/kg.
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg airf: Dry air fraction in humid air in kg/kg.
    :type airf: float or None
    :arg dhum: Humid air density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dhum: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg airf0: Initial guess for the dry fraction in kg/kg. If None
        (default) then `iceair4a._approx_tp` is used.
    :type airf0: float or None
    :arg dhum0: Initial guess for the humid air density in kg/m3. If
        None (default) then `iceair4a._approx_tp` is used.
    :type dhum0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Mass fraction in kg/kg.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    :raises RuntimeWarning: If air with the given parameters would be
        unsaturated.
    
    :Examples:
    
    >>> solidfraction(0.5,270.,1e5)
    0.498525089434
    """
    airf, __, __, dhum = _eq_atpe(temp=temp,pres=pres,airf=airf,dhum=dhum,
        chkvals=chkvals,chktol=chktol,airf0=airf0,dhum0=dhum0,chkbnd=chkbnd,
        mathargs=mathargs)
    if airf <= wair:
        warnmsg = 'Air with the given parameters is unsaturated'
        warnings.warn(warnmsg,RuntimeWarning)
    wice = max(1 - wair/airf, 0.)
    return wice

def vapourfraction(wair,temp,pres,airf=None,dhum=None,chkvals=False,
    chktol=_CHKTOL,airf0=None,dhum0=None,chkbnd=False,mathargs=None):
    """Calculate icy air vapour fraction.
    
    Calculate the mass fraction of water vapour in icy air.
    
    :arg float wair: Total dry air fraction in kg/kg.
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg airf: Dry air fraction in humid air in kg/kg.
    :type airf: float or None
    :arg dhum: Humid air density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dhum: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg airf0: Initial guess for the dry fraction in kg/kg. If None
        (default) then `iceair4a._approx_tp` is used.
    :type airf0: float or None
    :arg dhum0: Initial guess for the humid air density in kg/m3. If
        None (default) then `iceair4a._approx_tp` is used.
    :type dhum0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Mass fraction in kg/kg.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    :raises RuntimeWarning: If air with the given parameters would be
        unsaturated.
    
    :Examples:
    
    >>> vapourfraction(0.5,270.,1e5)
    1.47491056602e-3
    """
    airf, __, __, dhum = _eq_atpe(temp=temp,pres=pres,airf=airf,dhum=dhum,
        chkvals=chkvals,chktol=chktol,airf0=airf0,dhum0=dhum0,chkbnd=chkbnd,
        mathargs=mathargs)
    if airf <= wair:
        warnmsg = 'Air with the given parameters is unsaturated'
        warnings.warn(warnmsg,RuntimeWarning)
    wvap = min(wair * (1-airf)/airf, 1-wair)
    return wvap

