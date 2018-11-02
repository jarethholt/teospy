"""Wet air Gibbs energy and related properties.

This module provides the Gibbs function for liquid water-saturated (wet)
air and related thermodynamic quantities. The primary variables are the
total dry air fraction, temperature, and pressure. The 'total' fraction
here is the mass fraction of dry air in the total parcel (including
liquid) and uses the variable ``wair``. The dry air mass fraction in
humid air uses the variable ``airf``.

:Examples:

>>> liqair_g(0,0,0,0.5,300.,1e5)
-5396.77820137
>>> liqair_g(0,0,1,0.5,300.,1e5)
0.446729465555
>>> liqair_g(0,1,1,0.5,300.,1e5)
2.45335972867e-03
>>> cp(0.5,300.,1e5)
4267.95671050
>>> expansion(0.5,300.,1e5)
5.49182428703e-03
>>> lapserate(0.5,300.,1e5)
1.72449715057e-04

:Functions:

* :func:`liqair_g`: Wet air Gibbs free energy with derivatives.
* :func:`cp`: Wet air isobaric heat capacity.
* :func:`density`: Wet air density.
* :func:`enthalpy`: Wet air enthalpy.
* :func:`entropy`: Wet air entropy.
* :func:`expansion`: Wet air thermal expansion coefficient.
* :func:`kappa_t`: Wet air isothermal compressibility.
* :func:`lapserate`: Wet air adiabatic lapse rate.
* :func:`liquidfraction`: Total mass fraction of liquid water in wet
  air.
* :func:`vapourfraction`: Total mass fraction of water vapour in wet
  air.

"""

__all__ = ['liqair_g','cp','density','enthalpy','entropy','expansion','kappa_t',
    'lapserate','liquidfraction','vapourfraction']

import numpy
import warnings
import constants0
import flu1
import air2
import flu2
import maths3
import air3a
import liqair4a

_CHKTOL = constants0.CHKTOL
_chkhumbnds = constants0.chkhumbnds
_chkflubnds = constants0.chkflubnds
_flu_f = flu1.flu_f
_air_f = air2.air_f
_air_eq_pressure = air2.eq_pressure
_air_eq_vappot = air2.eq_vappot
_flu_eq_pressure = flu2.eq_pressure
_flu_eq_chempot = flu2.eq_chempot
_newton = maths3.newton
_eq_atpe = liqair4a.eq_atpe


## Gibbs function
def liqair_g(drvw,drvt,drvp,wair,temp,pres,airf=None,dhum=None,
    dliq=None,chkvals=False,chktol=_CHKTOL,airf0=None,dhum0=None,
    dliq0=None,chkbnd=False,mathargs=None):
    """Calculate wet air Gibbs free energy with derivatives.
    
    Calculate the specific Gibbs free energy of wet air or its
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
    :arg dliq: Liquid water density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dliq: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg airf0: Initial guess for the dry fraction in kg/kg. If None
        (default) then `iceair4a._approx_tp` is used.
    :type airf0: float or None
    :arg dhum0: Initial guess for the humid air density in kg/m3. If
        None (default) then `liqair4a._approx_tp` is used.
    :type dhum0: float or None
    :arg dliq0: Initial guess for the liquid water density in kg/m3. If
        None (default) then `liqair4a._approx_tp` is used.
    :type dliq0: float or None
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
    
    >>> liqair_g(0,0,0,0.5,300.,1e5)
    -5396.77820137
    >>> liqair_g(1,0,0,0.5,300.,1e5)
    -263.4554912
    >>> liqair_g(0,1,0,0.5,300.,1e5)
    -343.783393872
    >>> liqair_g(0,0,1,0.5,300.,1e5)
    0.446729465555
    >>> liqair_g(2,0,0,0.5,300.,1e5)
    0.
    >>> liqair_g(1,1,0,0.5,300.,1e5)
    98.5580798842
    >>> liqair_g(1,0,1,0.5,300.,1e5)
    0.891452019991
    >>> liqair_g(0,2,0,0.5,300.,1e5)
    -14.2265223683
    >>> liqair_g(0,1,1,0.5,300.,1e5)
    2.45335972867e-03
    >>> liqair_g(0,0,2,0.5,300.,1e5)
    -4.62725155875e-06
    """
    airf, __, __, dhum, dliq = _eq_atpe(temp=temp,pres=pres,airf=airf,
        dhum=dhum,dliq=dliq,chkvals=chkvals,chktol=chktol,airf0=airf0,
        dhum0=dhum0,dliq0=dliq0,chkbnd=chkbnd,mathargs=mathargs)
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
        fl = _flu_f(0,0,temp,dliq)
        fl_d = _flu_f(0,1,temp,dliq)
        g = w*(fh + dhum*fh_d) + (1-w)*(fl + dliq*fl_d)
        return g
    elif (drvw,drvt,drvp) == (1,0,0):
        fh_a = _air_f(1,0,0,airf,temp,dhum)
        g_w = fh_a
        return g_w
    elif (drvw,drvt,drvp) == (0,1,0):
        fh_t = _air_f(0,1,0,airf,temp,dhum)
        fl_t = _flu_f(1,0,temp,dliq)
        g_t = w*fh_t + (1-w)*fl_t
        return g_t
    elif (drvw,drvt,drvp) == (0,0,1):
        g_p = w/dhum + (1-w)/dliq
        return g_p
    elif (drvw,drvt,drvp) == (2,0,0):
        g_ww = 0.
        return g_ww
    elif (drvw,drvt,drvp) == (1,1,0):
        fh_t = _air_f(0,1,0,airf,temp,dhum)
        fl_t = _flu_f(1,0,temp,dliq)
        g_wt = (fh_t - fl_t) / airf
        return g_wt
    elif (drvw,drvt,drvp) == (1,0,1):
        g_wp = (dhum**(-1) - dliq**(-1)) / airf
        return g_wp
    
    # Higher-order derivatives require inversion
    __, __, dlhs, drhs = liqair4a._diff_tp(airf,dhum,dliq,temp,pres)
    ppg_x = drhs - dlhs
    
    if (drvw,drvt,drvp) == (0,2,0):
        ph_t = _air_eq_pressure(0,1,0,airf,temp,dhum)
        pl_t = _flu_eq_pressure(1,0,temp,dliq)
        muv_t = _air_eq_vappot(0,1,0,airf,temp,dhum)
        gl_t = _flu_eq_chempot(1,0,temp,dliq)
        ppg_t = numpy.array([ph_t,pl_t,muv_t-gl_t])
        x_t = numpy.linalg.solve(ppg_x,-ppg_t)
        
        fh_t = _air_f(0,1,0,airf,temp,dhum)
        fh_at = _air_f(1,1,0,airf,temp,dhum)
        fh_tt = _air_f(0,2,0,airf,temp,dhum)
        fh_td = _air_f(0,1,1,airf,temp,dhum)
        fl_t = _flu_f(1,0,temp,dliq)
        fl_tt = _flu_f(2,0,temp,dliq)
        fl_td = _flu_f(1,1,temp,dliq)
        g_ta = -w/airf*(fh_t - airf*fh_at - fl_t)
        g_th = w*fh_td
        g_tl = (1-w)*fl_td
        g_tx = numpy.array([g_ta,g_th,g_tl])
        g_tt = w*fh_tt + (1-w)*fl_tt + g_tx.dot(x_t)
        return g_tt
    elif (drvw,drvt,drvp) == (0,1,1):
        ppg_p = numpy.array([1.,1.,0.])
        x_p = numpy.linalg.solve(ppg_x,ppg_p)
        
        fh_t = _air_f(0,1,0,airf,temp,dhum)
        fh_at = _air_f(1,1,0,airf,temp,dhum)
        fh_td = _air_f(0,1,1,airf,temp,dhum)
        fl_t = _flu_f(1,0,temp,dliq)
        fl_td = _flu_f(1,1,temp,dliq)
        g_ta = -w/airf*(fh_t - airf*fh_at - fl_t)
        g_th = w*fh_td
        g_tl = (1-w)*fl_td
        g_tx = numpy.array([g_ta,g_th,g_tl])
        g_tp = g_tx.dot(x_p)
        return g_tp
    elif (drvw,drvt,drvp) == (0,0,2):
        ppg_p = numpy.array([1.,1.,0.])
        x_p = numpy.linalg.solve(ppg_x,ppg_p)
        
        g_pa = -w/airf*(dhum**(-1) - dliq**(-1))
        g_ph = -w/dhum**2
        g_pl = -(1-w)/dliq**2
        g_px = numpy.array([g_pa,g_ph,g_pl])
        g_pp = g_px.dot(x_p)
        return g_pp
    
    # Should not have made it this far!
    errmsg = 'Derivatives {0} not recognized'.format((drvw,drvt,drvp))
    raise ValueError(errmsg)


## Thermodynamic properties
def cp(wair,temp,pres,airf=None,dhum=None,dliq=None,chkvals=False,
    chktol=_CHKTOL,airf0=None,dhum0=None,dliq0=None,chkbnd=False,
    mathargs=None):
    """Calculate wet air isobaric heat capacity.
    
    Calculate the isobaric heat capacity of wet air.
    
    :arg float wair: Total dry air fraction in kg/kg.
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg airf: Dry air fraction in humid air in kg/kg.
    :type airf: float or None
    :arg dhum: Humid air density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dhum: float or None
    :arg dliq: Liquid water density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dliq: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg airf0: Initial guess for the dry fraction in kg/kg. If None
        (default) then `iceair4a._approx_tp` is used.
    :type airf0: float or None
    :arg dhum0: Initial guess for the humid air density in kg/m3. If
        None (default) then `liqair4a._approx_tp` is used.
    :type dhum0: float or None
    :arg dliq0: Initial guess for the liquid water density in kg/m3. If
        None (default) then `liqair4a._approx_tp` is used.
    :type dliq0: float or None
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
    
    >>> cp(0.5,300.,1e5)
    4267.95671050
    """
    g_tt = liqair_g(0,2,0,wair,temp,pres,airf=airf,dhum=dhum,dliq=dliq,
        chkvals=chkvals,chktol=chktol,airf0=airf0,dhum0=dhum0,dliq0=dliq0,
        chkbnd=chkbnd,mathargs=mathargs)
    cp = -temp * g_tt
    return cp

def density(wair,temp,pres,airf=None,dhum=None,dliq=None,chkvals=False,
    chktol=_CHKTOL,airf0=None,dhum0=None,dliq0=None,chkbnd=False,
    mathargs=None):
    """Calculate wet air density.
    
    Calculate the density of wet air.
    
    :arg float wair: Total dry air fraction in kg/kg.
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg airf: Dry air fraction in humid air in kg/kg.
    :type airf: float or None
    :arg dhum: Humid air density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dhum: float or None
    :arg dliq: Liquid water density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dliq: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg airf0: Initial guess for the dry fraction in kg/kg. If None
        (default) then `iceair4a._approx_tp` is used.
    :type airf0: float or None
    :arg dhum0: Initial guess for the humid air density in kg/m3. If
        None (default) then `liqair4a._approx_tp` is used.
    :type dhum0: float or None
    :arg dliq0: Initial guess for the liquid water density in kg/m3. If
        None (default) then `liqair4a._approx_tp` is used.
    :type dliq0: float or None
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
    
    >>> density(0.5,300.,1e5)
    2.23849125053
    """
    g_p = liqair_g(0,0,1,wair,temp,pres,airf=airf,dhum=dhum,dliq=dliq,
        chkvals=chkvals,chktol=chktol,airf0=airf0,dhum0=dhum0,dliq0=dliq0,
        chkbnd=chkbnd,mathargs=mathargs)
    dtot = g_p**(-1)
    return dtot

def enthalpy(wair,temp,pres,airf=None,dhum=None,dliq=None,chkvals=False,
    chktol=_CHKTOL,airf0=None,dhum0=None,dliq0=None,chkbnd=False,
    mathargs=None):
    """Calculate wet air enthalpy.
    
    Calculate the specific enthalpy of wet air.
    
    :arg float wair: Total dry air fraction in kg/kg.
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg airf: Dry air fraction in humid air in kg/kg.
    :type airf: float or None
    :arg dhum: Humid air density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dhum: float or None
    :arg dliq: Liquid water density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dliq: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg airf0: Initial guess for the dry fraction in kg/kg. If None
        (default) then `iceair4a._approx_tp` is used.
    :type airf0: float or None
    :arg dhum0: Initial guess for the humid air density in kg/m3. If
        None (default) then `liqair4a._approx_tp` is used.
    :type dhum0: float or None
    :arg dliq0: Initial guess for the liquid water density in kg/m3. If
        None (default) then `liqair4a._approx_tp` is used.
    :type dliq0: float or None
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
    
    >>> enthalpy(0.5,300.,1e5)
    97738.2399604
    """
    airf, __, __, dhum, dliq = _eq_atpe(temp=temp,pres=pres,airf=airf,
        dhum=dhum,dliq=dliq,chkvals=chkvals,chktol=chktol,airf0=airf0,
        dhum0=dhum0,dliq0=dliq0,chkbnd=chkbnd,mathargs=mathargs)
    if airf <= wair:
        warnmsg = 'Air with the given parameters is unsaturated'
        warnings.warn(warnmsg,RuntimeWarning)
        h = air3b.enthalpy(wair,temp,pres,dhum0=dhum0,mathargs=mathargs)
        return h
    g = liqair_g(0,0,0,wair,temp,pres,airf=airf,dhum=dhum,dliq=dliq)
    g_t = liqair_g(0,1,0,wair,temp,pres,airf=airf,dhum=dhum,dliq=dliq)
    h = g - temp*g_t
    return h

def entropy(wair,temp,pres,airf=None,dhum=None,dliq=None,chkvals=False,
    chktol=_CHKTOL,airf0=None,dhum0=None,dliq0=None,chkbnd=False,
    mathargs=None):
    """Calculate wet air entropy.
    
    Calculate the specific entropy of wet air.
    
    :arg float wair: Total dry air fraction in kg/kg.
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg airf: Dry air fraction in humid air in kg/kg.
    :type airf: float or None
    :arg dhum: Humid air density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dhum: float or None
    :arg dliq: Liquid water density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dliq: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg airf0: Initial guess for the dry fraction in kg/kg. If None
        (default) then `iceair4a._approx_tp` is used.
    :type airf0: float or None
    :arg dhum0: Initial guess for the humid air density in kg/m3. If
        None (default) then `liqair4a._approx_tp` is used.
    :type dhum0: float or None
    :arg dliq0: Initial guess for the liquid water density in kg/m3. If
        None (default) then `liqair4a._approx_tp` is used.
    :type dliq0: float or None
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
    
    >>> entropy(0.5,300.,1e5)
    343.783393872
    """
    g_t = liqair_g(0,1,0,wair,temp,pres,airf=airf,dhum=dhum,dliq=dliq,
        chkvals=chkvals,chktol=chktol,airf0=airf0,dhum0=dhum0,dliq0=dliq0,
        chkbnd=chkbnd,mathargs=mathargs)
    s = -g_t
    return s

def expansion(wair,temp,pres,airf=None,dhum=None,dliq=None,
    chkvals=False,chktol=_CHKTOL,airf0=None,dhum0=None,dliq0=None,
    chkbnd=False,mathargs=None):
    """Calculate wet air thermal expansion coefficient.
    
    Calculate the thermal expansion coefficient of wet air.
    
    :arg float wair: Total dry air fraction in kg/kg.
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg airf: Dry air fraction in humid air in kg/kg.
    :type airf: float or None
    :arg dhum: Humid air density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dhum: float or None
    :arg dliq: Liquid water density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dliq: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg airf0: Initial guess for the dry fraction in kg/kg. If None
        (default) then `iceair4a._approx_tp` is used.
    :type airf0: float or None
    :arg dhum0: Initial guess for the humid air density in kg/m3. If
        None (default) then `liqair4a._approx_tp` is used.
    :type dhum0: float or None
    :arg dliq0: Initial guess for the liquid water density in kg/m3. If
        None (default) then `liqair4a._approx_tp` is used.
    :type dliq0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Expansion coefficient in 1/K.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    :raises RuntimeWarning: If air with the given parameters would be
        unsaturated.
    
    :Examples:
    
    >>> expansion(0.5,300.,1e5)
    5.49182428703e-03
    """
    airf, __, __, dhum, dliq = _eq_atpe(temp=temp,pres=pres,airf=airf,
        dhum=dhum,dliq=dliq,chkvals=chkvals,chktol=chktol,airf0=airf0,
        dhum0=dhum0,dliq0=dliq0,chkbnd=chkbnd,mathargs=mathargs)
    if airf <= wair:
        warnmsg = 'Air with the given parameters is unsaturated'
        warnings.warn(warnmsg,RuntimeWarning)
        alpha = air3b.expansion(wair,temp,pres,dhum0=dhum0,mathargs=mathargs)
        return alpha
    g_p = liqair_g(0,0,1,wair,temp,pres,airf=airf,dhum=dhum,dliq=dliq)
    g_tp = liqair_g(0,1,1,wair,temp,pres,airf=airf,dhum=dhum,dliq=dliq)
    alpha = g_tp / g_p
    return alpha

def kappa_t(wair,temp,pres,airf=None,dhum=None,dliq=None,chkvals=False,
    chktol=_CHKTOL,airf0=None,dhum0=None,dliq0=None,chkbnd=False,
    mathargs=None):
    """Calculate wet air isothermal compressibility.
    
    Calculate the isothermal compressibility of wet air.
    
    :arg float wair: Total dry air fraction in kg/kg.
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg airf: Dry air fraction in humid air in kg/kg.
    :type airf: float or None
    :arg dhum: Humid air density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dhum: float or None
    :arg dliq: Liquid water density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dliq: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg airf0: Initial guess for the dry fraction in kg/kg. If None
        (default) then `iceair4a._approx_tp` is used.
    :type airf0: float or None
    :arg dhum0: Initial guess for the humid air density in kg/m3. If
        None (default) then `liqair4a._approx_tp` is used.
    :type dhum0: float or None
    :arg dliq0: Initial guess for the liquid water density in kg/m3. If
        None (default) then `liqair4a._approx_tp` is used.
    :type dliq0: float or None
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
    
    >>> kappa_t(0.5,300.,1e5)
    1.03580621283e-05
    """
    airf, __, __, dhum, dliq = _eq_atpe(temp=temp,pres=pres,airf=airf,
        dhum=dhum,dliq=dliq,chkvals=chkvals,chktol=chktol,airf0=airf0,
        dhum0=dhum0,dliq0=dliq0,chkbnd=chkbnd,mathargs=mathargs)
    if airf <= wair:
        warnmsg = 'Air with the given parameters is unsaturated'
        warnings.warn(warnmsg,RuntimeWarning)
        kappa = air3b.kappa_t(wair,temp,pres,dhum0=dhum0,mathargs=mathargs)
        return kappa
    g_p = liqair_g(0,0,1,wair,temp,pres,airf=airf,dhum=dhum,dliq=dliq)
    g_pp = liqair_g(0,0,2,wair,temp,pres,airf=airf,dhum=dhum,dliq=dliq)
    kappa = -g_pp / g_p
    return kappa

def lapserate(wair,temp,pres,airf=None,dhum=None,dliq=None,
    chkvals=False,chktol=_CHKTOL,airf0=None,dhum0=None,dliq0=None,
    chkbnd=False,mathargs=None):
    """Calculate wet air adiabatic lapse rate.
    
    Calculate the adiabatic lapse rate of wet air.
    
    :arg float wair: Total dry air fraction in kg/kg.
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg airf: Dry air fraction in humid air in kg/kg.
    :type airf: float or None
    :arg dhum: Humid air density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dhum: float or None
    :arg dliq: Liquid water density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dliq: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg airf0: Initial guess for the dry fraction in kg/kg. If None
        (default) then `iceair4a._approx_tp` is used.
    :type airf0: float or None
    :arg dhum0: Initial guess for the humid air density in kg/m3. If
        None (default) then `liqair4a._approx_tp` is used.
    :type dhum0: float or None
    :arg dliq0: Initial guess for the liquid water density in kg/m3. If
        None (default) then `liqair4a._approx_tp` is used.
    :type dliq0: float or None
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
    
    >>> lapserate(0.5,300.,1e5)
    1.72449715057e-04
    """
    airf, __, __, dhum, dliq = _eq_atpe(temp=temp,pres=pres,airf=airf,
        dhum=dhum,dliq=dliq,chkvals=chkvals,chktol=chktol,airf0=airf0,
        dhum0=dhum0,dliq0=dliq0,chkbnd=chkbnd,mathargs=mathargs)
    if airf <= wair:
        warnmsg = 'Air with the given parameters is unsaturated'
        warnings.warn(warnmsg,RuntimeWarning)
        gamma = air3b.lapserate(wair,temp,pres,dhum0=dhum0,mathargs=mathargs)
    g_tt = liqair_g(0,2,0,wair,temp,pres,airf=airf,dhum=dhum,dliq=dliq)
    g_tp = liqair_g(0,1,1,wair,temp,pres,airf=airf,dhum=dhum,dliq=dliq)
    gamma = -g_tp / g_tt
    return gamma

def liquidfraction(wair,temp,pres,airf=None,dhum=None,dliq=None,
    chkvals=False,chktol=_CHKTOL,airf0=None,dhum0=None,dliq0=None,
    chkbnd=False,mathargs=None):
    """Calculate wet air liquid water fraction.
    
    Calculate the mass fraction of liquid water in wet air.
    
    :arg float wair: Total dry air fraction in kg/kg.
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg airf: Dry air fraction in humid air in kg/kg.
    :type airf: float or None
    :arg dhum: Humid air density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dhum: float or None
    :arg dliq: Liquid water density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dliq: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg airf0: Initial guess for the dry fraction in kg/kg. If None
        (default) then `iceair4a._approx_tp` is used.
    :type airf0: float or None
    :arg dhum0: Initial guess for the humid air density in kg/m3. If
        None (default) then `liqair4a._approx_tp` is used.
    :type dhum0: float or None
    :arg dliq0: Initial guess for the liquid water density in kg/m3. If
        None (default) then `liqair4a._approx_tp` is used.
    :type dliq0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Liquid water mass fraction in kg/kg.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    :raises RuntimeWarning: If air with the given parameters would be
        unsaturated.
    
    :Examples:
    
    >>> liquidfraction(0.5,300.,1e5)
    0.488546404734
    """
    airf, __, __, dhum, dliq = _eq_atpe(temp=temp,pres=pres,airf=airf,
        dhum=dhum,dliq=dliq,chkvals=chkvals,chktol=chktol,airf0=airf0,
        dhum0=dhum0,dliq0=dliq0,chkbnd=chkbnd,mathargs=mathargs)
    if airf <= wair:
        warnmsg = 'Air with the given parameters is unsaturated'
        warnings.warn(warnmsg,RuntimeWarning)
    wliq = max(1 - wair/airf, 0.)
    return wliq

def vapourfraction(wair,temp,pres,airf=None,dhum=None,dliq=None,
    chkvals=False,chktol=_CHKTOL,airf0=None,dhum0=None,dliq0=None,
    chkbnd=False,mathargs=None):
    """Calculate wet air vapour fraction.
    
    Calculate the mass fraction of water vapour in wet air.
    
    :arg float wair: Total dry air fraction in kg/kg.
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg airf: Dry air fraction in humid air in kg/kg.
    :type airf: float or None
    :arg dhum: Humid air density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dhum: float or None
    :arg dliq: Liquid water density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dliq: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg airf0: Initial guess for the dry fraction in kg/kg. If None
        (default) then `iceair4a._approx_tp` is used.
    :type airf0: float or None
    :arg dhum0: Initial guess for the humid air density in kg/m3. If
        None (default) then `liqair4a._approx_tp` is used.
    :type dhum0: float or None
    :arg dliq0: Initial guess for the liquid water density in kg/m3. If
        None (default) then `liqair4a._approx_tp` is used.
    :type dliq0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Water vapour mass fraction in kg/kg.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    :raises RuntimeWarning: If air with the given parameters would be
        unsaturated.
    
    :Examples:
    
    >>> vapourfraction(0.5,300.,1e5)
    1.14535952655e-2
    """
    airf, __, __, dhum, dliq = _eq_atpe(temp=temp,pres=pres,airf=airf,
        dhum=dhum,dliq=dliq,chkvals=chkvals,chktol=chktol,airf0=airf0,
        dhum0=dhum0,dliq0=dliq0,chkbnd=chkbnd,mathargs=mathargs)
    if airf <= wair:
        warnmsg = 'Air with the given parameters is unsaturated'
        warnings.warn(warnmsg,RuntimeWarning)
    wvap = min(wair * (1-airf)/airf, 1-wair)
    return wvap

