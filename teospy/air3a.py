"""Humid air Gibbs free energy.

This module implements the Gibbs free energy of humid air (dry air with
water vapour) and its derivatives with respect to dry air mass fraction,
temperature, and pressure.

:Examples:

>>> eq_atp(0.9,300.,1e5)
1.09708772444
>>> air_g(0,0,0,0.9,300.,1e5)
4577.93065689
>>> air_g(0,0,1,0.9,300.,1e5)
0.911504137473
>>> air_g(1,1,0,0.9,300.,1e5)
7566.34779196
>>> air_g(0,2,0,0.9,300.,1e5)
-4.15449972148

:Functions:

* eq_atp: Calculate equilibrium quantities at dry fraction, temperature,
    and pressure.
* air_g: Humid air Gibbs free energy with derivatives.

"""

__all__ = ['eq_atp','air_g']

import warnings
import constants0
import air2
import maths3

_RUNIV = constants0.RUNIV
_MDRY = constants0.MDRY
_MWAT = constants0.MWAT
_CHKTOL = constants0.CHKTOL

_chkhumbnds = constants0.chkhumbnds
_air_f = air2.air_f
_eq_pressure = air2.eq_pressure
_newton = maths3.newton


## Auxiliary functions
def _approx_atp(airf,temp,pres):
    """Approximate humid air density from ATP.
    
    Approximate the density of humid air from the dry air mass fraction,
    temperature, and pressure using the ideal gas law.
    
    :arg float airf: Dry air mass fraction in kg/kg.
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :returns: Humid air density in kg/m3.
    """
    raw = _RUNIV * (airf/_MDRY + (1-airf)/_MWAT)
    dhum = pres / (raw*temp)
    return dhum

def _diff_atp(d,airf,temp,pres):
    """Calculate humid air disequilibrium at ATP.
    
    Calculate both sides of the equation
    
        given pressure = pressure of humid air
    
    and their derivative with respect to humid air density. Solving this
    equation gives the equilibrium humid air density for the given dry
    air mass fraction, temperature, and pressure.
    
    :arg float d: Humid air density in kg/m3.
    :arg float airf: Dry air mass fraction in kg/kg.
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :returns: Left-hand side of the equation, right-hand side,
        derivative of LHS, and derivative of RHS.
    :rtype: tuple(float)
    """
    phum = _eq_pressure(0,0,0,airf,temp,d)
    lhs = pres
    rhs = phum
    
    phum_d = _eq_pressure(0,0,1,airf,temp,d)
    dlhs = 0.
    drhs = phum_d
    return lhs, rhs, dlhs, drhs

def eq_atp(airf,temp,pres,dhum=None,chkvals=False,chktol=_CHKTOL,
    dhum0=None,chkbnd=False,mathargs=None):
    """Get primary variables at ATP.
    
    Get the value of the equilibrium humid air density for the given dry
    air mass fraction, temperature, and pressure.
    
    If the calculation has already been done, the results can be passed
    to avoid unnecessary repeat calculations. If enough values are
    passed, they will be checked for consistency if chkvals is True.
    
    :arg float airf: Dry air mass fraction in kg/kg.
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg dhum: Humid air density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dhum: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg dhum0: Initial guess for the humid air density in kg/m3. If
        None (default) then `_approx_atp` is used.
    :type dhum0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Humid air density in kg/m3.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> eq_atp(0.9,300.,1e5)
    1.09708772444
    """
    if dhum is None:
        if dhum0 is None:
            x0 = _approx_atp(airf,temp,pres)
        else:
            x0 = dhum0
        fargs = (airf,temp,pres)
        if mathargs is None:
            mathargs = dict()
        dhum = _newton(_diff_atp,x0,fargs=fargs,**mathargs)
    
    _chkhumbnds(airf,temp,dhum,chkbnd=chkbnd)
    if not chkvals:
        return dhum
    
    lhs, rhs, __, __ = _diff_atp(dhum,airf,temp,pres)
    errs = [abs(lhs/rhs-1)]
    if max(errs) > chktol:
        warnmsg = ('Given value {0} and solution {1} disagree to more than the '
            'tolerance {2}').format(lhs,rhs,chktol)
        warnings.warn(warnmsg,RuntimeWarning)
    return dhum
    

## Public functions
def air_g(drva,drvt,drvp,airf,temp,pres,dhum=None,chkvals=False,
    chktol=_CHKTOL,dhum0=None,chkbnd=False,mathargs=None):
    """Calculate humid air Gibbs free energy with derivatives.
    
    Calculate the specific Gibbs free energy of humid air or its
    derivatives with respect to dry air mass fraction, temperature, and
    pressure.
    
    :arg int drva: Number of dry fraction derivatives.
    :arg int drvt: Number of temperature derivatives.
    :arg int drvp: Number of pressure derivatives.
    :arg float airf: Dry air mass fraction in kg/kg.
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg dhum: Humid air density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dhum: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg dhum0: Initial guess for the humid air density in kg/m3. If
        None (default) then `_approx_atp` is used.
    :type dhum0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Gibbs free energy in units of
            (J/kg) / (kg/kg)^drva / K^drvt / Pa^drvp.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> air_g(0,0,0,0.9,300.,1e5)
    4577.93065689
    >>> air_g(0,0,1,0.9,300.,1e5)
    0.911504137473
    >>> air_g(1,1,0,0.9,300.,1e5)
    7566.34779196
    >>> air_g(0,2,0,0.9,300.,1e5)
    -4.15449972148
    """
    dhum = eq_atp(airf,temp,pres,dhum=dhum,chkvals=chkvals,chktol=chktol,
        dhum0=dhum0,chkbnd=chkbnd,mathargs=mathargs)
    
    # Simple cases: First-order derivatives
    if (drva,drvt,drvp) == (0,0,0):
        f = _air_f(0,0,0,airf,temp,dhum)
        f_d = _air_f(0,0,1,airf,temp,dhum)
        g = f + dhum*f_d
        return g
    elif (drva,drvt,drvp) == (1,0,0):
        f_a = _air_f(1,0,0,airf,temp,dhum)
        g_a = f_a
        return g_a
    elif (drva,drvt,drvp) == (0,1,0):
        f_t = _air_f(0,1,0,airf,temp,dhum)
        g_t = f_t
        return g_t
    elif (drva,drvt,drvp) == (0,0,1):
        g_p = 1. / dhum
        return g_p
    
    # Second-order derivatives require this denominator
    p_d = _eq_pressure(0,0,1,airf,temp,dhum)
    if (drva,drvt,drvp) == (2,0,0):
        f_aa = _air_f(2,0,0,airf,temp,dhum)
        f_ad = _air_f(1,0,1,airf,temp,dhum)
        p_a = _eq_pressure(1,0,0,airf,temp,dhum)
        d_a = -p_a / p_d
        g_aa = f_aa + f_ad*d_a
        return g_aa
    elif (drva,drvt,drvp) == (1,1,0):
        f_at = _air_f(1,1,0,airf,temp,dhum)
        f_ad = _air_f(1,0,1,airf,temp,dhum)
        f_td = _air_f(0,1,1,airf,temp,dhum)
        p_t = _eq_pressure(0,1,0,airf,temp,dhum)
        d_t = -p_t / p_d
        g_at = f_at + f_ad*d_t
        return g_at
    elif (drva,drvt,drvp) == (1,0,1):
        f_ad = _air_f(1,0,1,airf,temp,dhum)
        d_p = p_d**(-1)
        g_ap = f_ad*d_p
        return g_ap
    elif (drva,drvt,drvp) == (0,2,0):
        f_tt = _air_f(0,2,0,airf,temp,dhum)
        f_td = _air_f(0,1,1,airf,temp,dhum)
        p_t = _eq_pressure(0,1,0,airf,temp,dhum)
        d_t = -p_t / p_d
        g_tt = f_tt + f_td * d_t
        return g_tt
    elif (drva,drvt,drvp) == (0,1,1):
        f_td = _air_f(0,1,1,airf,temp,dhum)
        d_p = p_d**(-1)
        g_tp = f_td * d_p
        return g_tp
    elif (drva,drvt,drvp) == (0,0,2):
        d_p = p_d**(-1)
        g_pp = -d_p / dhum**2
        return g_pp
    
    # Should not have made it this far!
    errmsg = 'Derivatives {0} not recognized'.format((drva,drvt,drvp))
    raise ValueError(errmsg)

