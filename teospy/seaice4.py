"""Seawater-ice equilibrium functions.

This module provides thermodynamic functions for seawater in equilibrium
with ice (sea-ice), e.g. the enthalpy of melting. It also provides a
Gibbs free energy function for sea-ice parcels, with primary variables
being the total salinity (mass of salt per mass of salt, liquid, and
ice), temperature, and pressure.

:Examples:

>>> temperature(salt=0.035,pres=1e5)
271.240373585159
>>> enthalpymelt(salt=0.035,pres=1e5)
329942.976285
>>> volumemelt(salt=0.035,pres=1e5)
-9.10140854473e-5
>>> pressure(salt=0.035,temp=270.)
16132047.4385
>>> enthalpymelt(salt=0.035,temp=270.)
326829.393605
>>> volumemelt(salt=0.035,temp=270.)
-9.67135426848e-5
>>> salinity(temp=270.,pres=1e5)
0.0560264150322
>>> enthalpymelt(temp=270.,pres=1e5)
328249.119579
>>> volumemelt(temp=270.,pres=1e5)
-9.18186917900e-5
>>> brinesalinity(270.,1e5)
0.0560264150322
>>> meltingpressure(0.035,270.)
16132047.4385
>>> freezingtemperature(0.035,1e5)
271.240373585
>>> dtfdp(0.035,1e5)
7.48210942879e-8
>>> dtfds(0.035,1e5)
-56.8751336296
>>> sea_ice_g(0,0,0,0.035,270.,1e5)
-414.017574547
>>> sea_ice_g(0,1,0,0.035,270.,1e5)
500.445444181
>>> sea_ice_g(0,1,1,0.035,270.,1e5)
-1.65866446694e-5
>>> brinefraction(0.035,270.,1e5)
0.624705328368
>>> cp(0.035,270.,1e5)
62868.9015126
>>> density(0.035,270.,1e5)
993.156434117
>>> enthalpy(0.035,270.,1e5)
-135534.287503
>>> entropy(0.035,270.,1e5)
-500.445444181
>>> expansion(0.035,270.,1e5)
-1.64731328738
>>> kappa_t(0.035,270.,1e5)
1.56513441348e-9

:Functions:

* :func:`eq_stp`: Calculate primary variables for sea-ice at any two of
  the seawater salinity, temperature, and pressure.
* :func:`densityice`: Sea-ice ice density.
* :func:`densitysea`: Sea-ice seawater density.
* :func:`enthalpyice`: Sea-ice ice enthalpy.
* :func:`enthalpysea`: Sea-ice seawater enthalpy.
* :func:`entropyice`: Sea-ice ice entropy for sea ice.
* :func:`entropysea`: Sea-ice seawater entropy.
* :func:`pressure`: Sea-ice pressure.
* :func:`temperature`: Sea-ice temperature.
* :func:`salinity`: Sea-ice salinity.
* :func:`enthalpymelt`: Enthalpy of melting.
* :func:`volumemelt`: Specific volume of melting.
* :func:`brinesalinity`: Salinity of seawater in equilibrium with ice.
* :func:`meltingpressure`: Pressure of seawater in equilibrium with ice.
* :func:`freezingtemperature`: Temperature of seawater in equilibrium
  with ice.
* :func:`dtfdp`: Freezing point depression of seawater due to pressure.
* :func:`dtfds`: Freezing point depression of seawater due to salinity.
* :func:`eq_seaice`: Calculate primary variables for a sea-ice parcel at
  the given total salinity, temperature, and pressure.
* :func:`seaice_g`: Sea-ice Gibbs free energy with derivatives.
* :func:`brinefraction`: Sea-ice seawater mass fraction.
* :func:`cp`: Sea-ice isobaric heat capacity.
* :func:`density`: Sea-ice total density.
* :func:`enthalpy`: Sea-ice specific enthalpy.
* :func:`entropy`: Sea-ice specific entropy.
* :func:`expansion`: Sea-ice thermal expansion coefficient.
* :func:`kappa_t`: Sea-ice isothermal compressibility.

"""

__all__ = ['eq_stp','densityice','densitysea','enthalpyice','enthalpysea',
    'entropyice','entropysea','pressure','temperature','salinity',
    'enthalpymelt','volumemelt',
    'brinesalinity','meltingpressure','freezingtemperature','dtfdp','dtfds',
    'eq_seaice','seaice_g','brinefraction','cp','density','enthalpy','entropy',
    'expansion','kappa_t']

import warnings
import numpy
import constants0
import flu1
import ice1
import flu2
import ice2
import sal2
import maths3
import flu3a
import sea3a

_CHKTOL = constants0.CHKTOL
_MSAL = constants0.MSAL
_RUNIV = constants0.RUNIV
_RWAT = constants0.RWAT
_TTP = constants0.TTP
_PTPE = constants0.PTPE
_DLTP = constants0.DLTP
_DITP = constants0.DITP
_LILTP = constants0.LILTP
_CLIQ = constants0.CLIQ
_CICE = constants0.CICE
_SAL0 = constants0.SAL0
_RSAL = _RUNIV / _MSAL
_VLTP = _DLTP**(-1)
_VITP = _DITP**(-1)
_chkflubnds = constants0.chkflubnds
_chkicebnds = constants0.chkicebnds
_chksalbnds = constants0.chksalbnds
_flu_f = flu1.flu_f
_ice_g = ice1.ice_g
_eq_pressure = flu2.eq_pressure
_eq_chempot = flu2.eq_chempot
_sal_g = sal2.sal_g
_eq_liqpot = sal2.eq_liqpot
_newton = maths3.newton
_dliq_default = flu3a._dliq_default


## Equilibrium functions
def _approx_st(salt,temp):
    """Approximate PDl at ST.
    
    Approximate the pressure and liquid water density of sea-ice with
    the given salinity and temperature.
    
    :arg float salt: Salinity in kg/kg.
    :arg float temp: Temperature in K.
    :returns: Pressure and liquid water density (both in SI units).
    """
    
    dmu = ((_CLIQ-_CICE)*(temp - _TTP - temp*numpy.log(temp/_TTP))
        + -_LILTP/_TTP*(temp - _TTP) - _RSAL*temp*salt)
    pres = _PTPE + dmu/(_VITP-_VLTP)
    dliq = _dliq_default(temp,pres)
    return pres, dliq

def _approx_sp(salt,pres):
    """Approximate TDl at SP.
    
    Approximate the temperature and liquid water density of sea-ice with
    the given salinity and pressure.
    
    :arg float salt: Salinity in kg/kg.
    :arg float pres: Pressure in Pa.
    :returns: Temperature and liquid water density (both in SI units).
    """
    dt = _RSAL*_TTP*salt + (pres-_PTPE)*(_VITP-_VLTP)
    dt /= _LILTP/_TTP + _RSAL*_SAL0
    temp = _TTP - dt
    dliq = _dliq_default(temp,pres)
    return temp, dliq

def _approx_tp(temp,pres,dliq):
    """Approximate S at TP.
    
    Approximate the salinity of sea-ice with the given temperature and
    pressure.
    
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg float dliq: Liquid water density in kg/m3 (unused).
    :returns: Salinity in kg/kg.
    """
    dmu = ((_CLIQ-_CICE) * (temp-_TTP-temp*numpy.log(temp/_TTP))
        + -_LILTP/_TTP*(temp-_TTP) - (pres-_PTPE)*(_VITP-_VLTP))
    salt = dmu / (_RSAL*temp)
    return salt

def _diff_st(p,dl,salt,temp,useext=False):
    """Calculate sea-ice disequilibrium at ST.
    
    Calculate both sides of the equations
    
        given pressure = pressure of liquid water
        chemical potential of ice = potential of liquid water
    
    and their Jacobians with respect to pressure and liquid water
    density. Solving these equations gives equilibrium values at the
    given salinity and temperature.
    
    :arg float p: Pressure in Pa.
    :arg float dl: Liquid water density in kg/m3.
    :arg float salt: Salinity in kg/kg.
    :arg float temp: Temperature in K.
    :arg bool useext: If False (default) then the salt contribution is
        calculated from _GSCOEFFS; if True, from _GSCOEFFS_EXT.
    :returns: Left-hand side of the equation, right-hand side,
        Jacobian of LHS, and Jacobian of RHS.
    :rtype: tuple(array(float))
    """
    pl = _eq_pressure(0,0,temp,dl)
    gi = _ice_g(0,0,temp,p)
    gl = _eq_chempot(0,0,temp,dl)
    gl += _eq_liqpot(0,0,0,salt,temp,p,useext=useext)
    lhs = numpy.array([p, gi])
    rhs = numpy.array([pl, gl])
    
    pl_d = _eq_pressure(0,1,temp,dl)
    gi_p = _ice_g(0,1,temp,p)
    gl_d = _eq_chempot(0,1,temp,dl)
    gl_p = _eq_liqpot(0,0,1,salt,temp,p,useext=useext)
    dlhs = numpy.array([[1.,0.], [gi_p,0.]])
    drhs = numpy.array([[0.,pl_d], [gl_p,gl_d]])
    return lhs, rhs, dlhs, drhs

def _diff_sp(t,dl,salt,pres,useext=False):
    """Calculate sea-ice disequilibrium at SP.
    
    Calculate both sides of the equations
    
        given pressure = pressure of liquid water
        chemical potential of ice = potential of liquid water
    
    and their Jacobians with respect to temperature and liquid water
    density. Solving these equations gives equilibrium values at the
    given salinity and pressure.

    :arg float t: Temperature in K.
    :arg float dl: Liquid water density in kg/m3.
    :arg float salt: Salinity in kg/kg.
    :arg float pres: Pressure in Pa.
    :arg bool useext: If False (default) then the salt contribution is
        calculated from _GSCOEFFS; if True, from _GSCOEFFS_EXT.
    :returns: Left-hand side of the equation, right-hand side,
        Jacobian of LHS, and Jacobian of RHS.
    :rtype: tuple(array(float))
    """
    pl = _eq_pressure(0,0,t,dl)
    gi = _ice_g(0,0,t,pres)
    gl = _eq_chempot(0,0,t,dl)
    gl += _eq_liqpot(0,0,0,salt,t,pres,useext=useext)
    lhs = numpy.array([pres, gi])
    rhs = numpy.array([pl, gl])
    
    pl_t = _eq_pressure(1,0,t,dl)
    pl_d = _eq_pressure(0,1,t,dl)
    gi_t = _ice_g(1,0,t,pres)
    gl_t = _eq_chempot(1,0,t,dl)
    gl_t += _eq_liqpot(0,1,0,salt,t,pres,useext=useext)
    gl_d = _eq_chempot(0,1,t,dl)
    dlhs = numpy.array([[0.,0.], [gi_t,0.]])
    drhs = numpy.array([[pl_t,pl_d], [gl_t,gl_d]])
    return lhs, rhs, dlhs, drhs

def _diff_tp(s,temp,pres,dliq,useext=False):
    """Calculate sea-ice disequilibrium at TP.
    
    Calculate both sides of the equation
    
        chemical potential of ice = potential of liquid water
    
    and their derivatives with respect to salinity. Solving these
    equations gives the equilibrium salinity at the given temperature
    and pressure.

    :arg float s: Salinity in kg/kg.
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg float dliq: Liquid water density in kg/m3.
    :arg bool useext: If False (default) then the salt contribution is
        calculated from _GSCOEFFS; if True, from _GSCOEFFS_EXT.
    :returns: Left-hand side of the equation, right-hand side,
        derivative of LHS, and derivative of RHS.
    :rtype: tuple(float)
    """
    gi = _ice_g(0,0,temp,pres)
    gl = _eq_chempot(0,0,temp,dliq)
    gl += _eq_liqpot(0,0,0,s,temp,pres,useext=useext)
    lhs = gi
    rhs = gl
    
    gl_s = _eq_liqpot(1,0,0,s,temp,pres,useext=useext)
    dlhs = 0.
    drhs = gl_s
    return lhs, rhs, dlhs, drhs

def eq_stp(salt=None,temp=None,pres=None,dliq=None,chkvals=False,
    chktol=_CHKTOL,salt0=None,temp0=None,pres0=None,dliq0=None,chkbnd=False,
    useext=False,mathargs=None):
    """Get primary sea-ice variables at STP.
    
    Get the values of all primary variables for sea-ice in equilibrium.
    At least two of the salinity, temperature, and pressure must be
    provided.
    
    If the calculation has already been done, the results can be passed
    to avoid unnecessary repeat calculations. If enough values are
    passed, they will be checked for consistency if chkvals is True.
    
    :arg salt: Salinity in kg/kg. If unknown, pass None (default) and it
        will be calculated.
    :type salt: float or None
    :arg temp: Temperature in K. If unknown, pass None (default) and it
        will be calculated.
    :type temp: float or None
    :arg pres: Pressure in Pa. If unknown, pass None (default) and it
        will be calculated.
    :type pres: float or None
    :arg dliq: Density of liquid water in seawater in kg/m3. If
        unknown, pass None (default) and it will be calculated.
    :type dliq: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg salt0: Initial guess for the salinity in kg/kg. If None
        (default) then `_approx_tp` is used.
    :type salt0: float or None
    :arg temp0: Initial guess for the temperature in K. If None
        (default) then `_approx_sp` is used.
    :type temp0: float or None
    :arg pres0: Initial guess for the pressure in Pa. If None (default)
        then `_approx_st` is used.
    :type pres0: float or None
    :arg dliq0: Initial guess for the liquid water density in kg/m3. If
        None (default) then `flu3a._dliq_default` is used.
    :type dliq0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg bool useext: If False (default) then the salt contribution is
        calculated from _GSCOEFFS; if True, from _GSCOEFFS_EXT.
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Salinity, temperature, pressure, and seawater liquid
        density (all in SI units).
    :raises ValueError: If fewer than two of salt, temp, and pres are
        provided.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    """
    if sum(val is None for val in (salt,temp,pres)) > 1:
        errmsg = 'Must provide at least two of (salt,temp,pres)'
        raise ValueError(errmsg)
    if mathargs is None:
        mathargs = dict()
    fkwargs = {'useext': useext}
    
    if salt is None:
        dliq = flu3a.eq_tp_liq(temp,pres,dliq=dliq,dliq0=dliq0,
            mathargs=mathargs)
        fargs = (temp,pres,dliq)
        salt = _newton(_diff_tp,salt0,_approx_tp,fargs=fargs,fkwargs=fkwargs,
            **mathargs)
    elif temp is None:
        x0 = (temp0,dliq0)
        fargs = (salt,pres)
        x1 = _newton(_diff_sp,x0,_approx_sp,fargs=fargs,fkwargs=fkwargs,
            **mathargs)
        temp, dliq = x1
    elif pres is None:
        x0 = (pres0,dliq0)
        fargs = (salt,temp)
        x1 = _newton(_diff_st,x0,_approx_st,fargs=fargs,fkwargs=fkwargs,
            **mathargs)
        pres, dliq = x1
    elif dliq is None:
        dliq = flu3a.eq_tp_liq(temp,pres,dliq0=dliq0,mathargs=mathargs)
    
    _chkflubnds(temp,dliq,chkbnd=chkbnd)
    _chkicebnds(temp,pres,chkbnd=chkbnd)
    _chksalbnds(salt,temp,pres,chkbnd=chkbnd)
    if not chkvals:
        return salt, temp, pres, dliq
    
    lhs, rhs, __, __ = _diff_st(pres,dliq,salt,temp,useext=useext)
    errs = list()
    for (l,r) in zip(lhs,rhs):
        if abs(r) >= chktol:
            errs.append(abs(l/r-1))
        else:
            errs.append(abs(l-r))
    if max(errs) > chktol:
        warnmsg = ('Given values {0} and solutions {1} disagree to more than '
            'the tolerance {2}').format(lhs,rhs,chktol)
        warnings.warn(warnmsg,RuntimeWarning)
    return salt, temp, pres, dliq


## Thermodynamic functions
def densityice(salt=None,temp=None,pres=None,dliq=None,chkvals=False,
    chktol=_CHKTOL,salt0=None,temp0=None,pres0=None,dliq0=None,
    chkbnd=False,useext=False,mathargs=None):
    """Calculate sea-ice ice density.
    
    Calculate the density of ice in sea-ice.
    
    :arg salt: Salinity in kg/kg. If unknown, pass None (default) and it
        will be calculated.
    :type salt: float or None
    :arg temp: Temperature in K. If unknown, pass None (default) and it
        will be calculated.
    :type temp: float or None
    :arg pres: Pressure in Pa. If unknown, pass None (default) and it
        will be calculated.
    :type pres: float or None
    :arg dliq: Density of liquid water in seawater in kg/m3. If
        unknown, pass None (default) and it will be calculated.
    :type dliq: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg salt0: Initial guess for the salinity in kg/kg. If None
        (default) then `_approx_tp` is used.
    :type salt0: float or None
    :arg temp0: Initial guess for the temperature in K. If None
        (default) then `_approx_sp` is used.
    :type temp0: float or None
    :arg pres0: Initial guess for the pressure in Pa. If None (default)
        then `_approx_st` is used.
    :type pres0: float or None
    :arg dliq0: Initial guess for the liquid water density in kg/m3. If
        None (default) then `flu3a._dliq_default` is used.
    :type dliq0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg bool useext: If False (default) then the salt contribution is
        calculated from _GSCOEFFS; if True, from _GSCOEFFS_EXT.
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Density in kg/m3.
    :raises ValueError: If fewer than two of salt, temp, and pres are
        provided.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> densityice(salt=0.035,pres=1e5)
    917.000739687
    >>> densityice(salt=0.035,temp=270.)
    918.898527655
    >>> densityice(temp=270.,pres=1e5)
    917.181167192
    """
    salt, temp, pres, dliq = eq_stp(salt=salt,temp=temp,pres=pres,dliq=dliq,
        chkvals=chkvals,chktol=chktol,salt0=salt0,temp0=temp0,pres0=pres0,
        dliq0=dliq0,chkbnd=chkbnd,useext=useext,mathargs=mathargs)
    dice = ice2.density(temp,pres)
    return dice

def densitysea(salt=None,temp=None,pres=None,dliq=None,chkvals=False,
    chktol=_CHKTOL,salt0=None,temp0=None,pres0=None,dliq0=None,
    chkbnd=False,useext=False,mathargs=None):
    """Calculate sea-ice seawater density.
    
    Calculate the density of seawater in sea-ice.
    
    :arg salt: Salinity in kg/kg. If unknown, pass None (default) and it
        will be calculated.
    :type salt: float or None
    :arg temp: Temperature in K. If unknown, pass None (default) and it
        will be calculated.
    :type temp: float or None
    :arg pres: Pressure in Pa. If unknown, pass None (default) and it
        will be calculated.
    :type pres: float or None
    :arg dliq: Density of liquid water in seawater in kg/m3. If
        unknown, pass None (default) and it will be calculated.
    :type dliq: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg salt0: Initial guess for the salinity in kg/kg. If None
        (default) then `_approx_tp` is used.
    :type salt0: float or None
    :arg temp0: Initial guess for the temperature in K. If None
        (default) then `_approx_sp` is used.
    :type temp0: float or None
    :arg pres0: Initial guess for the pressure in Pa. If None (default)
        then `_approx_st` is used.
    :type pres0: float or None
    :arg dliq0: Initial guess for the liquid water density in kg/m3. If
        None (default) then `flu3a._dliq_default` is used.
    :type dliq0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg bool useext: If False (default) then the salt contribution is
        calculated from _GSCOEFFS; if True, from _GSCOEFFS_EXT.
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Density in kg/m3.
    :raises ValueError: If fewer than two of salt, temp, and pres are
        provided.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> densitysea(salt=0.035,pres=1e5)
    1028.05199645
    >>> densitysea(salt=0.035,temp=270.)
    1035.73670169
    >>> densitysea(temp=270.,pres=1e5)
    1045.16805918
    """
    salt, temp, pres, dliq = eq_stp(salt=salt,temp=temp,pres=pres,dliq=dliq,
        chkvals=chkvals,chktol=chktol,salt0=salt0,temp0=temp0,pres0=pres0,
        dliq0=dliq0,chkbnd=chkbnd,useext=useext,mathargs=mathargs)
    dsea = sea3a.density(salt,temp,pres,dliq=dliq,useext=useext)
    return dsea

def enthalpyice(salt=None,temp=None,pres=None,dliq=None,chkvals=False,
    chktol=_CHKTOL,salt0=None,temp0=None,pres0=None,dliq0=None,
    chkbnd=False,useext=False,mathargs=None):
    """Calculate sea-ice ice enthalpy.
    
    Calculate the specific enthalpy of ice in sea-ice.
    
    :arg salt: Salinity in kg/kg. If unknown, pass None (default) and it
        will be calculated.
    :type salt: float or None
    :arg temp: Temperature in K. If unknown, pass None (default) and it
        will be calculated.
    :type temp: float or None
    :arg pres: Pressure in Pa. If unknown, pass None (default) and it
        will be calculated.
    :type pres: float or None
    :arg dliq: Density of liquid water in seawater in kg/m3. If
        unknown, pass None (default) and it will be calculated.
    :type dliq: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg salt0: Initial guess for the salinity in kg/kg. If None
        (default) then `_approx_tp` is used.
    :type salt0: float or None
    :arg temp0: Initial guess for the temperature in K. If None
        (default) then `_approx_sp` is used.
    :type temp0: float or None
    :arg pres0: Initial guess for the pressure in Pa. If None (default)
        then `_approx_st` is used.
    :type pres0: float or None
    :arg dliq0: Initial guess for the liquid water density in kg/m3. If
        None (default) then `flu3a._dliq_default` is used.
    :type dliq0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg bool useext: If False (default) then the salt contribution is
        calculated from _GSCOEFFS; if True, from _GSCOEFFS_EXT.
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Enthalpy in J/kg.
    :raises ValueError: If fewer than two of salt, temp, and pres are
        provided.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> enthalpyice(salt=0.035,pres=1e5)
    -337351.999358
    >>> enthalpyice(salt=0.035,temp=270.)
    -323205.968289
    >>> enthalpyice(temp=270.,pres=1e5)
    -339929.555499
    """
    salt, temp, pres, dliq = eq_stp(salt=salt,temp=temp,pres=pres,dliq=dliq,
        chkvals=chkvals,chktol=chktol,salt0=salt0,temp0=temp0,pres0=pres0,
        dliq0=dliq0,chkbnd=chkbnd,useext=useext,mathargs=mathargs)
    hice = ice2.enthalpy(temp,pres)
    return hice

def enthalpysea(salt=None,temp=None,pres=None,dliq=None,chkvals=False,
    chktol=_CHKTOL,salt0=None,temp0=None,pres0=None,dliq0=None,
    chkbnd=False,useext=False,mathargs=None):
    """Calculate sea-ice seawater enthalpy.
    
    Calculate the specific enthalpy of seawater in sea-ice.
    
    :arg salt: Salinity in kg/kg. If unknown, pass None (default) and it
        will be calculated.
    :type salt: float or None
    :arg temp: Temperature in K. If unknown, pass None (default) and it
        will be calculated.
    :type temp: float or None
    :arg pres: Pressure in Pa. If unknown, pass None (default) and it
        will be calculated.
    :type pres: float or None
    :arg dliq: Density of liquid water in seawater in kg/m3. If
        unknown, pass None (default) and it will be calculated.
    :type dliq: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg salt0: Initial guess for the salinity in kg/kg. If None
        (default) then `_approx_tp` is used.
    :type salt0: float or None
    :arg temp0: Initial guess for the temperature in K. If None
        (default) then `_approx_sp` is used.
    :type temp0: float or None
    :arg pres0: Initial guess for the pressure in Pa. If None (default)
        then `_approx_st` is used.
    :type pres0: float or None
    :arg dliq0: Initial guess for the liquid water density in kg/m3. If
        None (default) then `flu3a._dliq_default` is used.
    :type dliq0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg bool useext: If False (default) then the salt contribution is
        calculated from _GSCOEFFS; if True, from _GSCOEFFS_EXT.
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Enthalpy in J/kg.
    :raises ValueError: If fewer than two of salt, temp, and pres are
        provided.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> enthalpysea(salt=0.035,pres=1e5)
    -7613.19337919
    >>> enthalpysea(salt=0.035,temp=270.)
    2832.94910407
    >>> enthalpysea(temp=270.,pres=1e5)
    -12742.8664892
    """
    salt, temp, pres, dliq = eq_stp(salt=salt,temp=temp,pres=pres,dliq=dliq,
        chkvals=chkvals,chktol=chktol,salt0=salt0,temp0=temp0,pres0=pres0,
        dliq0=dliq0,chkbnd=chkbnd,useext=useext,mathargs=mathargs)
    hsea = sea3a.enthalpy(salt,temp,pres,dliq=dliq,useext=useext)
    return hsea

def entropyice(salt=None,temp=None,pres=None,dliq=None,chkvals=False,
    chktol=_CHKTOL,salt0=None,temp0=None,pres0=None,dliq0=None,
    chkbnd=False,useext=False,mathargs=None):
    """Calculate sea-ice ice entropy.
    
    Calculate the specific entropy of ice in sea-ice.
    
    :arg salt: Salinity in kg/kg. If unknown, pass None (default) and it
        will be calculated.
    :type salt: float or None
    :arg temp: Temperature in K. If unknown, pass None (default) and it
        will be calculated.
    :type temp: float or None
    :arg pres: Pressure in Pa. If unknown, pass None (default) and it
        will be calculated.
    :type pres: float or None
    :arg dliq: Density of liquid water in seawater in kg/m3. If
        unknown, pass None (default) and it will be calculated.
    :type dliq: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg salt0: Initial guess for the salinity in kg/kg. If None
        (default) then `_approx_tp` is used.
    :type salt0: float or None
    :arg temp0: Initial guess for the temperature in K. If None
        (default) then `_approx_sp` is used.
    :type temp0: float or None
    :arg pres0: Initial guess for the pressure in Pa. If None (default)
        then `_approx_st` is used.
    :type pres0: float or None
    :arg dliq0: Initial guess for the liquid water density in kg/m3. If
        None (default) then `flu3a._dliq_default` is used.
    :type dliq0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg bool useext: If False (default) then the salt contribution is
        calculated from _GSCOEFFS; if True, from _GSCOEFFS_EXT.
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Entropy in J/kg/K.
    :raises ValueError: If fewer than two of salt, temp, and pres are
        provided.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> entropyice(salt=0.035,pres=1e5)
    -1235.44872812
    >>> entropyice(salt=0.035,temp=270.)
    -1247.71314646
    >>> entropyice(temp=270.,pres=1e5)
    -1244.97335506
    """
    salt, temp, pres, dliq = eq_stp(salt=salt,temp=temp,pres=pres,dliq=dliq,
        chkvals=chkvals,chktol=chktol,salt0=salt0,temp0=temp0,pres0=pres0,
        dliq0=dliq0,chkbnd=chkbnd,useext=useext,mathargs=mathargs)
    sice = ice2.entropy(temp,pres)
    return sice

def entropysea(salt=None,temp=None,pres=None,dliq=None,chkvals=False,
    chktol=_CHKTOL,salt0=None,temp0=None,pres0=None,dliq0=None,
    chkbnd=False,useext=False,mathargs=None):
    """Calculate sea-ice seawater entropy.
    
    Calculate the specific entropy of seawater in sea-ice.
    
    :arg salt: Salinity in kg/kg. If unknown, pass None (default) and it
        will be calculated.
    :type salt: float or None
    :arg temp: Temperature in K. If unknown, pass None (default) and it
        will be calculated.
    :type temp: float or None
    :arg pres: Pressure in Pa. If unknown, pass None (default) and it
        will be calculated.
    :type pres: float or None
    :arg dliq: Density of liquid water in seawater in kg/m3. If
        unknown, pass None (default) and it will be calculated.
    :type dliq: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg salt0: Initial guess for the salinity in kg/kg. If None
        (default) then `_approx_tp` is used.
    :type salt0: float or None
    :arg temp0: Initial guess for the temperature in K. If None
        (default) then `_approx_sp` is used.
    :type temp0: float or None
    :arg pres0: Initial guess for the pressure in Pa. If None (default)
        then `_approx_st` is used.
    :type pres0: float or None
    :arg dliq0: Initial guess for the liquid water density in kg/m3. If
        None (default) then `flu3a._dliq_default` is used.
    :type dliq0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg bool useext: If False (default) then the salt contribution is
        calculated from _GSCOEFFS; if True, from _GSCOEFFS_EXT.
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Entropy in J/kg/K.
    :raises ValueError: If fewer than two of salt, temp, and pres are
        provided.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> entropysea(salt=0.035,pres=1e5)
    -27.9264598103
    >>> entropysea(salt=0.035,temp=270.)
    -46.7361169560
    >>> entropysea(temp=270.,pres=1e5)
    -53.1667911144
    """
    salt, temp, pres, dliq = eq_stp(salt=salt,temp=temp,pres=pres,dliq=dliq,
        chkvals=chkvals,chktol=chktol,salt0=salt0,temp0=temp0,pres0=pres0,
        dliq0=dliq0,chkbnd=chkbnd,useext=useext,mathargs=mathargs)
    ssea = sea3a.entropy(salt,temp,pres,dliq=dliq,useext=useext)
    return ssea

def pressure(salt=None,temp=None,pres=None,dliq=None,chkvals=False,
    chktol=_CHKTOL,salt0=None,temp0=None,pres0=None,dliq0=None,
    chkbnd=False,useext=False,mathargs=None):
    """Calculate sea-ice pressure.
    
    Calculate the pressure of sea-ice.
    
    :arg salt: Salinity in kg/kg. If unknown, pass None (default) and it
        will be calculated.
    :type salt: float or None
    :arg temp: Temperature in K. If unknown, pass None (default) and it
        will be calculated.
    :type temp: float or None
    :arg pres: Pressure in Pa. If unknown, pass None (default) and it
        will be calculated.
    :type pres: float or None
    :arg dliq: Density of liquid water in seawater in kg/m3. If
        unknown, pass None (default) and it will be calculated.
    :type dliq: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg salt0: Initial guess for the salinity in kg/kg. If None
        (default) then `_approx_tp` is used.
    :type salt0: float or None
    :arg temp0: Initial guess for the temperature in K. If None
        (default) then `_approx_sp` is used.
    :type temp0: float or None
    :arg pres0: Initial guess for the pressure in Pa. If None (default)
        then `_approx_st` is used.
    :type pres0: float or None
    :arg dliq0: Initial guess for the liquid water density in kg/m3. If
        None (default) then `flu3a._dliq_default` is used.
    :type dliq0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg bool useext: If False (default) then the salt contribution is
        calculated from _GSCOEFFS; if True, from _GSCOEFFS_EXT.
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Pressure in Pa.
    :raises ValueError: If fewer than two of salt, temp, and pres are
        provided.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> pressure(salt=0.035,temp=270.)
    16132047.4385
    """
    salt, temp, pres, dliq = eq_stp(salt=salt,temp=temp,pres=pres,dliq=dliq,
        chkvals=chkvals,chktol=chktol,salt0=salt0,temp0=temp0,pres0=pres0,
        dliq0=dliq0,chkbnd=chkbnd,useext=useext,mathargs=mathargs)
    return pres

def temperature(salt=None,temp=None,pres=None,dliq=None,chkvals=False,
    chktol=_CHKTOL,salt0=None,temp0=None,pres0=None,dliq0=None,
    chkbnd=False,useext=False,mathargs=None):
    """Calculate sea-ice temperature.
    
    Calculate the temperature of sea-ice.
    
    :arg salt: Salinity in kg/kg. If unknown, pass None (default) and it
        will be calculated.
    :type salt: float or None
    :arg temp: Temperature in K. If unknown, pass None (default) and it
        will be calculated.
    :type temp: float or None
    :arg pres: Pressure in Pa. If unknown, pass None (default) and it
        will be calculated.
    :type pres: float or None
    :arg dliq: Density of liquid water in seawater in kg/m3. If
        unknown, pass None (default) and it will be calculated.
    :type dliq: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg salt0: Initial guess for the salinity in kg/kg. If None
        (default) then `_approx_tp` is used.
    :type salt0: float or None
    :arg temp0: Initial guess for the temperature in K. If None
        (default) then `_approx_sp` is used.
    :type temp0: float or None
    :arg pres0: Initial guess for the pressure in Pa. If None (default)
        then `_approx_st` is used.
    :type pres0: float or None
    :arg dliq0: Initial guess for the liquid water density in kg/m3. If
        None (default) then `flu3a._dliq_default` is used.
    :type dliq0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg bool useext: If False (default) then the salt contribution is
        calculated from _GSCOEFFS; if True, from _GSCOEFFS_EXT.
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Temperature in K.
    :raises ValueError: If fewer than two of salt, temp, and pres are
        provided.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> temperature(salt=0.035,pres=1e5)
    271.240373585159
    """
    salt, temp, pres, dliq = eq_stp(salt=salt,temp=temp,pres=pres,dliq=dliq,
        chkvals=chkvals,chktol=chktol,salt0=salt0,temp0=temp0,pres0=pres0,
        dliq0=dliq0,chkbnd=chkbnd,useext=useext,mathargs=mathargs)
    return temp

def salinity(salt=None,temp=None,pres=None,dliq=None,chkvals=False,
    chktol=_CHKTOL,salt0=None,temp0=None,pres0=None,dliq0=None,
    chkbnd=False,useext=False,mathargs=None):
    """Calculate sea-ice salinity.
    
    Calculate the salinity of sea-ice.
    
    :arg salt: Salinity in kg/kg. If unknown, pass None (default) and it
        will be calculated.
    :type salt: float or None
    :arg temp: Temperature in K. If unknown, pass None (default) and it
        will be calculated.
    :type temp: float or None
    :arg pres: Pressure in Pa. If unknown, pass None (default) and it
        will be calculated.
    :type pres: float or None
    :arg dliq: Density of liquid water in seawater in kg/m3. If
        unknown, pass None (default) and it will be calculated.
    :type dliq: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg salt0: Initial guess for the salinity in kg/kg. If None
        (default) then `_approx_tp` is used.
    :type salt0: float or None
    :arg temp0: Initial guess for the temperature in K. If None
        (default) then `_approx_sp` is used.
    :type temp0: float or None
    :arg pres0: Initial guess for the pressure in Pa. If None (default)
        then `_approx_st` is used.
    :type pres0: float or None
    :arg dliq0: Initial guess for the liquid water density in kg/m3. If
        None (default) then `flu3a._dliq_default` is used.
    :type dliq0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg bool useext: If False (default) then the salt contribution is
        calculated from _GSCOEFFS; if True, from _GSCOEFFS_EXT.
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Salinity in kg/kg.
    :raises ValueError: If fewer than two of salt, temp, and pres are
        provided.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> salinity(temp=270.,pres=1e5)
    0.0560264150322
    """
    salt, temp, pres, dliq = eq_stp(salt=salt,temp=temp,pres=pres,dliq=dliq,
        chkvals=chkvals,chktol=chktol,salt0=salt0,temp0=temp0,pres0=pres0,
        dliq0=dliq0,chkbnd=chkbnd,useext=useext,mathargs=mathargs)
    return salt

def enthalpymelt(salt=None,temp=None,pres=None,dliq=None,chkvals=False,
    chktol=_CHKTOL,salt0=None,temp0=None,pres0=None,dliq0=None,
    chkbnd=False,useext=False,mathargs=None):
    """Calculate the enthalpy of melting.
    
    Calculate the specific enthalpy of melting of sea-ice.
    
    :arg salt: Salinity in kg/kg. If unknown, pass None (default) and it
        will be calculated.
    :type salt: float or None
    :arg temp: Temperature in K. If unknown, pass None (default) and it
        will be calculated.
    :type temp: float or None
    :arg pres: Pressure in Pa. If unknown, pass None (default) and it
        will be calculated.
    :type pres: float or None
    :arg dliq: Density of liquid water in seawater in kg/m3. If
        unknown, pass None (default) and it will be calculated.
    :type dliq: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg salt0: Initial guess for the salinity in kg/kg. If None
        (default) then `_approx_tp` is used.
    :type salt0: float or None
    :arg temp0: Initial guess for the temperature in K. If None
        (default) then `_approx_sp` is used.
    :type temp0: float or None
    :arg pres0: Initial guess for the pressure in Pa. If None (default)
        then `_approx_st` is used.
    :type pres0: float or None
    :arg dliq0: Initial guess for the liquid water density in kg/m3. If
        None (default) then `flu3a._dliq_default` is used.
    :type dliq0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg bool useext: If False (default) then the salt contribution is
        calculated from _GSCOEFFS; if True, from _GSCOEFFS_EXT.
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Enthalpy in J/kg.
    :raises ValueError: If fewer than two of salt, temp, and pres are
        provided.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> enthalpymelt(salt=0.035,pres=1e5)
    329942.976285
    >>> enthalpymelt(salt=0.035,temp=270.)
    326829.393605
    >>> enthalpymelt(temp=270.,pres=1e5)
    328249.119579
    """
    salt, temp, pres, dliq = eq_stp(salt=salt,temp=temp,pres=pres,dliq=dliq,
        chkvals=chkvals,chktol=chktol,salt0=salt0,temp0=temp0,pres0=pres0,
        dliq0=dliq0,chkbnd=chkbnd,useext=useext,mathargs=mathargs)
    fl_t = _flu_f(1,0,temp,dliq)
    gs_t = _sal_g(0,1,0,salt,temp,pres,useext=useext)
    gs_st = _sal_g(1,1,0,salt,temp,pres,useext=useext)
    gi_t = _ice_g(1,0,temp,pres)
    hmelt = temp * (gi_t - (fl_t + gs_t - salt*gs_st))
    return hmelt

def volumemelt(salt=None,temp=None,pres=None,dliq=None,chkvals=False,
    chktol=_CHKTOL,salt0=None,temp0=None,pres0=None,dliq0=None,
    chkbnd=False,useext=False,mathargs=None):
    """Calculate the volume of melting.
    
    Calculate the specific volume of melting of sea-ice.
    
    :arg salt: Salinity in kg/kg. If unknown, pass None (default) and it
        will be calculated.
    :type salt: float or None
    :arg temp: Temperature in K. If unknown, pass None (default) and it
        will be calculated.
    :type temp: float or None
    :arg pres: Pressure in Pa. If unknown, pass None (default) and it
        will be calculated.
    :type pres: float or None
    :arg dliq: Density of liquid water in seawater in kg/m3. If
        unknown, pass None (default) and it will be calculated.
    :type dliq: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg salt0: Initial guess for the salinity in kg/kg. If None
        (default) then `_approx_tp` is used.
    :type salt0: float or None
    :arg temp0: Initial guess for the temperature in K. If None
        (default) then `_approx_sp` is used.
    :type temp0: float or None
    :arg pres0: Initial guess for the pressure in Pa. If None (default)
        then `_approx_st` is used.
    :type pres0: float or None
    :arg dliq0: Initial guess for the liquid water density in kg/m3. If
        None (default) then `flu3a._dliq_default` is used.
    :type dliq0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg bool useext: If False (default) then the salt contribution is
        calculated from _GSCOEFFS; if True, from _GSCOEFFS_EXT.
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Specific volume in m3/kg.
    :raises ValueError: If fewer than two of salt, temp, and pres are
        provided.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> volumemelt(salt=0.035,pres=1e5)
    -9.10140854473e-5
    >>> volumemelt(salt=0.035,temp=270.)
    -9.67135426848e-5
    >>> volumemelt(temp=270.,pres=1e5)
    -9.18186917900e-5
    """
    salt, temp, pres, dliq = eq_stp(salt=salt,temp=temp,pres=pres,dliq=dliq,
        chkvals=chkvals,chktol=chktol,salt0=salt0,temp0=temp0,pres0=pres0,
        dliq0=dliq0,chkbnd=chkbnd,useext=useext,mathargs=mathargs)
    gs_p = _sal_g(0,0,1,salt,temp,pres,useext=useext)
    gs_sp = _sal_g(1,0,1,salt,temp,pres,useext=useext)
    gi_p = _ice_g(0,1,temp,pres)
    vmelt = dliq**(-1) + gs_p - salt*gs_sp - gi_p
    return vmelt


## Thermodynamic functions of two variables
def brinesalinity(temp,pres,salt=None,dliq=None,chkvals=False,
    chktol=_CHKTOL,salt0=None,temp0=None,pres0=None,dliq0=None,
    chkbnd=False,useext=False,mathargs=None):
    """Calculate sea-ice brine salinity.
    
    Calculate the salinity of seawater (brine) in equilibrium with ice
    of the given temperature and pressure.
    
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg salt: Salinity in kg/kg. If unknown, pass None (default) and it
        will be calculated.
    :type salt: float or None
    :arg dliq: Density of liquid water in seawater in kg/m3. If
        unknown, pass None (default) and it will be calculated.
    :type dliq: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg salt0: Initial guess for the salinity in kg/kg. If None
        (default) then `_approx_tp` is used.
    :type salt0: float or None
    :arg dliq0: Initial guess for the liquid water density in kg/m3. If
        None (default) then `flu3a._dliq_default` is used.
    :type dliq0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg bool useext: If False (default) then the salt contribution is
        calculated from _GSCOEFFS; if True, from _GSCOEFFS_EXT.
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Salinity in kg/kg.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> brinesalinity(270.,1e5)
    0.0560264150322
    """
    salt, __, __, dliq = eq_stp(temp=temp,pres=pres,salt=salt,dliq=dliq,
        chkvals=chkvals,chktol=chktol,salt0=salt0,dliq0=dliq0,chkbnd=chkbnd,
        useext=useext,mathargs=mathargs)
    return salt

def meltingpressure(salt,temp,pres=None,dliq=None,chkvals=False,
    chktol=_CHKTOL,pres0=None,dliq0=None,chkbnd=False,useext=False,
    mathargs=None):
    """Calculate sea-ice melting pressure.
    
    Calculate the pressure required to melt ice into seawater at the
    given salinity and temperature.
    
    :arg float salt: Salinity in kg/kg.
    :arg float temp: Temperature in K.
    :arg pres: Pressure in Pa. If unknown, pass None (default) and it
        will be calculated.
    :type pres: float or None
    :arg dliq: Density of liquid water in seawater in kg/m3. If
        unknown, pass None (default) and it will be calculated.
    :type dliq: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg pres0: Initial guess for the pressure in Pa. If None (default)
        then `_approx_st` is used.
    :type pres0: float or None
    :arg dliq0: Initial guess for the liquid water density in kg/m3. If
        None (default) then `flu3a._dliq_default` is used.
    :type dliq0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg bool useext: If False (default) then the salt contribution is
        calculated from _GSCOEFFS; if True, from _GSCOEFFS_EXT.
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Pressure in Pa.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> meltingpressure(0.035,270.)
    16132047.4385
    """
    __, __, pres, dliq = eq_stp(temp=temp,pres=pres,salt=salt,dliq=dliq,
        chkvals=chkvals,chktol=chktol,pres0=pres0,dliq0=dliq0,chkbnd=chkbnd,
        useext=useext,mathargs=mathargs)
    return pres

def freezingtemperature(salt,pres,temp=None,dliq=None,chkvals=False,
    chktol=_CHKTOL,temp0=None,dliq0=None,chkbnd=False,useext=False,
    mathargs=None):
    """Calculate sea-ice freezing temperature.
    
    Calculate the temperature required to freeze seawater at the given
    salinity and pressure.
    
    :arg float salt: Salinity in kg/kg.
    :arg float pres: Pressure in Pa.
    :arg temp: Temperature in K. If unknown, pass None (default) and it
        will be calculated.
    :type temp: float or None
    :arg dliq: Density of liquid water in seawater in kg/m3. If
        unknown, pass None (default) and it will be calculated.
    :type dliq: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg temp0: Initial guess for the temperature in K. If None
        (default) then `_approx_sp` is used.
    :type temp0: float or None
    :arg dliq0: Initial guess for the liquid water density in kg/m3. If
        None (default) then `flu3a._dliq_default` is used.
    :type dliq0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg bool useext: If False (default) then the salt contribution is
        calculated from _GSCOEFFS; if True, from _GSCOEFFS_EXT.
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Temperature in K.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> freezingtemperature(0.035,1e5)
    271.240373585
    """
    __, temp, __, dliq = eq_stp(temp=temp,pres=pres,salt=salt,dliq=dliq,
        chkvals=chkvals,chktol=chktol,temp0=temp0,dliq0=dliq0,chkbnd=chkbnd,
        useext=useext,mathargs=mathargs)
    return temp

def dtfdp(salt,pres,temp=None,dliq=None,chkvals=False,chktol=_CHKTOL,
    temp0=None,dliq0=None,chkbnd=False,useext=False,mathargs=None):
    """Calculate sea-ice freezing point pressure lowering.
    
    Calculate the effect of pressure on lowering the freezing point of
    sea-ice.
    
    :arg float salt: Salinity in kg/kg.
    :arg float pres: Pressure in Pa.
    :arg temp: Temperature in K. If unknown, pass None (default) and it
        will be calculated.
    :type temp: float or None
    :arg dliq: Density of liquid water in seawater in kg/m3. If
        unknown, pass None (default) and it will be calculated.
    :type dliq: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg temp0: Initial guess for the temperature in K. If None
        (default) then `_approx_sp` is used.
    :type temp0: float or None
    :arg dliq0: Initial guess for the liquid water density in kg/m3. If
        None (default) then `flu3a._dliq_default` is used.
    :type dliq0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg bool useext: If False (default) then the salt contribution is
        calculated from _GSCOEFFS; if True, from _GSCOEFFS_EXT.
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Freezing point lowering in K/Pa.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> dtfdp(0.035,1e5)
    7.48210942879e-8
    """
    __, temp, __, dliq = eq_stp(temp=temp,pres=pres,salt=salt,dliq=dliq,
        chkvals=chkvals,chktol=chktol,temp0=temp0,dliq0=dliq0,chkbnd=chkbnd,
        useext=useext,mathargs=mathargs)
    fl_t = _flu_f(1,0,temp,dliq)
    gs_t = _sal_g(0,1,0,salt,temp,pres,useext=useext)
    gs_p = _sal_g(0,0,1,salt,temp,pres,useext=useext)
    gs_st = _sal_g(1,1,0,salt,temp,pres,useext=useext)
    gs_sp = _sal_g(1,0,1,salt,temp,pres,useext=useext)
    gi_t = _ice_g(1,0,temp,pres)
    gi_p = _ice_g(0,1,temp,pres)
    dent = fl_t + gs_t - salt*gs_st - gi_t
    dvol = dliq**(-1) + gs_p - salt*gs_sp - gi_p
    dtfdp = dvol/dent
    return dtfdp

def dtfds(salt,pres,temp=None,dliq=None,chkvals=False,chktol=_CHKTOL,
    temp0=None,dliq0=None,chkbnd=False,useext=False,mathargs=None):
    """Calculate sea-ice freezing point salt lowering.
    
    Calculate the effect of salinity on lowering the freezing point of
    sea-ice.
    
    :arg float salt: Salinity in kg/kg.
    :arg float pres: Pressure in Pa.
    :arg temp: Temperature in K. If unknown, pass None (default) and it
        will be calculated.
    :type temp: float or None
    :arg dliq: Density of liquid water in seawater in kg/m3. If
        unknown, pass None (default) and it will be calculated.
    :type dliq: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg temp0: Initial guess for the temperature in K. If None
        (default) then `_approx_sp` is used.
    :type temp0: float or None
    :arg dliq0: Initial guess for the liquid water density in kg/m3. If
        None (default) then `flu3a._dliq_default` is used.
    :type dliq0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg bool useext: If False (default) then the salt contribution is
        calculated from _GSCOEFFS; if True, from _GSCOEFFS_EXT.
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Freezing point lowering in K/(kg/kg).
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> dtfds(0.035,1e5)
    -56.8751336296
    """
    __, temp, __, dliq = eq_stp(temp=temp,pres=pres,salt=salt,dliq=dliq,
        chkvals=chkvals,chktol=chktol,temp0=temp0,dliq0=dliq0,chkbnd=chkbnd,
        useext=useext,mathargs=mathargs)
    fl_t = _flu_f(1,0,temp,dliq)
    gs_t = _sal_g(0,1,0,salt,temp,pres,useext=useext)
    gs_ss = _sal_g(2,0,0,salt,temp,pres,useext=useext)
    gs_st = _sal_g(1,1,0,salt,temp,pres,useext=useext)
    gi_t = _ice_g(1,0,temp,pres)
    dent = fl_t + gs_t - salt*gs_st - gi_t
    dtfds = salt*gs_ss / dent
    return dtfds


## Seawater-ice combined system
def eq_seaice(sisal,temp,pres,salt=None,dliq=None,chkvals=False,
    chktol=_CHKTOL,salt0=None,dliq0=None,chkbnd=False,useext=False,
    mathargs=None):
    """Get primary sea-ice variables at SsiTP.
    
    Get the values of all primary variables for a seawater-ice parcel at
    the given total salinity, temperature, and pressure. Total salinity
    here is the ratio of the mass of salt to the total parcel mass (salt
    + liquid water + ice).
    
    If the calculation has already been done, the results can be passed
    to avoid unnecessary repeat calculations. If enough values are
    passed, they will be checked for consistency if chkvals is True.
    
    :arg float sisal: Total sea-ice salinity in kg/kg.
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg salt: Seawater salinity in kg/kg. If unknown, pass None
        (default) and it will be calculated.
    :type salt: float or None
    :arg dliq: Seawater liquid water density in kg/m3. If unknown, pass
        None (default) and it will be calculated.
    :type dliq: float or None
    :arg dvap: Water vapour density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dvap: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg salt0: Initial guess for the seawater salinity in kg/kg. If
        None (default) then `_approx_tp` is used.
    :type salt0: float or None
    :arg dliq0: Initial guess for the seawater liquid water density in
        kg/m3. If None (default) then `flu3a._dliq_default` is used.
    :type dliq0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg bool useext: If False (default) then the salt contribution is
        calculated from _GSCOEFFS; if True, from _GSCOEFFS_EXT.
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Seawater salinity and liquid water density (both in SI
        units).
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    :raises RuntimeWarning: If the equilibrium seawater salinity is
        lower than the total parcel salinity.
    """
    if salt is None or dliq is None:
        salt, __, __, dliq = eq_stp(temp=temp,pres=pres,salt=salt,dliq=dliq,
            chkvals=chkvals,chktol=chktol,salt0=salt0,dliq0=dliq0,chkbnd=chkbnd,
            useext=useext,mathargs=mathargs)
    if salt < sisal:
        warnmsg = ('Equilibrium salinity {0} is lower than the total parcel '
            'salinity {1}').format(salt,sisal)
        warnings.warn(warnmsg,RuntimeWarning)
        salt = sisal
    return salt, dliq

def seaice_g(drvs,drvt,drvp,sisal,temp,pres,salt=None,dliq=None,
    chkvals=False,chktol=_CHKTOL,salt0=None,dliq0=None,chkbnd=False,
    useext=False,mathargs=None):
    """Calculate sea-ice Gibbs free energy with derivatives.
    
    Calculate the specific Gibbs free energy of a sea-ice parcel or its
    derivatives with respect to total salinity, temperature, and
    pressure.
    
    :arg int drvs: Number of total salinity derivatives.
    :arg int drvt: Number of temperature derivatives.
    :arg int drvp: Number of pressure derivatives.
    :arg float sisal: Total sea-ice salinity in kg/kg.
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg salt: Seawater salinity in kg/kg. If unknown, pass None
        (default) and it will be calculated.
    :type salt: float or None
    :arg dliq: Seawater liquid water density in kg/m3. If unknown, pass
        None (default) and it will be calculated.
    :type dliq: float or None
    :arg dvap: Water vapour density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dvap: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg salt0: Initial guess for the seawater salinity in kg/kg. If
        None (default) then `_approx_tp` is used.
    :type salt0: float or None
    :arg dliq0: Initial guess for the seawater liquid water density in
        kg/m3. If None (default) then `flu3a._dliq_default` is used.
    :type dliq0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg bool useext: If False (default) then the salt contribution is
        calculated from _GSCOEFFS; if True, from _GSCOEFFS_EXT.
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Gibbs free energy in units of
        (J/kg) / (kg/kg)^drvs / K^drvt / Pa^drvp.
    :raises ValueError: If any of (drvs,drvt,drvp) are negative or if
        (drvs+drvt+drvp) > 2.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    :raises RuntimeWarning: If the equilibrium seawater salinity is
        lower than the total parcel salinity.
    
    :Examples:
    
    >>> seaice_g(0,0,0,0.035,270.,1e5)
    -414.017574547
    >>> seaice_g(1,0,0,0.035,270.,1e5)
    96363.7730495
    >>> seaice_g(0,1,0,0.035,270.,1e5)
    500.445444181
    >>> seaice_g(0,0,1,0.035,270.,1e5)
    1.00689072300e-3
    >>> seaice_g(2,0,0,0.035,270.,1e5)
    0.
    >>> seaice_g(1,1,0,0.035,270.,1e5)
    -1144.02883419
    >>> seaice_g(1,0,1,0.035,270.,1e5)
    -8.62856321467e-4
    >>> seaice_g(0,2,0,0.035,270.,1e5)
    -232.847783380
    >>> seaice_g(0,1,1,0.035,270.,1e5)
    -1.65866446694e-5
    >>> seaice_g(0,0,2,0.035,270.,1e5)
    -1.57591932118e-12
    """
    drvtup = (drvs,drvt,drvp)
    if any(drv < 0 for drv in drvtup) or sum(drvtup) > 2:
        errmsg = 'Derivatives {0} not recognized'.format(drvtup)
        raise ValueError(errmsg)
    salt, dliq = eq_seaice(sisal,temp,pres,salt=salt,dliq=dliq,chkvals=chkvals,
        chktol=chktol,salt0=salt0,dliq0=dliq0,chkbnd=chkbnd,useext=useext,
        mathargs=mathargs)
    seaf = sisal/salt
    
    # Straightforward derivatives
    if (drvs,drvt,drvp) == (0,0,0):
        gl = _eq_chempot(0,0,temp,dliq)
        gs = _sal_g(0,0,0,salt,temp,pres,useext=useext)
        gi = _ice_g(0,0,temp,pres)
        g = seaf*(gl + gs) + (1-seaf)*gi
        return g
    elif (drvs,drvt,drvp) == (1,0,0):
        gs_s = _sal_g(1,0,0,salt,temp,pres,useext=useext)
        g_s = gs_s
        return g_s
    elif (drvs,drvt,drvp) == (0,1,0):
        fl_t = _flu_f(1,0,temp,dliq)
        gs_t = _sal_g(0,1,0,salt,temp,pres,useext=useext)
        gi_t = _ice_g(1,0,temp,pres)
        g_t = seaf*(fl_t + gs_t) + (1-seaf)*gi_t
        return g_t
    elif (drvs,drvt,drvp) == (0,0,1):
        gs_p = _sal_g(0,0,1,salt,temp,pres,useext=useext)
        gi_p = _ice_g(0,1,temp,pres)
        g_p = seaf*(dliq**(-1) + gs_p) + (1-seaf)*gi_p
        return g_p
    elif (drvs,drvt,drvp) == (2,0,0):
        g_ss = 0.0
        return g_ss
    elif (drvs,drvt,drvp) == (1,1,0):
        fl_t = _flu_f(1,0,temp,dliq)
        gs_t = _sal_g(0,1,0,salt,temp,pres,useext=useext)
        gi_t = _ice_g(1,0,temp,pres)
        g_st = (fl_t + gs_t - gi_t)/salt
        return g_st
    elif (drvs,drvt,drvp) == (1,0,1):
        gs_p = _sal_g(0,0,1,salt,temp,pres,useext=useext)
        gi_p = _ice_g(0,1,temp,pres)
        g_sp = (dliq**(-1) + gs_p - gi_p)/salt
        return g_sp
    
    # Other derivatives require inversion
    cl = _eq_pressure(0,1,temp,dliq)
    gs_ss = _sal_g(2,0,0,salt,temp,pres,useext=useext)
    if drvt > 0:
        fl_t = _flu_f(1,0,temp,dliq)
        gs_t = _sal_g(0,1,0,salt,temp,pres,useext=useext)
        gs_st = _sal_g(1,1,0,salt,temp,pres,useext=useext)
        gi_t = _ice_g(1,0,temp,pres)
        dentr = fl_t + gs_t - salt*gs_st - gi_t
    if drvp > 0:
        gs_p = _sal_g(0,0,1,salt,temp,pres,useext=useext)
        gs_sp = _sal_g(1,0,1,salt,temp,pres,useext=useext)
        gi_p = _ice_g(0,1,temp,pres)
        dvol = dliq**(-1) + gs_p - salt*gs_sp - gi_p
        s_p = dvol / (salt*gs_ss)
        dl_p = cl**(-1)
    
    if (drvs,drvt,drvp) == (0,2,0):
        fl_tt = _flu_f(2,0,temp,dliq)
        fl_td = _flu_f(1,1,temp,dliq)
        gs_tt = _sal_g(0,2,0,salt,temp,pres,useext=useext)
        gi_tt = _ice_g(2,0,temp,pres)
        
        s_t = dentr / (salt*gs_ss)
        dl_t = -dliq**2*fl_td/cl
        gb_tt = fl_tt + fl_td*dl_t + gs_tt
        g_tt = -seaf/salt*dentr*s_t + seaf*gb_tt + (1-seaf)*gi_tt
        return g_tt
    elif (drvs,drvt,drvp) == (0,1,1):
        fl_td = _flu_f(1,1,temp,dliq)
        gs_tp = _sal_g(0,1,1,salt,temp,pres,useext=useext)
        gi_tp = _ice_g(1,1,temp,pres)
        gb_tp = fl_td*dl_p + gs_tp
        g_tp = -seaf/salt*s_p + seaf*gb_tp + (1-seaf)*gi_tp
        return g_tp
    elif (drvs,drvt,drvp) == (0,0,2):
        gs_pp = _sal_g(0,0,2,salt,temp,pres,useext=useext)
        gi_pp = _ice_g(0,2,temp,pres)
        gb_pp = -dl_p/dliq**2 + gs_pp
        g_pp = -seaf/salt*dvol*s_p + seaf*gb_pp + (1-seaf)*gi_pp
        return g_pp
    
    # Should not have made it this far!
    errmsg = 'Derivatives {0} not recognized'.format((drvs,drvt,drvp))
    raise ValueError(errmsg)

def brinefraction(sisal,temp,pres,salt=None,dliq=None,chkvals=False,
    chktol=_CHKTOL,salt0=None,dliq0=None,chkbnd=False,useext=False,
    mathargs=None):
    """Calculate sea-ice brine fraction.
    
    Calculate the mass fraction of seawater (brine) in a sea-ice parcel,
    the ratio of the mass of seawater (salt + liquid water) to the total
    mass (salt + liquid water + ice).
    
    :arg float sisal: Total sea-ice salinity in kg/kg.
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg salt: Seawater salinity in kg/kg. If unknown, pass None
        (default) and it will be calculated.
    :type salt: float or None
    :arg dliq: Seawater liquid water density in kg/m3. If unknown, pass
        None (default) and it will be calculated.
    :type dliq: float or None
    :arg dvap: Water vapour density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dvap: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg salt0: Initial guess for the seawater salinity in kg/kg. If
        None (default) then `_approx_tp` is used.
    :type salt0: float or None
    :arg dliq0: Initial guess for the seawater liquid water density in
        kg/m3. If None (default) then `flu3a._dliq_default` is used.
    :type dliq0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg bool useext: If False (default) then the salt contribution is
        calculated from _GSCOEFFS; if True, from _GSCOEFFS_EXT.
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Brine fraction in kg/kg.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    :raises RuntimeWarning: If the equilibrium seawater salinity is
        lower than the total parcel salinity.
    
    :Examples:
    
    >>> brinefraction(0.035,270.,1e5)
    0.624705328368
    """
    salt, dliq = eq_seaice(sisal,temp,pres,salt=salt,dliq=dliq,chkvals=chkvals,
        chktol=chktol,salt0=salt0,dliq0=dliq0,chkbnd=chkbnd,useext=useext,
        mathargs=mathargs)
    seaf = sisal/salt
    return seaf

def cp(sisal,temp,pres,salt=None,dliq=None,chkvals=False,chktol=_CHKTOL,
    salt0=None,dliq0=None,chkbnd=False,useext=False,mathargs=None):
    """Calculate sea-ice isobaric heat capacity.
    
    Calculate the isobaric heat capacity of sea-ice.
    
    :arg float sisal: Total sea-ice salinity in kg/kg.
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg salt: Seawater salinity in kg/kg. If unknown, pass None
        (default) and it will be calculated.
    :type salt: float or None
    :arg dliq: Seawater liquid water density in kg/m3. If unknown, pass
        None (default) and it will be calculated.
    :type dliq: float or None
    :arg dvap: Water vapour density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dvap: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg salt0: Initial guess for the seawater salinity in kg/kg. If
        None (default) then `_approx_tp` is used.
    :type salt0: float or None
    :arg dliq0: Initial guess for the seawater liquid water density in
        kg/m3. If None (default) then `flu3a._dliq_default` is used.
    :type dliq0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg bool useext: If False (default) then the salt contribution is
        calculated from _GSCOEFFS; if True, from _GSCOEFFS_EXT.
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Heat capacity in J/kg/K.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    :raises RuntimeWarning: If the equilibrium seawater salinity is
        lower than the total parcel salinity.
    
    :Examples:
    
    >>> cp(0.035,270.,1e5)
    62868.9015126
    """
    g_tt = seaice_g(0,2,0,sisal,temp,pres,salt=salt,dliq=dliq,chkvals=chkvals,
        chktol=chktol,salt0=salt0,dliq0=dliq0,chkbnd=chkbnd,useext=useext,
        mathargs=mathargs)
    cp = -temp * g_tt
    return cp

def density(sisal,temp,pres,salt=None,dliq=None,chkvals=False,
    chktol=_CHKTOL,salt0=None,dliq0=None,chkbnd=False,useext=False,
    mathargs=None):
    """Calculate sea-ice total density.
    
    Calculate the total density of a sea-ice parcel.
    
    :arg float sisal: Total sea-ice salinity in kg/kg.
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg salt: Seawater salinity in kg/kg. If unknown, pass None
        (default) and it will be calculated.
    :type salt: float or None
    :arg dliq: Seawater liquid water density in kg/m3. If unknown, pass
        None (default) and it will be calculated.
    :type dliq: float or None
    :arg dvap: Water vapour density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dvap: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg salt0: Initial guess for the seawater salinity in kg/kg. If
        None (default) then `_approx_tp` is used.
    :type salt0: float or None
    :arg dliq0: Initial guess for the seawater liquid water density in
        kg/m3. If None (default) then `flu3a._dliq_default` is used.
    :type dliq0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg bool useext: If False (default) then the salt contribution is
        calculated from _GSCOEFFS; if True, from _GSCOEFFS_EXT.
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Density in kg/m3.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    :raises RuntimeWarning: If the equilibrium seawater salinity is
        lower than the total parcel salinity.
    
    :Examples:
    
    >>> density(0.035,270.,1e5)
    993.156434117
    """
    g_p = seaice_g(0,0,1,sisal,temp,pres,salt=salt,dliq=dliq,chkvals=chkvals,
        chktol=chktol,salt0=salt0,dliq0=dliq0,chkbnd=chkbnd,useext=useext,
        mathargs=mathargs)
    rho = g_p**(-1)
    return rho

def enthalpy(sisal,temp,pres,salt=None,dliq=None,chkvals=False,
    chktol=_CHKTOL,salt0=None,dliq0=None,chkbnd=False,useext=False,
    mathargs=None):
    """Calculate sea-ice enthalpy.
    
    Calculate the specific enthalpy of a sea-ice parcel.
    
    :arg float sisal: Total sea-ice salinity in kg/kg.
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg salt: Seawater salinity in kg/kg. If unknown, pass None
        (default) and it will be calculated.
    :type salt: float or None
    :arg dliq: Seawater liquid water density in kg/m3. If unknown, pass
        None (default) and it will be calculated.
    :type dliq: float or None
    :arg dvap: Water vapour density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dvap: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg salt0: Initial guess for the seawater salinity in kg/kg. If
        None (default) then `_approx_tp` is used.
    :type salt0: float or None
    :arg dliq0: Initial guess for the seawater liquid water density in
        kg/m3. If None (default) then `flu3a._dliq_default` is used.
    :type dliq0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg bool useext: If False (default) then the salt contribution is
        calculated from _GSCOEFFS; if True, from _GSCOEFFS_EXT.
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Enthalpy in J/kg.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    :raises RuntimeWarning: If the equilibrium seawater salinity is
        lower than the total parcel salinity.
    
    :Examples:
    
    >>> enthalpy(0.035,270.,1e5)
    -135534.287503
    """
    salt, dliq = eq_seaice(sisal,temp,pres,salt=salt,dliq=dliq,chkvals=chkvals,
        chktol=chktol,salt0=salt0,dliq0=dliq0,chkbnd=chkbnd,useext=useext,
        mathargs=mathargs)
    g = seaice_g(0,0,0,sisal,temp,pres,salt=salt,dliq=dliq,useext=useext)
    g_t = seaice_g(0,1,0,sisal,temp,pres,salt=salt,dliq=dliq,useext=useext)
    h = g - temp*g_t
    return h

def entropy(sisal,temp,pres,salt=None,dliq=None,chkvals=False,
    chktol=_CHKTOL,salt0=None,dliq0=None,chkbnd=False,useext=False,
    mathargs=None):
    """Calculate sea-ice entropy.
    
    Calculate the specific entropy of a sea-ice parcel.
    
    :arg float sisal: Total sea-ice salinity in kg/kg.
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg salt: Seawater salinity in kg/kg. If unknown, pass None
        (default) and it will be calculated.
    :type salt: float or None
    :arg dliq: Seawater liquid water density in kg/m3. If unknown, pass
        None (default) and it will be calculated.
    :type dliq: float or None
    :arg dvap: Water vapour density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dvap: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg salt0: Initial guess for the seawater salinity in kg/kg. If
        None (default) then `_approx_tp` is used.
    :type salt0: float or None
    :arg dliq0: Initial guess for the seawater liquid water density in
        kg/m3. If None (default) then `flu3a._dliq_default` is used.
    :type dliq0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg bool useext: If False (default) then the salt contribution is
        calculated from _GSCOEFFS; if True, from _GSCOEFFS_EXT.
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Entropy in J/kg/K.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    :raises RuntimeWarning: If the equilibrium seawater salinity is
        lower than the total parcel salinity.
    
    :Examples:
    
    >>> entropy(0.035,270.,1e5)
    -500.445444181
    """
    g_t = seaice_g(0,1,0,sisal,temp,pres,salt=salt,dliq=dliq,chkvals=chkvals,
        chktol=chktol,salt0=salt0,dliq0=dliq0,chkbnd=chkbnd,useext=useext,
        mathargs=mathargs)
    s = -g_t
    return s

def expansion(sisal,temp,pres,salt=None,dliq=None,chkvals=False,
    chktol=_CHKTOL,salt0=None,dliq0=None,chkbnd=False,useext=False,
    mathargs=None):
    """Calculate sea-ice thermal expansion coefficient.
    
    Calculate the thermal expansion coefficient of a sea-ice parcel.
    
    :arg float sisal: Total sea-ice salinity in kg/kg.
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg salt: Seawater salinity in kg/kg. If unknown, pass None
        (default) and it will be calculated.
    :type salt: float or None
    :arg dliq: Seawater liquid water density in kg/m3. If unknown, pass
        None (default) and it will be calculated.
    :type dliq: float or None
    :arg dvap: Water vapour density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dvap: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg salt0: Initial guess for the seawater salinity in kg/kg. If
        None (default) then `_approx_tp` is used.
    :type salt0: float or None
    :arg dliq0: Initial guess for the seawater liquid water density in
        kg/m3. If None (default) then `flu3a._dliq_default` is used.
    :type dliq0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg bool useext: If False (default) then the salt contribution is
        calculated from _GSCOEFFS; if True, from _GSCOEFFS_EXT.
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Expansion coefficient in 1/K.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    :raises RuntimeWarning: If the equilibrium seawater salinity is
        lower than the total parcel salinity.
    
    :Examples:
    
    >>> expansion(0.035,270.,1e5)
    -1.64731328738
    """
    salt, dliq = eq_seaice(sisal,temp,pres,salt=salt,dliq=dliq,chkvals=chkvals,
        chktol=chktol,salt0=salt0,dliq0=dliq0,chkbnd=chkbnd,useext=useext,
        mathargs=mathargs)
    g_p = seaice_g(0,0,1,sisal,temp,pres,salt=salt,dliq=dliq,useext=useext)
    g_tp = seaice_g(0,1,1,sisal,temp,pres,salt=salt,dliq=dliq,useext=useext)
    alpha = g_tp / g_p
    return alpha

def kappa_t(sisal,temp,pres,salt=None,dliq=None,chkvals=False,
    chktol=_CHKTOL,salt0=None,dliq0=None,chkbnd=False,useext=False,
    mathargs=None):
    """Calculate sea-ice isothermal compressibility.
    
    Calculate the isothermal compressibility of a sea-ice parcel.
    
    :arg float sisal: Total sea-ice salinity in kg/kg.
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg salt: Seawater salinity in kg/kg. If unknown, pass None
        (default) and it will be calculated.
    :type salt: float or None
    :arg dliq: Seawater liquid water density in kg/m3. If unknown, pass
        None (default) and it will be calculated.
    :type dliq: float or None
    :arg dvap: Water vapour density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dvap: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg salt0: Initial guess for the seawater salinity in kg/kg. If
        None (default) then `_approx_tp` is used.
    :type salt0: float or None
    :arg dliq0: Initial guess for the seawater liquid water density in
        kg/m3. If None (default) then `flu3a._dliq_default` is used.
    :type dliq0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg bool useext: If False (default) then the salt contribution is
        calculated from _GSCOEFFS; if True, from _GSCOEFFS_EXT.
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Compressibility in 1/Pa.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    :raises RuntimeWarning: If the equilibrium seawater salinity is
        lower than the total parcel salinity.
    
    :Examples:
    
    >>> kappa_t(0.035,270.,1e5)
    1.56513441348e-9
    """
    salt, dliq = eq_seaice(sisal,temp,pres,salt=salt,dliq=dliq,chkvals=chkvals,
        chktol=chktol,salt0=salt0,dliq0=dliq0,chkbnd=chkbnd,useext=useext,
        mathargs=mathargs)
    g_p = seaice_g(0,0,1,sisal,temp,pres,salt=salt,dliq=dliq,useext=useext)
    g_pp = seaice_g(0,0,2,sisal,temp,pres,salt=salt,dliq=dliq,useext=useext)
    kappa = -g_pp / g_p
    return kappa

