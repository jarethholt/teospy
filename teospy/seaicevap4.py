"""Seawater-ice-water vapour equilibrium functions.

This module provides properties of seawater in equilibrium with ice and
water vapour (sea-ice-vapour).

:Examples:

>>> densityvap(salt=0.035)
4.17156419318e-3
>>> pressure(salt=0.035)
521.950349225
>>> temperature(salt=0.035)
271.247815057
>>> densityvap(temp=270.)
3.77406140772e-3
>>> pressure(temp=270.)
470.059067981
>>> salinity(temp=270.)
0.05614892885
>>> densityvap(pres=500.)
4.00364833230e-3
>>> salinity(pres=500.)
0.0438955766482
>>> temperature(pres=500.)
270.734430917

:Functions:

* :func:`eq_stp`: Calculate primary variables for sea-ice-vapour at the
  seawater salinity, temperature, or pressure.
* :func:`densityvap`: Sea ice vapour water vapour density.
* :func:`pressure`: Sea ice vapour pressure.
* :func:`salinity`: Sea ice vapour salinity.
* :func:`temperature`: Sea ice vapour temperature1

"""

__all__ = ['eq_stp','densityvap','pressure','salinity','temperature']

import numpy
import warnings
import constants0
import flu1
import ice1
import flu2
import sal2
import maths3
import flu3a

_CHKTOL = constants0.CHKTOL
_RWAT = constants0.RWAT
_MSAL = constants0.MSAL
_RUNIV = constants0.RUNIV
_TTP = constants0.TTP
_PTPE = constants0.PTPE
_LLVTP = constants0.LLVTP
_LILTP = constants0.LILTP
_RSAL = _RUNIV / _MSAL
_EPSS = _RSAL / _RWAT
_AVL = _LLVTP/(_RWAT*_TTP)
_ALI = _LILTP/(_RWAT*_TTP)
_chkflubnds = constants0.chkflubnds
_chkicebnds = constants0.chkicebnds
_chksalbnds = constants0.chksalbnds
_flu_f = flu1.flu_f
_ice_g = ice1.ice_g
_eq_chempot = flu2.eq_chempot
_eq_pressure = flu2.eq_pressure
_sal_g = sal2.sal_g
_eq_liqpot = sal2.eq_liqpot
_newton = maths3.newton
_dliq_default = flu3a._dliq_default
_dvap_default = flu3a._dvap_default


## Equilibrium functions
def _approx_s(salt):
    """Approximate TPDlDv at S.
    
    Approximate the temperature, pressure, liquid water density, and
    water vapour density for sea-ice-vapour at the given salinity.
    
    :arg float salt: Salinity in kg/kg.
    :returns: Temperature, pressure, liquid water density, and water
        vapour density (all in SI units).
    """
    x = _EPSS/_ALI * salt
    y = (_AVL+_ALI) * x
    temp = _TTP*(1 - x)
    pres = _PTPE*(1 - y)
    dliq = _dliq_default(temp,pres)
    dvap = _dvap_default(temp,pres)
    return temp, pres, dliq, dvap

def _approx_t(temp):
    """Approximate SPDlDv at T.
    
    Approximate the salinity, pressure, liquid water density, and water
    vapour density for sea-ice-vapour at the given temperature.
    
    :arg float temp: Temperature in K.
    :returns: Salinity, pressure, liquid water density, and water vapour
        density (all in SI units).
    """
    x = 1 - temp/_TTP
    salt = _ALI/_EPSS*x
    y = (_AVL+_ALI) * x
    pres = _PTPE*(1 - y)
    dliq = _dliq_default(temp,pres)
    dvap = _dvap_default(temp,pres)
    return salt, pres, dliq, dvap

def _approx_p(pres):
    """Approximate STDlDv at P.
    
    Approximate the salinity, temperature, liquid water density, and
    water vapour density for sea-ice-vapour at the given pressure.
    
    :arg float salt: Pressure in Pa.
    :returns: Salinity, temperature, liquid water density, and water
        vapour density (all in SI units).
    """
    y = 1 - pres/_PTPE
    x = y / (_AVL+_ALI)
    temp = _TTP*(1 - x)
    salt = _ALI/_EPSS*x
    dliq = _dliq_default(temp,pres)
    dvap = _dvap_default(temp,pres)
    return salt, temp, dliq, dvap

def _diff_s(t,p,dl,dv,salt,useext=False):
    """Calculate sea-ice-vapour disequilibrium at S.
    
    Calculate both sides of the equations
    
        given pressure = pressure in liquid water
        given pressure = pressure in water vapour
        chemical potential of ice = potential of liquid in seawater
        chemical potential of ice = potential of water vapour
    
    and their Jacobians with respect to temperature, pressure, liquid
    water density, and water vapour density. Solving these equations
    gives equilibrium values at the given salinity.
    
    :arg float t: Temperature in K.
    :arg float p: Pressure in Pa.
    :arg float dl: Seawater liquid water density in kg/m3.
    :arg float dv: Water vapour density in kg/m3.
    :arg float salt: Salinity in kg/kg.
    :arg bool useext: If False (default) then the salt contribution is
        calculated from _GSCOEFFS; if True, from _GSCOEFFS_EXT.
    :returns: Left-hand side of the equation, right-hand side,
        Jacobian of LHS, and Jacobian of RHS.
    :rtype: tuple(array(float))
    """
    pl = _eq_pressure(0,0,t,dl)
    pv = _eq_pressure(0,0,t,dv)
    muv = _eq_chempot(0,0,t,dv)
    mui = _ice_g(0,0,t,p)
    mul = _eq_chempot(0,0,t,dl)
    mul += _eq_liqpot(0,0,0,salt,t,p,useext=useext)
    lhs = numpy.array([p, p, mui, mui])
    rhs = numpy.array([pl, pv, mul, muv])
    
    pl_t = _eq_pressure(1,0,t,dl)
    pl_d = _eq_pressure(0,1,t,dl)
    pv_t = _eq_pressure(1,0,t,dv)
    pv_d = _eq_pressure(0,1,t,dv)
    mui_t = _ice_g(1,0,t,p)
    mui_p = _ice_g(0,1,t,p)
    mul_t = _eq_chempot(1,0,t,dl)
    mul_t += _eq_liqpot(0,1,0,salt,t,p,useext=useext)
    mul_p = _eq_liqpot(0,0,1,salt,t,p,useext=useext)
    mul_d = _eq_chempot(0,1,t,dl)
    muv_t = _eq_chempot(1,0,t,dv)
    muv_d = _eq_chempot(0,1,t,dv)
    dlhs = numpy.array([[0.,1.,0.,0.], [0.,1.,0.,0.], [mui_t,mui_p,0.,0.],
        [mui_t,mui_p,0.,0.]])
    drhs = numpy.array([[pl_t,0.,pl_d,0.], [pv_t,0.,0.,pv_d],
        [mul_t,mul_p,mul_d,0.], [muv_t,0.,0.,muv_d]])
    return lhs, rhs, dlhs, drhs

def _diff_t(s,p,dl,dv,temp,useext=False):
    """Calculate sea-ice-vapour disequilibrium at T.
    
    Calculate both sides of the equations
    
        pressure in liquid water = given pressure
        pressure in water vapour = given pressure
        chemical potential of water vapour = potential of ice
        chemical potential of water vapour = potential of liquid water
    
    and their Jacobians with respect to salinity, pressure, liquid water
    density, and water vapour density. Solving these equations gives
    equilibrium values at the given temperature.

    :arg float s: Salinity in kg/kg.
    :arg float p: Pressure in Pa.
    :arg float dl: Seawater liquid water density in kg/m3.
    :arg float dv: Water vapour density in kg/m3.
    :arg float temp: Temperature in K.
    :arg bool useext: If False (default) then the salt contribution is
        calculated from _GSCOEFFS; if True, from _GSCOEFFS_EXT.
    :returns: Left-hand side of the equation, right-hand side,
        Jacobian of LHS, and Jacobian of RHS.
    :rtype: tuple(array(float))
    """
    pl = _eq_pressure(0,0,temp,dl)
    pv = _eq_pressure(0,0,temp,dv)
    muv = _eq_chempot(0,0,temp,dv)
    mui = _ice_g(0,0,temp,p)
    mul = _eq_chempot(0,0,temp,dl)
    mul += _eq_liqpot(0,0,0,s,temp,p,useext=useext)
    lhs = numpy.array([pl, pv, muv, muv])
    rhs = numpy.array([p, p, mui, mul])
    
    pl_d = _eq_pressure(0,1,temp,dl)
    pv_d = _eq_pressure(0,1,temp,dv)
    muv_d = _eq_chempot(0,1,temp,dv)
    mui_p = _ice_g(0,1,temp,p)
    mul_s = _eq_liqpot(1,0,0,s,temp,p,useext=useext)
    mul_p = _eq_liqpot(0,0,1,s,temp,p,useext=useext)
    mul_d = _eq_chempot(0,1,temp,dl)
    dlhs = numpy.array([[0.,0.,pl_d,0.], [0.,0.,0.,pv_d], [0.,0.,0.,muv_d],
        [0.,0.,0.,muv_d]])
    drhs = numpy.array([[1.,1.,0.,0.], [1.,1.,0.,0.], [0.,mui_p,0.,0.],
        [mul_s,mul_p,mul_d,0.]])
    return lhs, rhs, dlhs, drhs

def _diff_p(s,t,dl,dv,pres,useext=False):
    """Calculate sea-ice-vapour disequilibrium at P.
    
    Calculate both sides of the equations
    
        pressure in liquid water = given pressure
        pressure in water vapour = given pressure
        chemical potential of water vapour = potential of ice
        chemical potential of water vapour = potential of liquid water
    
    and their Jacobians with respect to salinity, temperature, liquid
    water density, and water vapour density. Solving these equations
    gives equilibrium values at the given pressure.
    
    :arg float s: Salinity in kg/kg.
    :arg float t: Temperature in K.
    :arg float dl: Seawater liquid water density in kg/m3.
    :arg float dv: Water vapour density in kg/m3.
    :arg float pres: Pressure in Pa.
    :arg bool useext: If False (default) then the salt contribution is
        calculated from _GSCOEFFS; if True, from _GSCOEFFS_EXT.
    :returns: Left-hand side of the equation, right-hand side,
        Jacobian of LHS, and Jacobian of RHS.
    :rtype: tuple(array(float))
    """
    pl = _eq_pressure(0,0,t,dl)
    pv = _eq_pressure(0,0,t,dv)
    muv = _eq_chempot(0,0,t,dv)
    mui = _ice_g(0,0,t,pres)
    mul = _eq_chempot(0,0,t,dl)
    mul += _eq_liqpot(0,0,0,s,t,pres,useext=useext)
    lhs = numpy.array([pl, pv, muv, muv])
    rhs = numpy.array([pres, pres, mui, mul])
    
    pl_t = _eq_pressure(1,0,t,dl)
    pl_d = _eq_pressure(0,1,t,dl)
    pv_t = _eq_pressure(1,0,t,dv)
    pv_d = _eq_pressure(0,1,t,dv)
    muv_t = _eq_chempot(1,0,t,dv)
    muv_d = _eq_chempot(0,1,t,dv)
    mui_t = _ice_g(1,0,t,pres)
    mul_s = _eq_liqpot(1,0,0,s,t,pres,useext=useext)
    mul_t = _eq_chempot(1,0,t,dl)
    mul_t += _eq_liqpot(0,1,0,s,t,pres,useext=useext)
    mul_d = _eq_chempot(0,1,t,dl)
    dlhs = numpy.array([[0.,pl_t,pl_d,0.], [0.,pv_t,0.,pv_d],
        [0.,muv_t,0.,muv_d], [0.,muv_t,0.,muv_d]])
    drhs = numpy.array([[0.,0.,0.,0.], [0.,0.,0.,0.], [0.,mui_t,0.,0.],
        [mul_s,mul_t,mul_d,0.]])
    return lhs, rhs, dlhs, drhs

def eq_stp(salt=None,temp=None,pres=None,dliq=None,dvap=None,
    chkvals=False,chktol=_CHKTOL,salt0=None,temp0=None,pres0=None,
    dliq0=None,dvap0=None,chkbnd=False,useext=False,mathargs=None):
    """Get primary sea-ice-vapour variables at STP.
    
    Get the values of all primary variables for seawater, ice, and pure
    water vapour in equilibrium at any of the seawater salinity,
    temperature, and pressure.
    
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
    :arg dliq: Seawater liquid water density in kg/m3. If unknown,
        pass None (default) and it will be calculated.
    :type dliq: float or None
    :arg dvap: Water vapour density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dvap: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg salt0: Initial guess for the salinity in kg/kg. If None
        (default) then the appropriate `_approx_*` is used.
    :type salt0: float or None
    :arg temp0: Initial guess for the temperature in K. If None
        (default) then the appropriate `_approx_*` is used.
    :type temp0: float or None
    :arg pres0: Initial guess for the pressure in Pa. If None (default)
        then the appropriate `_approx_*` is used.
    :type pres0: float or None
    :arg dliq0: Initial guess for the seawater liquid water density in
        kg/m3. If None (default) then `flu3a._dliq_default` is used.
    :type dliq0: float or None
    :arg dvap0: Initial guess for the salinity in kg/kg. If None
        (default) then `flu3a._dvap_default` is used.
    :type dvap0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg bool useext: If False (default) then the salt contribution is
        calculated from _GSCOEFFS; if True, from _GSCOEFFS_EXT.
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Salinity, temperature, pressure, seawater liquid density,
        and water vapour density (all in SI units).
    :raises ValueError: If no values are provided.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    """
    if all(val is None for val in (salt,temp,pres)):
        errmsg = 'Must provide at least one of (salt,temp,pres)'
        raise ValueError(errmsg)
    if mathargs is None:
        mathargs = dict()
    fkwargs = {'useext': useext}
    
    if salt is not None:
        if any(val is None for val in (temp,pres)):
            x0 = (temp0,pres0,dliq0,dvap0)
            fargs = (salt,)
            x1 = _newton(_diff_s,x0,_approx_s,fargs=fargs,fkwargs=fkwargs,
                **mathargs)
            temp, pres, dliq, dvap = x1
    elif temp is not None:
        if any(val is None for val in (salt,pres)):
            x0 = (salt0,pres0,dliq0,dvap0)
            fargs = (temp,)
            x1 = _newton(_diff_t,x0,_approx_t,fargs=fargs,fkwargs=fkwargs,
                **mathargs)
            salt, pres, dliq, dvap = x1
    elif pres is not None:
        x0 = (salt0,temp0,dliq0,dvap0)
        fargs = (pres,)
        x1 = _newton(_diff_p,x0,_approx_p,fargs=fargs,fkwargs=fkwargs,
            **mathargs)
        salt, temp, dliq, dvap = x1
    if dliq is None:
        dliq = flu3a.eq_tp_liq(temp,pres,dliq0=dliq0,mathargs=mathargs)
    if dvap is None:
        dvap = flu3a.eq_tp_vap(temp,pres,dvap0=dvap0,mathargs=mathargs)
    
    _chkflubnds(temp,dliq,chkbnd=chkbnd)
    _chkflubnds(temp,dvap,chkbnd=chkbnd)
    _chkicebnds(temp,pres,chkbnd=chkbnd)
    _chksalbnds(salt,temp,pres,chkbnd=chkbnd)
    if not chkvals:
        return salt, temp, pres, dliq, dvap
    
    lhs, rhs, __, __ = _diff_s(temp,pres,dliq,dvap,salt,useext=useext)
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
    return salt, temp, pres, dliq, dvap


## Thermodynamic properties
def densityvap(salt=None,temp=None,pres=None,dliq=None,dvap=None,
    chkvals=False,chktol=_CHKTOL,salt0=None,temp0=None,pres0=None,
    dliq0=None,dvap0=None,chkbnd=False,useext=False,mathargs=None):
    """Calculate sea-ice-vapour water vapour density.
    
    Calculate the density of water vapour in sea-ice-vapour.
    
    :arg salt: Salinity in kg/kg. If unknown, pass None (default) and it
        will be calculated.
    :type salt: float or None
    :arg temp: Temperature in K. If unknown, pass None (default) and it
        will be calculated.
    :type temp: float or None
    :arg pres: Pressure in Pa. If unknown, pass None (default) and it
        will be calculated.
    :type pres: float or None
    :arg dliq: Seawater liquid water density in kg/m3. If unknown,
        pass None (default) and it will be calculated.
    :type dliq: float or None
    :arg dvap: Water vapour density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dvap: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg salt0: Initial guess for the salinity in kg/kg. If None
        (default) then the appropriate `_approx_*` is used.
    :type salt0: float or None
    :arg temp0: Initial guess for the temperature in K. If None
        (default) then the appropriate `_approx_*` is used.
    :type temp0: float or None
    :arg pres0: Initial guess for the pressure in Pa. If None (default)
        then the appropriate `_approx_*` is used.
    :type pres0: float or None
    :arg dliq0: Initial guess for the seawater liquid water density in
        kg/m3. If None (default) then `flu3a._dliq_default` is used.
    :type dliq0: float or None
    :arg dvap0: Initial guess for the salinity in kg/kg. If None
        (default) then `flu3a._dvap_default` is used.
    :type dvap0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg bool useext: If False (default) then the salt contribution is
        calculated from _GSCOEFFS; if True, from _GSCOEFFS_EXT.
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Water vapour density in kg/m3.
    :raises ValueError: If no values are provided.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> densityvap(salt=0.035)
    4.17156419318e-3
    >>> densityvap(temp=270.)
    3.77406140772e-3
    >>> densityvap(pres=500.)
    4.00364833230e-3
    """
    salt, temp, pres, dliq, dvap = eq_stp(salt=salt,temp=temp,pres=pres,
        dliq=dliq,dvap=dvap,chkvals=chkvals,chktol=chktol,salt0=salt0,
        temp0=temp0,pres0=pres0,dliq0=dliq0,dvap0=dvap0,chkbnd=chkbnd,
        useext=useext,mathargs=mathargs)
    return dvap

def pressure(salt=None,temp=None,pres=None,dliq=None,dvap=None,
    chkvals=False,chktol=_CHKTOL,salt0=None,temp0=None,pres0=None,
    dliq0=None,dvap0=None,chkbnd=False,useext=False,mathargs=None):
    """Calculate sea-ice-vapour pressure.
    
    Calculate the pressure in sea-ice-vapour.
    
    :arg salt: Salinity in kg/kg. If unknown, pass None (default) and it
        will be calculated.
    :type salt: float or None
    :arg temp: Temperature in K. If unknown, pass None (default) and it
        will be calculated.
    :type temp: float or None
    :arg pres: Pressure in Pa. If unknown, pass None (default) and it
        will be calculated.
    :type pres: float or None
    :arg dliq: Seawater liquid water density in kg/m3. If unknown,
        pass None (default) and it will be calculated.
    :type dliq: float or None
    :arg dvap: Water vapour density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dvap: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg salt0: Initial guess for the salinity in kg/kg. If None
        (default) then the appropriate `_approx_*` is used.
    :type salt0: float or None
    :arg temp0: Initial guess for the temperature in K. If None
        (default) then the appropriate `_approx_*` is used.
    :type temp0: float or None
    :arg pres0: Initial guess for the pressure in Pa. If None (default)
        then the appropriate `_approx_*` is used.
    :type pres0: float or None
    :arg dliq0: Initial guess for the seawater liquid water density in
        kg/m3. If None (default) then `flu3a._dliq_default` is used.
    :type dliq0: float or None
    :arg dvap0: Initial guess for the salinity in kg/kg. If None
        (default) then `flu3a._dvap_default` is used.
    :type dvap0: float or None
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
    :raises ValueError: If no values are provided.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> pressure(salt=0.035)
    521.950349225
    >>> pressure(temp=270.)
    470.059067981
    """
    salt, temp, pres, dliq, dvap = eq_stp(salt=salt,temp=temp,pres=pres,
        dliq=dliq,dvap=dvap,chkvals=chkvals,chktol=chktol,salt0=salt0,
        temp0=temp0,pres0=pres0,dliq0=dliq0,dvap0=dvap0,chkbnd=chkbnd,
        useext=useext,mathargs=mathargs)
    return pres

def salinity(salt=None,temp=None,pres=None,dliq=None,dvap=None,
    chkvals=False,chktol=_CHKTOL,salt0=None,temp0=None,pres0=None,
    dliq0=None,dvap0=None,chkbnd=False,useext=False,mathargs=None):
    """Calculate sea-ice-vapour salinity.
    
    Calculate the salinity of sea-ice-vapour.
    
    :arg salt: Salinity in kg/kg. If unknown, pass None (default) and it
        will be calculated.
    :type salt: float or None
    :arg temp: Temperature in K. If unknown, pass None (default) and it
        will be calculated.
    :type temp: float or None
    :arg pres: Pressure in Pa. If unknown, pass None (default) and it
        will be calculated.
    :type pres: float or None
    :arg dliq: Seawater liquid water density in kg/m3. If unknown,
        pass None (default) and it will be calculated.
    :type dliq: float or None
    :arg dvap: Water vapour density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dvap: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg salt0: Initial guess for the salinity in kg/kg. If None
        (default) then the appropriate `_approx_*` is used.
    :type salt0: float or None
    :arg temp0: Initial guess for the temperature in K. If None
        (default) then the appropriate `_approx_*` is used.
    :type temp0: float or None
    :arg pres0: Initial guess for the pressure in Pa. If None (default)
        then the appropriate `_approx_*` is used.
    :type pres0: float or None
    :arg dliq0: Initial guess for the seawater liquid water density in
        kg/m3. If None (default) then `flu3a._dliq_default` is used.
    :type dliq0: float or None
    :arg dvap0: Initial guess for the salinity in kg/kg. If None
        (default) then `flu3a._dvap_default` is used.
    :type dvap0: float or None
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
    :raises ValueError: If no values are provided.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> salinity(temp=270.)
    0.05614892885
    >>> salinity(pres=500.)
    0.0438955766482
    """
    salt, temp, pres, dliq, dvap = eq_stp(salt=salt,temp=temp,pres=pres,
        dliq=dliq,dvap=dvap,chkvals=chkvals,chktol=chktol,salt0=salt0,
        temp0=temp0,pres0=pres0,dliq0=dliq0,dvap0=dvap0,chkbnd=chkbnd,
        useext=useext,mathargs=mathargs)
    return salt

def temperature(salt=None,temp=None,pres=None,dliq=None,dvap=None,
    chkvals=False,chktol=_CHKTOL,salt0=None,temp0=None,pres0=None,
    dliq0=None,dvap0=None,chkbnd=False,useext=False,mathargs=None):
    """Calculate sea-ice-vapour temperature.
    
    Calculate the temperature of sea-ice-vapour.
    
    :arg salt: Salinity in kg/kg. If unknown, pass None (default) and it
        will be calculated.
    :type salt: float or None
    :arg temp: Temperature in K. If unknown, pass None (default) and it
        will be calculated.
    :type temp: float or None
    :arg pres: Pressure in Pa. If unknown, pass None (default) and it
        will be calculated.
    :type pres: float or None
    :arg dliq: Seawater liquid water density in kg/m3. If unknown,
        pass None (default) and it will be calculated.
    :type dliq: float or None
    :arg dvap: Water vapour density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dvap: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg salt0: Initial guess for the salinity in kg/kg. If None
        (default) then the appropriate `_approx_*` is used.
    :type salt0: float or None
    :arg temp0: Initial guess for the temperature in K. If None
        (default) then the appropriate `_approx_*` is used.
    :type temp0: float or None
    :arg pres0: Initial guess for the pressure in Pa. If None (default)
        then the appropriate `_approx_*` is used.
    :type pres0: float or None
    :arg dliq0: Initial guess for the seawater liquid water density in
        kg/m3. If None (default) then `flu3a._dliq_default` is used.
    :type dliq0: float or None
    :arg dvap0: Initial guess for the salinity in kg/kg. If None
        (default) then `flu3a._dvap_default` is used.
    :type dvap0: float or None
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
    :raises ValueError: If no values are provided.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> temperature(salt=0.035)
    271.247815057
    >>> temperature(pres=500.)
    270.734430917
    """
    salt, temp, pres, dliq, dvap = eq_stp(salt=salt,temp=temp,pres=pres,
        dliq=dliq,dvap=dvap,chkvals=chkvals,chktol=chktol,salt0=salt0,
        temp0=temp0,pres0=pres0,dliq0=dliq0,dvap0=dvap0,chkbnd=chkbnd,
        useext=useext,mathargs=mathargs)
    return temp

