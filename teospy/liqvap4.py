"""Liquid water-water vapour equilibrium functions.

This module provides thermodynamic properties for pure liquid water and
water vapour in equilibrium, e.g. the saturation vapour pressure. This
module can also be called as a function,

    python liqvap4.py

which will compare results from this module against reference values in
IAPWS 1995, table 8.

:Examples:

>>> pressure(temp=300.)
3536.80675227
>>> densityvap(temp=300.)
2.55896736829e-2
>>> densityliq(temp=300.)
996.513027468
>>> enthalpyevap(temp=300.)
2437289.24124
>>> entropyevap(temp=300.)
8124.29747080
>>> volumeevap(temp=300.)
39.0772595686
>>> temperature(pres=1e4)
318.956328924
>>> densityvap(pres=1e4)
6.81657223094e-2
>>> densityliq(pres=1e4)
989.833275365
>>> enthalpyevap(pres=1e4)
2392052.72723
>>> entropyevap(pres=1e4)
7499.62458905
>>> volumeevap(pres=1e4)
14.6691196141

:Functions:

* :func:`eq_tp`: Calculate liquid water-water vapour equilibrium
  properties at either temperature or pressure.
* :func:`temperature`: Liquid-vapour temperature.
* :func:`pressure`: Liquid-vapour temperature.
* :func:`densityvap`: Water vapour density at liquid-vapour equilibrium.
* :func:`densityliq`: Liquid water density at liquid-vapour equilibrium.
* :func:`chempot`: Liquid-vapour chemical potential.
* :func:`enthalpyevap`: Enthalpy of evaporation of water.
* :func:`enthalpyliq`: Liquid water enthalpy at liquid-vapour
  equilibrium.
* :func:`enthalpyvap`: Water vapour enthalpy at liquid-vapour
  equilibrium.
* :func:`entropyevap`: Entropy of evaporation of water.
* :func:`entropyliq`: Liquid water entropy at liquid-vapour equilibrium.
* :func:`entropyvap`: Water vapour entropy at liquid-vapour equilibrium.
* :func:`volumeevap`: Specific volume of evaporation of water.
* :func:`chkiapws95table8`: Check module against IAPWS 1995, table 8.

"""

__all__ = ['eq_tp','temperature','pressure','densityvap','densityliq','chempot',
    'enthalpyevap','enthalpyliq','enthalpyvap','entropyevap','entropyliq',
    'entropyvap','volumeevap','chkiapws95table8']

import warnings
import numpy
import constants0
import flu2
import maths3

_CHKTOL = constants0.CHKTOL
_RWAT = constants0.RWAT
_DCP = constants0.DCP
_TCP = constants0.TCP
_TTP = constants0.TTP
_PTPI = constants0.PTPI
_DVTP = constants0.DVTP
_DLTP = constants0.DLTP
_chkflubnds = constants0.chkflubnds
_eq_chempot = flu2.eq_chempot
_eq_pressure = flu2.eq_pressure
_newton = maths3.newton
_C_APPS = ((350,550),
    (1.80066818428501e-2,-0.648994409718973,1.56594764908347,-3.18116999660964,
        2.98590977093295),
    (-7.34017329598858e-2,5.70516487711065e-3,-4.31313846955949e-4),
    (-19.2235086866063,-6.15770193302955,-4.965736126494),
    (-0.237216002118091,0.186593118426901,-0.258472040504799),
    (-19.8731005709116,-3.08975437352998))


## Equilibrium functions
def _approx_t(temp):
    """Approximate DvDl at T.
    
    Approximate the densities of both water vapour and liquid water from
    the given temperature. This approximation is based on separate
    empirical equations for high and low temperatures.
    
    :arg float temp: Temperature in K.
    :returns: Vapour and liquid densities in kg/m3.
    :raises ValueError: If the temperature is above the critical point.
    """
    if temp > _TCP:
        errmsg = ('No liquid-vapour equilibrium exists above the critical '
            'point temperature {0} K').format(_TCP)
        raise ValueError(errmsg)
    
    TLIML, TLIMV = _C_APPS[0]
    if temp < TLIML:
        tau = temp/_TTP - 1
        dta = 1.
        for (i,a) in enumerate(_C_APPS[1]):
            dta += a * tau**(i+1)
        dliq = _DLTP*dta
    else:
        a1, a2, a3 = _C_APPS[2]
        r = a2 / a3
        s = a1 / a3
        t = (1 - temp/_TCP) / a3
        p = s - r**2/3
        q = 2*r**3/27 - r*s/3 + t
        rtp3 = (p/3)**.5
        t0 = -2*rtp3 * numpy.sinh(numpy.arcsinh(3*q/(2*p)/rtp3)/3)
        x0 = t0 - r/3
        x1 = numpy.copysign(abs(x0)**(1./3.),x0)
        dliq = _DCP * (1 + x1)
    
    if temp < TLIMV:
        tau = _TTP/temp - 1
        earg = 0.
        for (i,a) in enumerate(_C_APPS[3]):
            earg += a * tau**(i+1)
        dvap = _DVTP * numpy.exp(earg)
    else:
        a1, a2, a3 = _C_APPS[4]
        r = a2/a3
        s = a1/a3
        t = (1 - temp/_TCP) / a3
        p = s - r**2/3
        q = 2*r**3/27 - r*s/3 + t
        rtp3 = (p/3)**.5
        t0 = -2*rtp3*numpy.sinh(1./3.*numpy.arcsinh(3*q/(2*p)/rtp3))
        x0 = t0 - r/3
        x1 = abs(x0)**.25
        dvap = _DCP * (1 - x1)
    return dvap, dliq

def _approx_p(pres):
    """Approximate TDvDl at P.
    
    Approximate the temperature, water vapour density, and liquid water
    density for liquid-vapour equilibrium at the given pressure. This
    approximation is based on an empirical equation.
    
    :arg float pres: Pressure in Pa.
    :returns: Temperature, water vapour density, and liquid water
        density (all in SI units).
    """
    a1, a2 = _C_APPS[5]
    p = a1 / a2
    q = -numpy.log(pres / _PTPI) / a2
    tau = -p/2 + (p**2/4 - q)**.5
    temp = _TTP / (1 + tau)
    dvap, dliq = _approx_t(temp)
    return temp, dvap, dliq

def _diff_t(dv,dl,temp):
    """Calculate fluid water disequilibrium at T.
    
    Calculate both sides of the equations
    
        pressure in water vapour = pressure in liquid water
        chemical potential of water vapour = potential of liquid water
    
    and their Jacobians with respect to water vapour and liquid water
    density. Solving these equations gives the vapour and liquid
    densities for the given temperature.
    
    :arg float dv: Water vapour density in kg/m3.
    :arg float dl: Liquid water density in kg/m3.
    :arg float temp: Temperature in K.
    :returns: Left-hand side of the equation, right-hand side,
        Jacobian of LHS, and Jacobian of RHS.
    :rtype: tuple(array(float))
    """
    pv = _eq_pressure(0,0,temp,dv)
    pl = _eq_pressure(0,0,temp,dl)
    gv = _eq_chempot(0,0,temp,dv)
    gl = _eq_chempot(0,0,temp,dl)
    lhs = numpy.array([pv,gv])
    rhs = numpy.array([pl,gl])
    
    pv_d = _eq_pressure(0,1,temp,dv)
    pl_d = _eq_pressure(0,1,temp,dl)
    gv_d = _eq_chempot(0,1,temp,dv)
    gl_d = _eq_chempot(0,1,temp,dl)
    dlhs = numpy.array([[pv_d,0.], [gv_d,0.]])
    drhs = numpy.array([[0.,pl_d], [0.,gl_d]])
    return lhs, rhs, dlhs, drhs

def _diff_p(t,dv,dl,pres):
    """Calculate fluid water disequilibrium at P.
    
    Calculate both sides of the equations
    
        given pressure = pressure in water vapour
        given pressure = pressure in liquid water
        chemical potential of water vapour = potential of liquid water
    
    and their Jacobians with respect to temperature, water vapour
    density, and liquid water density. Solving these equations gives the
    temperature, vapour density, and liquid density for the given
    pressure.
    
    :arg float t: Temperature in K.
    :arg float dv: Water vapour density in kg/m3.
    :arg float dl: Liquid water density in kg/m3.
    :arg float pres: Pressure in Pa.
    :returns: Left-hand side of the equation, right-hand side,
        Jacobian of LHS, and Jacobian of RHS.
    :rtype: tuple(array(float))
    """
    pv = _eq_pressure(0,0,t,dv)
    pl = _eq_pressure(0,0,t,dl)
    gv = _eq_chempot(0,0,t,dv)
    gl = _eq_chempot(0,0,t,dl)
    lhs = numpy.array([pres, pres, gv])
    rhs = numpy.array([pv, pl, gl])
    
    pv_t = _eq_pressure(1,0,t,dv)
    pv_d = _eq_pressure(0,1,t,dv)
    pl_t = _eq_pressure(1,0,t,dl)
    pl_d = _eq_pressure(0,1,t,dl)
    gv_t = _eq_chempot(1,0,t,dv)
    gv_d = _eq_chempot(0,1,t,dv)
    gl_t = _eq_chempot(1,0,t,dl)
    gl_d = _eq_chempot(0,1,t,dl)
    dlhs = numpy.array([[0.,0.,0.], [0.,0.,0.], [gv_t,gv_d,0.]])
    drhs = numpy.array([[pv_t,pv_d,0.], [pl_t,0.,pl_d], [gl_t,0.,gl_d]])
    return lhs, rhs, dlhs, drhs

def eq_tp(temp=None,pres=None,dvap=None,dliq=None,chkvals=False,
    chktol=_CHKTOL,temp0=None,dvap0=None,dliq0=None,chkbnd=False,
    mathargs=None):
    """Get primary liquid-vapour variables at T or P.
    
    Get the values of all primary variables for liquid water and water
    vapour in equilibrium at either of a given temperature or pressure.
    
    If the calculation has already been done, the results can be passed
    to avoid unnecessary repeat calculations. If enough values are
    passed, they will be checked for consistency if chkvals is True.
    
    :arg temp: Temperature in K.
    :type temp: float or None
    :arg pres: Pressure in Pa.
    :type pres: float or None
    :arg dvap: Water vapour density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dvap: float or None
    :arg dliq: Liquid water density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dliq: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg temp0: Initial guess for the temperature in K. If None (default)
        then `_approx_p` is used.
    :type temp0: float or None
    :arg dvap0: Initial guess for the water vapour density in kg/m3. If
        None (default) then `_approx_t` is used.
    :type dvap0: float or None
    :arg dliq0: Initial guess for the liquid water density in kg/m3. If
        None (default) then `_approx_t` is used.
    :type dliq0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Temperature, pressure, water vapour density, and liquid
        water density (all in SI units).
    :raises ValueError: If neither of temp or pres is provided.
    :raises ValueError: If the temperature is above the critical point.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    """
    if temp is None and pres is None:
        errmsg = 'One of temp or pres must be given'
        raise ValueError(errmsg)
    if temp is not None:
        if any(val is None for val in (dvap,dliq)):
            x0 = (dvap0,dliq0)
            fargs = (temp,)
            if mathargs is None:
                mathargs = dict()
            x1 = _newton(_diff_t,x0,_approx_t,fargs=fargs,**mathargs)
            dvap, dliq = x1
        if pres is None:
            pres = flu2.pressure(temp,dvap)
    else:
        x0 = (temp0,dvap0,dliq0)
        fargs = (pres,)
        if mathargs is None:
            mathargs = dict()
        x1 = _newton(_diff_p,x0,_approx_p,fargs=fargs,**mathargs)
        temp, dvap, dliq = x1
    
    _chkflubnds(temp,dvap,chkbnd=chkbnd)
    _chkflubnds(temp,dliq,chkbnd=chkbnd)
    if not chkvals:
        return temp, pres, dvap, dliq
    
    lhs, rhs, __, __ = _diff_p(temp,dvap,dliq,pres)
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
    return temp, pres, dvap, dliq


## Equilibrium properties
def temperature(temp=None,pres=None,dvap=None,dliq=None,chkvals=False,
    chktol=_CHKTOL,temp0=None,dvap0=None,dliq0=None,chkbnd=False,
    mathargs=None):
    """Calculate liquid-vapour temperature.
    
    Calculate the temperature of liquid water and water vapour in
    equilibrium.
    
    :arg temp: Temperature in K.
    :type temp: float or None
    :arg pres: Pressure in Pa.
    :type pres: float or None
    :arg dvap: Water vapour density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dvap: float or None
    :arg dliq: Liquid water density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dliq: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg temp0: Initial guess for the temperature in K. If None (default)
        then `_approx_p` is used.
    :type temp0: float or None
    :arg dvap0: Initial guess for the water vapour density in kg/m3. If
        None (default) then `_approx_t` is used.
    :type dvap0: float or None
    :arg dliq0: Initial guess for the liquid water density in kg/m3. If
        None (default) then `_approx_t` is used.
    :type dliq0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Temperature in K.
    :raises ValueError: If neither of temp or pres is provided.
    :raises ValueError: If the temperature is above the critical point.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> temperature(pres=1e4)
    318.956328924
    """
    temp, pres, dvap, dliq = eq_tp(temp=temp,pres=pres,dvap=dvap,dliq=dliq,
        chkvals=chkvals,chktol=chktol,temp0=temp0,dvap0=dvap0,dliq0=dliq0,
        chkbnd=chkbnd,mathargs=mathargs)
    return temp

def pressure(temp=None,pres=None,dvap=None,dliq=None,chkvals=False,
    chktol=_CHKTOL,temp0=None,dvap0=None,dliq0=None,chkbnd=False,
    mathargs=None):
    """Calculate liquid-vapour pressure.
    
    Calculate the pressure of liquid water and water vapour in
    equilibrium.
    
    :arg temp: Temperature in K.
    :type temp: float or None
    :arg pres: Pressure in Pa.
    :type pres: float or None
    :arg dvap: Water vapour density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dvap: float or None
    :arg dliq: Liquid water density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dliq: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg temp0: Initial guess for the temperature in K. If None (default)
        then `_approx_p` is used.
    :type temp0: float or None
    :arg dvap0: Initial guess for the water vapour density in kg/m3. If
        None (default) then `_approx_t` is used.
    :type dvap0: float or None
    :arg dliq0: Initial guess for the liquid water density in kg/m3. If
        None (default) then `_approx_t` is used.
    :type dliq0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Pressure in Pa.
    :raises ValueError: If neither of temp or pres is provided.
    :raises ValueError: If the temperature is above the critical point.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> pressure(temp=300.)
    3536.80675227
    """
    temp, pres, dvap, dliq = eq_tp(temp=temp,pres=pres,dvap=dvap,dliq=dliq,
        chkvals=chkvals,chktol=chktol,temp0=temp0,dvap0=dvap0,dliq0=dliq0,
        chkbnd=chkbnd,mathargs=mathargs)
    return pres

def densityvap(temp=None,pres=None,dvap=None,dliq=None,chkvals=False,
    chktol=_CHKTOL,temp0=None,dvap0=None,dliq0=None,chkbnd=False,
    mathargs=None):
    """Calculate water vapour density.
    
    Calculate the density of water vapour for liquid water and water
    vapour in equilibrium.
    
    :arg temp: Temperature in K.
    :type temp: float or None
    :arg pres: Pressure in Pa.
    :type pres: float or None
    :arg dvap: Water vapour density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dvap: float or None
    :arg dliq: Liquid water density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dliq: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg temp0: Initial guess for the temperature in K. If None (default)
        then `_approx_p` is used.
    :type temp0: float or None
    :arg dvap0: Initial guess for the water vapour density in kg/m3. If
        None (default) then `_approx_t` is used.
    :type dvap0: float or None
    :arg dliq0: Initial guess for the liquid water density in kg/m3. If
        None (default) then `_approx_t` is used.
    :type dliq0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Water vapour density in kg/m3.
    :raises ValueError: If neither of temp or pres is provided.
    :raises ValueError: If the temperature is above the critical point.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> densityvap(temp=300.)
    2.55896736829e-2
    >>> densityvap(pres=1e4)
    6.81657223094e-2
    """
    temp, pres, dvap, dliq = eq_tp(temp=temp,pres=pres,dvap=dvap,dliq=dliq,
        chkvals=chkvals,chktol=chktol,temp0=temp0,dvap0=dvap0,dliq0=dliq0,
        chkbnd=chkbnd,mathargs=mathargs)
    return dvap

def densityliq(temp=None,pres=None,dvap=None,dliq=None,chkvals=False,
    chktol=_CHKTOL,temp0=None,dvap0=None,dliq0=None,chkbnd=False,
    mathargs=None):
    """Calculate liquid water density.
    
    Calculate the density of liquid water for liquid water and water
    vapour in equilibrium.
    
    :arg temp: Temperature in K.
    :type temp: float or None
    :arg pres: Pressure in Pa.
    :type pres: float or None
    :arg dvap: Water vapour density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dvap: float or None
    :arg dliq: Liquid water density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dliq: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg temp0: Initial guess for the temperature in K. If None (default)
        then `_approx_p` is used.
    :type temp0: float or None
    :arg dvap0: Initial guess for the water vapour density in kg/m3. If
        None (default) then `_approx_t` is used.
    :type dvap0: float or None
    :arg dliq0: Initial guess for the liquid water density in kg/m3. If
        None (default) then `_approx_t` is used.
    :type dliq0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Liquid water density in kg/m3.
    :raises ValueError: If neither of temp or pres is provided.
    :raises ValueError: If the temperature is above the critical point.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> densityliq(temp=300.)
    996.513027468
    >>> densityliq(pres=1e4)
    989.833275365
    """
    temp, pres, dvap, dliq = eq_tp(temp=temp,pres=pres,dvap=dvap,dliq=dliq,
        chkvals=chkvals,chktol=chktol,temp0=temp0,dvap0=dvap0,dliq0=dliq0,
        chkbnd=chkbnd,mathargs=mathargs)
    return dliq

def chempot(temp=None,pres=None,dvap=None,dliq=None,chkvals=False,
    chktol=_CHKTOL,temp0=None,dvap0=None,dliq0=None,chkbnd=False,
    mathargs=None):
    """Calculate liquid-vapour chemical potential.
    
    Calculate the chemical potential of liquid water and water vapour in
    equilibrium.
    
    :arg temp: Temperature in K.
    :type temp: float or None
    :arg pres: Pressure in Pa.
    :type pres: float or None
    :arg dvap: Water vapour density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dvap: float or None
    :arg dliq: Liquid water density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dliq: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg temp0: Initial guess for the temperature in K. If None (default)
        then `_approx_p` is used.
    :type temp0: float or None
    :arg dvap0: Initial guess for the water vapour density in kg/m3. If
        None (default) then `_approx_t` is used.
    :type dvap0: float or None
    :arg dliq0: Initial guess for the liquid water density in kg/m3. If
        None (default) then `_approx_t` is used.
    :type dliq0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Chemical potential in J/kg.
    :raises ValueError: If neither of temp or pres is provided.
    :raises ValueError: If the temperature is above the critical point.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> chempot(temp=300.)
    -5361.84908682
    >>> chempot(pres=1e4)
    -15259.1024273
    """
    temp, pres, dvap, dliq = eq_tp(temp=temp,pres=pres,dvap=dvap,dliq=dliq,
        chkvals=chkvals,chktol=chktol,temp0=temp0,dvap0=dvap0,dliq0=dliq0,
        chkbnd=chkbnd,mathargs=mathargs)
    g = _eq_chempot(0,0,temp,dvap)
    return g

def enthalpyevap(temp=None,pres=None,dvap=None,dliq=None,chkvals=False,
    chktol=_CHKTOL,temp0=None,dvap0=None,dliq0=None,chkbnd=False,
    mathargs=None):
    """Calculate enthalpy of evaporation.
    
    Calculate the enthalpy of evaporation for liquid water and water
    vapour in equilibrium.
    
    :arg temp: Temperature in K.
    :type temp: float or None
    :arg pres: Pressure in Pa.
    :type pres: float or None
    :arg dvap: Water vapour density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dvap: float or None
    :arg dliq: Liquid water density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dliq: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg temp0: Initial guess for the temperature in K. If None (default)
        then `_approx_p` is used.
    :type temp0: float or None
    :arg dvap0: Initial guess for the water vapour density in kg/m3. If
        None (default) then `_approx_t` is used.
    :type dvap0: float or None
    :arg dliq0: Initial guess for the liquid water density in kg/m3. If
        None (default) then `_approx_t` is used.
    :type dliq0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Enthalpy in J/kg.
    :raises ValueError: If neither of temp or pres is provided.
    :raises ValueError: If the temperature is above the critical point.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> enthalpyevap(temp=300.)
    2437289.24124
    >>> enthalpyevap(pres=1e4)
    2392052.72723
    """
    temp, pres, dvap, dliq = eq_tp(temp=temp,pres=pres,dvap=dvap,dliq=dliq,
        chkvals=chkvals,chktol=chktol,temp0=temp0,dvap0=dvap0,dliq0=dliq0,
        chkbnd=chkbnd,mathargs=mathargs)
    hv = flu2.enthalpy(temp,dvap)
    hl = flu2.enthalpy(temp,dliq)
    hevap = hv - hl
    return hevap

def enthalpyliq(temp=None,pres=None,dvap=None,dliq=None,chkvals=False,
    chktol=_CHKTOL,temp0=None,dvap0=None,dliq0=None,chkbnd=False,
    mathargs=None):
    """Calculate liquid water enthalpy.
    
    Calculate the enthalpy of liquid water for liquid water and water
    vapour in equilibrium.
    
    :arg temp: Temperature in K.
    :type temp: float or None
    :arg pres: Pressure in Pa.
    :type pres: float or None
    :arg dvap: Water vapour density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dvap: float or None
    :arg dliq: Liquid water density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dliq: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg temp0: Initial guess for the temperature in K. If None (default)
        then `_approx_p` is used.
    :type temp0: float or None
    :arg dvap0: Initial guess for the water vapour density in kg/m3. If
        None (default) then `_approx_t` is used.
    :type dvap0: float or None
    :arg dliq0: Initial guess for the liquid water density in kg/m3. If
        None (default) then `_approx_t` is used.
    :type dliq0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Enthalpy in J/kg.
    :raises ValueError: If neither of temp or pres is provided.
    :raises ValueError: If the temperature is above the critical point.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> enthalpyliq(temp=300.)
    112564.859854
    >>> enthalpyliq(pres=1e4)
    191805.944559
    """
    temp, pres, dvap, dliq = eq_tp(temp=temp,pres=pres,dvap=dvap,dliq=dliq,
        chkvals=chkvals,chktol=chktol,temp0=temp0,dvap0=dvap0,dliq0=dliq0,
        chkbnd=chkbnd,mathargs=mathargs)
    hl = flu2.enthalpy(temp,dliq)
    return hl

def enthalpyvap(temp=None,pres=None,dvap=None,dliq=None,chkvals=False,
    chktol=_CHKTOL,temp0=None,dvap0=None,dliq0=None,chkbnd=False,
    mathargs=None):
    """Calculate water vapour enthalpy.
    
    Calculate the enthalpy of water vapour for liquid water and water
    vapour in equilibrium.
    
    :arg temp: Temperature in K.
    :type temp: float or None
    :arg pres: Pressure in Pa.
    :type pres: float or None
    :arg dvap: Water vapour density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dvap: float or None
    :arg dliq: Liquid water density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dliq: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg temp0: Initial guess for the temperature in K. If None (default)
        then `_approx_p` is used.
    :type temp0: float or None
    :arg dvap0: Initial guess for the water vapour density in kg/m3. If
        None (default) then `_approx_t` is used.
    :type dvap0: float or None
    :arg dliq0: Initial guess for the liquid water density in kg/m3. If
        None (default) then `_approx_t` is used.
    :type dliq0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Enthalpy in J/kg.
    :raises ValueError: If neither of temp or pres is provided.
    :raises ValueError: If the temperature is above the critical point.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> enthalpyvap(temp=300.)
    2549854.10109
    >>> enthalpyvap(pres=1e4)
    2583858.67179
    """
    temp, pres, dvap, dliq = eq_tp(temp=temp,pres=pres,dvap=dvap,dliq=dliq,
        chkvals=chkvals,chktol=chktol,temp0=temp0,dvap0=dvap0,dliq0=dliq0,
        chkbnd=chkbnd,mathargs=mathargs)
    hv = flu2.enthalpy(temp,dvap)
    return hv

def entropyevap(temp=None,pres=None,dvap=None,dliq=None,chkvals=False,
    chktol=_CHKTOL,temp0=None,dvap0=None,dliq0=None,chkbnd=False,
    mathargs=None):
    """Calculate entropy of evaporation.
    
    Calculate the entropy of evaporation for liquid water and water
    vapour in equilibrium.
    
    :arg temp: Temperature in K.
    :type temp: float or None
    :arg pres: Pressure in Pa.
    :type pres: float or None
    :arg dvap: Water vapour density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dvap: float or None
    :arg dliq: Liquid water density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dliq: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg temp0: Initial guess for the temperature in K. If None (default)
        then `_approx_p` is used.
    :type temp0: float or None
    :arg dvap0: Initial guess for the water vapour density in kg/m3. If
        None (default) then `_approx_t` is used.
    :type dvap0: float or None
    :arg dliq0: Initial guess for the liquid water density in kg/m3. If
        None (default) then `_approx_t` is used.
    :type dliq0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Entropy in J/kg/K.
    :raises ValueError: If neither of temp or pres is provided.
    :raises ValueError: If the temperature is above the critical point.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> entropyevap(temp=300.)
    8124.29747080
    >>> entropyevap(pres=1e4)
    7499.62458905
    """
    temp, pres, dvap, dliq = eq_tp(temp=temp,pres=pres,dvap=dvap,dliq=dliq,
        chkvals=chkvals,chktol=chktol,temp0=temp0,dvap0=dvap0,dliq0=dliq0,
        chkbnd=chkbnd,mathargs=mathargs)
    sv = flu2.entropy(temp,dvap)
    sl = flu2.entropy(temp,dliq)
    sevap = sv - sl
    return sevap

def entropyliq(temp=None,pres=None,dvap=None,dliq=None,chkvals=False,
    chktol=_CHKTOL,temp0=None,dvap0=None,dliq0=None,chkbnd=False,
    mathargs=None):
    """Calculate liquid water entropy.
    
    Calculate the entropy of liquid water for liquid water and water
    vapour in equilibrium.
    
    :arg temp: Temperature in K.
    :type temp: float or None
    :arg pres: Pressure in Pa.
    :type pres: float or None
    :arg dvap: Water vapour density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dvap: float or None
    :arg dliq: Liquid water density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dliq: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg temp0: Initial guess for the temperature in K. If None (default)
        then `_approx_p` is used.
    :type temp0: float or None
    :arg dvap0: Initial guess for the water vapour density in kg/m3. If
        None (default) then `_approx_t` is used.
    :type dvap0: float or None
    :arg dliq0: Initial guess for the liquid water density in kg/m3. If
        None (default) then `_approx_t` is used.
    :type dliq0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Entropy in J/kg/K.
    :raises ValueError: If neither of temp or pres is provided.
    :raises ValueError: If the temperature is above the critical point.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> entropyliq(temp=300.)
    393.089029801
    >>> entropyliq(pres=1e4)
    649.195605196
    """
    temp, pres, dvap, dliq = eq_tp(temp=temp,pres=pres,dvap=dvap,dliq=dliq,
        chkvals=chkvals,chktol=chktol,temp0=temp0,dvap0=dvap0,dliq0=dliq0,
        chkbnd=chkbnd,mathargs=mathargs)
    sl = flu2.entropy(temp,dliq)
    return sl

def entropyvap(temp=None,pres=None,dvap=None,dliq=None,chkvals=False,
    chktol=_CHKTOL,temp0=None,dvap0=None,dliq0=None,chkbnd=False,
    mathargs=None):
    """Calculate water vapour entropy.
    
    Calculate the entropy of water vapour for liquid water and water
    vapour in equilibrium.
    
    :arg temp: Temperature in K.
    :type temp: float or None
    :arg pres: Pressure in Pa.
    :type pres: float or None
    :arg dvap: Water vapour density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dvap: float or None
    :arg dliq: Liquid water density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dliq: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg temp0: Initial guess for the temperature in K. If None (default)
        then `_approx_p` is used.
    :type temp0: float or None
    :arg dvap0: Initial guess for the water vapour density in kg/m3. If
        None (default) then `_approx_t` is used.
    :type dvap0: float or None
    :arg dliq0: Initial guess for the liquid water density in kg/m3. If
        None (default) then `_approx_t` is used.
    :type dliq0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Entropy in J/kg/K.
    :raises ValueError: If neither of temp or pres is provided.
    :raises ValueError: If the temperature is above the critical point.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> entropyvap(temp=300.)
    8517.38650061
    >>> entropyvap(pres=1e4)
    8148.82019425
    """
    temp, pres, dvap, dliq = eq_tp(temp=temp,pres=pres,dvap=dvap,dliq=dliq,
        chkvals=chkvals,chktol=chktol,temp0=temp0,dvap0=dvap0,dliq0=dliq0,
        chkbnd=chkbnd,mathargs=mathargs)
    sv = flu2.entropy(temp,dvap)
    return sv

def volumeevap(temp=None,pres=None,dvap=None,dliq=None,chkvals=False,
    chktol=_CHKTOL,temp0=None,dvap0=None,dliq0=None,chkbnd=False,
    mathargs=None):
    """Calculate specific volume of evaporation.
    
    Calculate the specific volume of evaporation for liquid water and
    water vapour in equilibrium.
    
    :arg temp: Temperature in K.
    :type temp: float or None
    :arg pres: Pressure in Pa.
    :type pres: float or None
    :arg dvap: Water vapour density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dvap: float or None
    :arg dliq: Liquid water density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dliq: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg temp0: Initial guess for the temperature in K. If None (default)
        then `_approx_p` is used.
    :type temp0: float or None
    :arg dvap0: Initial guess for the water vapour density in kg/m3. If
        None (default) then `_approx_t` is used.
    :type dvap0: float or None
    :arg dliq0: Initial guess for the liquid water density in kg/m3. If
        None (default) then `_approx_t` is used.
    :type dliq0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Specific volume in m3/kg.
    :raises ValueError: If neither of temp or pres is provided.
    :raises ValueError: If the temperature is above the critical point.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> volumeevap(temp=300.)
    39.0772595686
    >>> volumeevap(pres=1e4)
    14.6691196141
    """
    temp, pres, dvap, dliq = eq_tp(temp=temp,pres=pres,dvap=dvap,dliq=dliq,
        chkvals=chkvals,chktol=chktol,temp0=temp0,dvap0=dvap0,dliq0=dliq0,
        chkbnd=chkbnd,mathargs=mathargs)
    vv = dvap**(-1)
    vl = dliq**(-1)
    vevap = vv - vl
    return vevap


## Functions to check results
def chkiapws95table8(printresult=True,chktol=_CHKTOL):
    """Check accuracy against IAPWS 1995 table 8.
    
    Evaluate the functions in this module and compare to reference
    values from IAPWS 1995, table 8 for thermodynamic properties of
    liquid-vapour equilibrium.
    
    :arg bool printresult: If True (default) and any results are outside
        of the given tolerance, then the function name, reference value,
        result value, and relative error are printed.
    :arg float chktol: Tolerance to use when choosing to print results
        (default _LEMMONTOL). The default tolerance is lowered due to
        the low number of significant figures for the reference values.
    :returns: :class:`~tester.Tester` instance containing the
        functions, arguments, reference values, results, and relative
        errors from the tests.
    """
    from tester import Tester
    funs = [pressure,densityvap,densityliq,enthalpyvap,enthalpyliq,entropyvap,
        entropyliq]
    fargs = [(275.,), (450.,), (625.,)]
    refs = [
        [6.9845117e2,9.32203563628e5,1.6908269318578e7],
        [5.506649185041e-3,4.812003601257,1.182902804512e2],
        [999.887406120,890.34124976167,5.670903851464e2],
        [2.5042899500405e6,2.77441077988962e6,2.5507162456235e6],
        [7.759722016e3,7.4916158501217e5,1.6862697594697e6],
        [9.1066012052322e3,6.6092122132788e3,5.1850612079574e3],
        [2.83094669595e1,2.10865844688447e3,3.8019468301114e3]
    ]
    fnames = ['pressure','densityvap','densityliq','enthalpyvap','enthalpyliq',
        'entropyvap','entropyliq']
    argfmt = '(temp={0:3g})'
    header = 'Liq-vap equilibrium'
    test = Tester(funs,fargs,refs,fnames,argfmt,header=header)
    test.run()
    if printresult:
        test.printresults(chktol=chktol)
    return test


## Main function: Check tables
if __name__ == '__main__':
    test = chkiapws95table8()

