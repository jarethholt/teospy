"""Ice-liquid water equilibrium functions.

This module provides thermodynamic properties of ice and liquid water in
equilibrium, e.g. the enthalpy of melting.

:Examples:

>>> pressure(temp=270.)
39313338.8825
>>> densityliq(temp=270.)
1019.05568894
>>> enthalpymelt(temp=270.)
325166.686739
>>> entropymelt(temp=270.)
1204.32106199
>>> volumemelt(temp=270.)
-1.04052121182e-4
>>> temperature(pres=1e7)
272.401648869
>>> densityliq(pres=1e7)
1004.79353660
>>> enthalpymelt(pres=1e7)
331548.910817
>>> entropymelt(pres=1e7)
1217.13254011
>>> volumemelt(pres=1e7)
-9.42178903288e-5

:Functions:

* :func:`eq_tp`: Calculate ice-liquid water equilibrium properties at
  either temperature or pressure.
* :func:`temperature`: Temperature at ice-liquid water equilibrium.
* :func:`pressure`: Pressure at ice-liquid water equilibrium.
* :func:`densityliq`: Liquid water density at ice-liquid water
  equilibrium.
* :func:`chempot`: Chemical potential at ice-liquid water equilibrium.
* :func:`densityice`: Ice density at ice-liquid water equilibrium.
* :func:`enthalpyice`: Ice enthalpy at ice-liquid water equilibrium.
* :func:`enthalpyliq`: Liquid water enthalpy at ice-liquid water
  equilibrium.
* :func:`enthalpymelt`: Enthalpy of melting.
* :func:`entropyice`: Ice entropy at ice-liquid water equilibrium.
* :func:`entropyliq`: Liquid water entropy at ice-liquid water
  equilibrium.
* :func:`entropymelt`: Entropy of melting.
* :func:`volumemelt`: Specific volume of melting.

"""

__all__ = ['eq_tp','temperature','pressure','densityliq','chempot','densityice',
    'enthalpyice','enthalpyliq','enthalpymelt','entropyice','entropyliq',
    'entropymelt','volumemelt']

import warnings
import numpy
import constants0
import ice1
import flu2
import ice2
import maths3

_CHKTOL = constants0.CHKTOL
_TTP = constants0.TTP
_PTPI = constants0.PTPI
_DLTP = constants0.DLTP
_LILTP = constants0.LILTP
_chkflubnds = constants0.chkflubnds
_chkicebnds = constants0.chkicebnds
_ice_g = ice1.ice_g
_eq_chempot = flu2.eq_chempot
_eq_pressure = flu2.eq_pressure
_newton = maths3.newton
_C_APPS = ((-1.78582981492113,-12.2325084306734,-52.8236936433529),
    (-1.67329759176351e-7,-2.02262929999658e-13))


## Equilibrium functions
def _approx_t(temp):
    """Approximate PDl at T.
    
    Approximate the pressure and liquid water density for ice and liquid
    water in equilibrium at the given temperature. This approximation is
    based on an empirical polynomial for density.
    
    :arg float temp: Temperature in K.
    :returns: Pressure in Pa and liquid water density in kg/m3.
    """
    tau = temp/_TTP - 1
    dta = 0.
    for (i,a) in enumerate(_C_APPS[0]):
        dta += a * tau**(i+1)
    dliq = _DLTP * (1 + dta)
    pres = flu2.pressure(temp,dliq)
    return pres, dliq

def _approx_p(pres):
    """Approximate TDl at P.
    
    Approximate the temperature and liquid water density for ice and
    liquid water in equilibrium at the given pressure. This
    approximation is based on empirical polynomials for temperature and
    density.
    
    :arg float pres: Pressure in Pa.
    :returns: Temperature in K and liquid water density in kg/m3.
    """
    a1, a2 = _C_APPS[1]
    psi = pres/_PTPI - 1
    tau = a1*psi + a2*psi**2
    temp = _TTP * (1 + tau)
    
    dta = 0.
    for (i,a) in enumerate(_C_APPS[0]):
        dta += a * tau**(i+1)
    dliq = _DLTP * (1 + dta)
    return temp, dliq

def _diff_t(p,dl,temp):
    """Calculate ice-liquid disequilibrium at T.
    
    Calculate both sides of the equations
    
        given pressure = pressure of liquid water
        chemical potential of ice = potential of liquid water
    
    and their Jacobians with respect to pressure and liquid water
    density. Solving these equations gives the pressure and liquid water
    density at the given temperature.
    
    :arg float p: Pressure in Pa.
    :arg float dl: Liquid water density in kg/m3.
    :arg float temp: Temperature in K.
    :returns: Left-hand side of the equation, right-hand side,
        Jacobian of LHS, and Jacobian of RHS.
    :rtype: tuple(array(float))
    """
    pl = _eq_pressure(0,0,temp,dl)
    gi = _ice_g(0,0,temp,p)
    gl = _eq_chempot(0,0,temp,dl)
    lhs = numpy.array([p, gi])
    rhs = numpy.array([pl, gl])
    
    pl_d = _eq_pressure(0,1,temp,dl)
    gi_p = _ice_g(0,1,temp,p)
    gl_d = _eq_chempot(0,1,temp,dl)
    dlhs = numpy.array([[1.,0.], [gi_p,0.]])
    drhs = numpy.array([[0.,pl_d], [0.,gl_d]])
    return lhs, rhs, dlhs, drhs

def _diff_p(t,dl,pres):
    """Calculate ice-liquid disequilibrium at P.
    
    Calculate both sides of the equations
    
        given pressure = pressure of liquid water
        chemical potential of ice = potential of liquid water
    
    and their Jacobians with respect to temperature and liquid water
    density. Solving these equations gives the temperature and liquid
    water density at the given temperature.
    
    :arg float t: Temperature in K.
    :arg float dl: Liquid water density in kg/m3.
    :arg float pres: Pressure in Pa.
    :returns: Left-hand side of the equation, right-hand side,
        Jacobian of LHS, and Jacobian of RHS.
    :rtype: tuple(array(float))
    """
    pl = _eq_pressure(0,0,t,dl)
    gi = _ice_g(0,0,t,pres)
    gl = _eq_chempot(0,0,t,dl)
    lhs = numpy.array([pres, gi])
    rhs = numpy.array([pl, gl])
    
    pl_t = _eq_pressure(1,0,t,dl)
    pl_d = _eq_pressure(0,1,t,dl)
    gi_t = _ice_g(1,0,t,pres)
    gl_t = _eq_chempot(1,0,t,dl)
    gl_d = _eq_chempot(0,1,t,dl)
    dlhs = numpy.array([[0.,0.], [gi_t,0.]])
    drhs = numpy.array([[pl_t,pl_d], [gl_t,gl_d]])
    return lhs, rhs, dlhs, drhs

def eq_tp(temp=None,pres=None,dliq=None,chkvals=False,chktol=_CHKTOL,
    temp0=None,pres0=None,dliq0=None,chkbnd=False,mathargs=None):
    """Get primary ice-liquid variables at T or P.
    
    Get the values of all primary variables for ice and liquid water in
    equilibrium at either of a given temperature or pressure.
    
    If the calculation has already been done, the results can be passed
    to avoid unnecessary repeat calculations. If enough values are
    passed, they will be checked for consistency if chkvals is True.
    
    :arg temp: Temperature in K.
    :type temp: float or None
    :arg pres: Pressure in Pa.
    :type pres: float or None
    :arg dliq: Liquid water density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dliq: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg temp0: Initial guess for the temperature in K. If None
        (default) then `_approx_p` is used.
    :type temp0: float or None
    :arg pres0: Initial guess for the pressure in Pa. If None (default)
        then `_approx_t` is used.
    :type pres0: float or None
    :arg dliq0: Initial guess for the liquid water density in kg/m3. If
        None (default) then `_approx_t` or `_approx_p` is used.
    :type dliq0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Temperature, pressure, and liquid water density (all in SI
        units).
    :raises ValueError: If neither of temp or pres is provided.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    """
    if temp is None and pres is None:
        errmsg = 'One of temp or pres must be provided'
        raise ValueError(errmsg)
    if temp is not None:
        if any(val is None for val in (pres,dliq)):
            vals1 = (pres0,dliq0)
            vals2 = (pres0,dliq0)
            if any(val1 is None for val1 in vals1):
                vals2 = _approx_t(temp)
            x0 = numpy.array([val2 if val1 is None else val1
                for (val1,val2) in zip(vals1,vals2)])
            fargs = (temp,)
            if mathargs is None:
                mathargs = dict()
            x1 = _newton(_diff_t,x0,fargs=fargs,**mathargs)
            pres, dliq = x1
    else:
        vals1 = (temp0,dliq0)
        vals2 = (temp0,dliq0)
        if any(val1 is None for val1 in vals1):
            vals2 = _approx_p(pres)
        x0 = numpy.array([val2 if val1 is None else val1
            for (val1,val2) in zip(vals1,vals2)])
        fargs = (pres,)
        if mathargs is None:
            mathargs = dict()
        x1 = _newton(_diff_p,x0,fargs=fargs,**mathargs)
        temp, dliq = x1
    
    _chkflubnds(temp,dliq,chkbnd=chkbnd)
    _chkicebnds(temp,pres,chkbnd=chkbnd)
    if not chkvals:
        return temp, pres, dliq
    
    lhs, rhs, __, __ = _diff_p(temp,dliq,pres)
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
    return temp, pres, dliq


## Thermodynamic properties
def temperature(temp=None,pres=None,dliq=None,chkvals=False,
    chktol=_CHKTOL,temp0=None,pres0=None,dliq0=None,chkbnd=False,
    mathargs=None):
    """Calculate ice-liquid temperature.
    
    Calculate the temperature of ice and liquid water in equilibrium.
    
    :arg temp: Temperature in K.
    :type temp: float or None
    :arg pres: Pressure in Pa.
    :type pres: float or None
    :arg dliq: Liquid water density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dliq: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg temp0: Initial guess for the temperature in K. If None
        (default) then `_approx_p` is used.
    :type temp0: float or None
    :arg pres0: Initial guess for the pressure in Pa. If None (default)
        then `_approx_t` is used.
    :type pres0: float or None
    :arg dliq0: Initial guess for the liquid water density in kg/m3. If
        None (default) then `_approx_t` or `_approx_p` is used.
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
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> temperature(pres=1e7)
    272.401648869
    """
    temp, pres, dliq = eq_tp(temp=temp,pres=pres,dliq=dliq,chkvals=chkvals,
        chktol=chktol,temp0=temp0,pres0=pres0,dliq0=dliq0,chkbnd=chkbnd,
        mathargs=mathargs)
    return temp

def pressure(temp=None,pres=None,dliq=None,chkvals=False,
    chktol=_CHKTOL,temp0=None,pres0=None,dliq0=None,chkbnd=False,
    mathargs=None):
    """Calculate ice-liquid pressure.
    
    Calculate the pressure of ice and liquid water in equilibrium.
    
    :arg temp: Temperature in K.
    :type temp: float or None
    :arg pres: Pressure in Pa.
    :type pres: float or None
    :arg dliq: Liquid water density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dliq: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg temp0: Initial guess for the temperature in K. If None
        (default) then `_approx_p` is used.
    :type temp0: float or None
    :arg pres0: Initial guess for the pressure in Pa. If None (default)
        then `_approx_t` is used.
    :type pres0: float or None
    :arg dliq0: Initial guess for the liquid water density in kg/m3. If
        None (default) then `_approx_t` or `_approx_p` is used.
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
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> pressure(temp=270.)
    39313338.8825
    """
    temp, pres, dliq = eq_tp(temp=temp,pres=pres,dliq=dliq,chkvals=chkvals,
        chktol=chktol,temp0=temp0,pres0=pres0,dliq0=dliq0,chkbnd=chkbnd,
        mathargs=mathargs)
    return pres

def densityliq(temp=None,pres=None,dliq=None,chkvals=False,
    chktol=_CHKTOL,temp0=None,pres0=None,dliq0=None,chkbnd=False,
    mathargs=None):
    """Calculate ice-liquid liquid water density.
    
    Calculate the density of liquid water for ice and liquid water in
    equilibrium.
    
    :arg temp: Temperature in K.
    :type temp: float or None
    :arg pres: Pressure in Pa.
    :type pres: float or None
    :arg dliq: Liquid water density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dliq: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg temp0: Initial guess for the temperature in K. If None
        (default) then `_approx_p` is used.
    :type temp0: float or None
    :arg pres0: Initial guess for the pressure in Pa. If None (default)
        then `_approx_t` is used.
    :type pres0: float or None
    :arg dliq0: Initial guess for the liquid water density in kg/m3. If
        None (default) then `_approx_t` or `_approx_p` is used.
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
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> densityliq(pres=1e7)
    1004.79353660
    >>> densityliq(temp=270.)
    1019.05568894
    """
    temp, pres, dliq = eq_tp(temp=temp,pres=pres,dliq=dliq,chkvals=chkvals,
        chktol=chktol,temp0=temp0,pres0=pres0,dliq0=dliq0,chkbnd=chkbnd,
        mathargs=mathargs)
    return dliq

def chempot(temp=None,pres=None,dliq=None,chkvals=False,chktol=_CHKTOL,
    temp0=None,pres0=None,dliq0=None,chkbnd=False,mathargs=None):
    """Calculate ice-liquid chemical potential.
    
    Calculate the chemical potential of ice and liquid water in
    equilibrium.
    
    :arg temp: Temperature in K.
    :type temp: float or None
    :arg pres: Pressure in Pa.
    :type pres: float or None
    :arg dliq: Liquid water density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dliq: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg temp0: Initial guess for the temperature in K. If None
        (default) then `_approx_p` is used.
    :type temp0: float or None
    :arg pres0: Initial guess for the pressure in Pa. If None (default)
        then `_approx_t` is used.
    :type pres0: float or None
    :arg dliq0: Initial guess for the liquid water density in kg/m3. If
        None (default) then `_approx_t` or `_approx_p` is used.
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
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> chempot(pres=1e7)
    9972.88171381
    >>> chempot(temp=270.)
    38870.0605192
    """
    temp, pres, dliq = eq_tp(temp=temp,pres=pres,dliq=dliq,chkvals=chkvals,
        chktol=chktol,temp0=temp0,pres0=pres0,dliq0=dliq0,chkbnd=chkbnd,
        mathargs=mathargs)
    g = _ice_g(0,0,temp,pres)
    return g

def densityice(temp=None,pres=None,dliq=None,chkvals=False,
    chktol=_CHKTOL,temp0=None,pres0=None,dliq0=None,chkbnd=False,
    mathargs=None):
    """Calculate ice-liquid ice density.
    
    Calculate the density of ice for ice and liquid water in
    equilibrium.
    
    :arg temp: Temperature in K.
    :type temp: float or None
    :arg pres: Pressure in Pa.
    :type pres: float or None
    :arg dliq: Liquid water density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dliq: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg temp0: Initial guess for the temperature in K. If None
        (default) then `_approx_p` is used.
    :type temp0: float or None
    :arg pres0: Initial guess for the pressure in Pa. If None (default)
        then `_approx_t` is used.
    :type pres0: float or None
    :arg dliq0: Initial guess for the liquid water density in kg/m3. If
        None (default) then `_approx_t` or `_approx_p` is used.
    :type dliq0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Ice density in kg/m3.
    :raises ValueError: If neither of temp or pres is provided.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> densityice(pres=1e7)
    917.896690831
    >>> densityice(temp=270.)
    921.359428514
    """
    temp, pres, dliq = eq_tp(temp=temp,pres=pres,dliq=dliq,chkvals=chkvals,
        chktol=chktol,temp0=temp0,pres0=pres0,dliq0=dliq0,chkbnd=chkbnd,
        mathargs=mathargs)
    dice = ice2.density(temp,pres)
    return dice

def enthalpyice(temp=None,pres=None,dliq=None,chkvals=False,
    chktol=_CHKTOL,temp0=None,pres0=None,dliq0=None,chkbnd=False,
    mathargs=None):
    """Calculate liquid-ice ice enthalpy.
    
    Calculate the specific enthalpy of ice for ice and liquid water in
    equilibrium.
    
    :arg temp: Temperature in K.
    :type temp: float or None
    :arg pres: Pressure in Pa.
    :type pres: float or None
    :arg dliq: Liquid water density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dliq: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg temp0: Initial guess for the temperature in K. If None
        (default) then `_approx_p` is used.
    :type temp0: float or None
    :arg pres0: Initial guess for the pressure in Pa. If None (default)
        then `_approx_t` is used.
    :type pres0: float or None
    :arg dliq0: Initial guess for the liquid water density in kg/m3. If
        None (default) then `_approx_t` or `_approx_p` is used.
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
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> enthalpyice(pres=1e7)
    -324602.983812
    >>> enthalpyice(temp=270.)
    -299055.938629
    """
    temp, pres, dliq = eq_tp(temp=temp,pres=pres,dliq=dliq,chkvals=chkvals,
        chktol=chktol,temp0=temp0,pres0=pres0,dliq0=dliq0,chkbnd=chkbnd,
        mathargs=mathargs)
    hi = ice2.enthalpy(temp,pres)
    return hi

def enthalpyliq(temp=None,pres=None,dliq=None,chkvals=False,
    chktol=_CHKTOL,temp0=None,pres0=None,dliq0=None,chkbnd=False,
    mathargs=None):
    """Calculate ice-liquid liquid water enthalpy.
    
    Calculate the specific enthalpy of liquid water for ice and liquid
    water in equilibrium.
    
    :arg temp: Temperature in K.
    :type temp: float or None
    :arg pres: Pressure in Pa.
    :type pres: float or None
    :arg dliq: Liquid water density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dliq: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg temp0: Initial guess for the temperature in K. If None
        (default) then `_approx_p` is used.
    :type temp0: float or None
    :arg pres0: Initial guess for the pressure in Pa. If None (default)
        then `_approx_t` is used.
    :type pres0: float or None
    :arg dliq0: Initial guess for the liquid water density in kg/m3. If
        None (default) then `_approx_t` or `_approx_p` is used.
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
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> enthalpyliq(pres=1e7)
    6945.92700483
    >>> enthalpyliq(temp=270.)
    26110.7481094
    """
    temp, pres, dliq = eq_tp(temp=temp,pres=pres,dliq=dliq,chkvals=chkvals,
        chktol=chktol,temp0=temp0,pres0=pres0,dliq0=dliq0,chkbnd=chkbnd,
        mathargs=mathargs)
    hl = flu2.enthalpy(temp,dliq)
    return hl

def enthalpymelt(temp=None,pres=None,dliq=None,chkvals=False,
    chktol=_CHKTOL,temp0=None,pres0=None,dliq0=None,chkbnd=False,
    mathargs=None):
    """Calculate enthalpy of melting.
    
    Calculate the specific enthalpy of melting.
    
    :arg temp: Temperature in K.
    :type temp: float or None
    :arg pres: Pressure in Pa.
    :type pres: float or None
    :arg dliq: Liquid water density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dliq: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg temp0: Initial guess for the temperature in K. If None
        (default) then `_approx_p` is used.
    :type temp0: float or None
    :arg pres0: Initial guess for the pressure in Pa. If None (default)
        then `_approx_t` is used.
    :type pres0: float or None
    :arg dliq0: Initial guess for the liquid water density in kg/m3. If
        None (default) then `_approx_t` or `_approx_p` is used.
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
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> enthalpymelt(pres=1e7)
    331548.910817
    >>> enthalpymelt(temp=270.)
    325166.686739
    """
    temp, pres, dliq = eq_tp(temp=temp,pres=pres,dliq=dliq,chkvals=chkvals,
        chktol=chktol,temp0=temp0,pres0=pres0,dliq0=dliq0,chkbnd=chkbnd,
        mathargs=mathargs)
    hl = flu2.enthalpy(temp,dliq)
    hi = ice2.enthalpy(temp,pres)
    hmelt = hl - hi
    return hmelt

def entropyice(temp=None,pres=None,dliq=None,chkvals=False,
    chktol=_CHKTOL,temp0=None,pres0=None,dliq0=None,chkbnd=False,
    mathargs=None):
    """Calculate ice-liquid ice entropy.
    
    Calculate the specific entropy of ice for ice and liquid water in
    equilibrium.
    
    :arg temp: Temperature in K.
    :type temp: float or None
    :arg pres: Pressure in Pa.
    :type pres: float or None
    :arg dliq: Liquid water density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dliq: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg temp0: Initial guess for the temperature in K. If None
        (default) then `_approx_p` is used.
    :type temp0: float or None
    :arg pres0: Initial guess for the pressure in Pa. If None (default)
        then `_approx_t` is used.
    :type pres0: float or None
    :arg dliq0: Initial guess for the liquid water density in kg/m3. If
        None (default) then `_approx_t` or `_approx_p` is used.
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
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> entropyice(pres=1e7)
    -1228.24464139
    >>> entropyice(temp=270.)
    -1251.57777462
    """
    temp, pres, dliq = eq_tp(temp=temp,pres=pres,dliq=dliq,chkvals=chkvals,
        chktol=chktol,temp0=temp0,pres0=pres0,dliq0=dliq0,chkbnd=chkbnd,
        mathargs=mathargs)
    si = ice2.entropy(temp,pres)
    return si

def entropyliq(temp=None,pres=None,dliq=None,chkvals=False,
    chktol=_CHKTOL,temp0=None,pres0=None,dliq0=None,chkbnd=False,
    mathargs=None):
    """Calculate ice-liquid liquid entropy.
    
    Calculate the specific entropy of liquid water for ice and liquid
    water in equilibrium.
    
    :arg temp: Temperature in K.
    :type temp: float or None
    :arg pres: Pressure in Pa.
    :type pres: float or None
    :arg dliq: Liquid water density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dliq: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg temp0: Initial guess for the temperature in K. If None
        (default) then `_approx_p` is used.
    :type temp0: float or None
    :arg pres0: Initial guess for the pressure in Pa. If None (default)
        then `_approx_t` is used.
    :type pres0: float or None
    :arg dliq0: Initial guess for the liquid water density in kg/m3. If
        None (default) then `_approx_t` or `_approx_p` is used.
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
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> entropyliq(pres=1e7)
    -11.1121012723
    >>> entropyliq(temp=270.)
    -47.2567126291
    """
    temp, pres, dliq = eq_tp(temp=temp,pres=pres,dliq=dliq,chkvals=chkvals,
        chktol=chktol,temp0=temp0,pres0=pres0,dliq0=dliq0,chkbnd=chkbnd,
        mathargs=mathargs)
    sl = flu2.entropy(temp,dliq)
    return sl

def entropymelt(temp=None,pres=None,dliq=None,chkvals=False,
    chktol=_CHKTOL,temp0=None,pres0=None,dliq0=None,chkbnd=False,
    mathargs=None):
    """Calculate entropy of melting.
    
    Calculate the specific entropy of melting.
    
    :arg temp: Temperature in K.
    :type temp: float or None
    :arg pres: Pressure in Pa.
    :type pres: float or None
    :arg dliq: Liquid water density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dliq: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg temp0: Initial guess for the temperature in K. If None
        (default) then `_approx_p` is used.
    :type temp0: float or None
    :arg pres0: Initial guess for the pressure in Pa. If None (default)
        then `_approx_t` is used.
    :type pres0: float or None
    :arg dliq0: Initial guess for the liquid water density in kg/m3. If
        None (default) then `_approx_t` or `_approx_p` is used.
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
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> entropymelt(pres=1e7)
    1217.13254011
    >>> entropymelt(temp=270.)
    1204.32106199
    """
    temp, pres, dliq = eq_tp(temp=temp,pres=pres,dliq=dliq,chkvals=chkvals,
        chktol=chktol,temp0=temp0,pres0=pres0,dliq0=dliq0,chkbnd=chkbnd,
        mathargs=mathargs)
    sl = flu2.entropy(temp,dliq)
    si = ice2.entropy(temp,pres)
    smelt = sl - si
    return smelt

def volumemelt(temp=None,pres=None,dliq=None,chkvals=False,
    chktol=_CHKTOL,temp0=None,pres0=None,dliq0=None,chkbnd=False,
    mathargs=None):
    """Calculate specific volume of melting.
    
    Calculate the specific volume of melting.
    
    :arg temp: Temperature in K.
    :type temp: float or None
    :arg pres: Pressure in Pa.
    :type pres: float or None
    :arg dliq: Liquid water density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dliq: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg temp0: Initial guess for the temperature in K. If None
        (default) then `_approx_p` is used.
    :type temp0: float or None
    :arg pres0: Initial guess for the pressure in Pa. If None (default)
        then `_approx_t` is used.
    :type pres0: float or None
    :arg dliq0: Initial guess for the liquid water density in kg/m3. If
        None (default) then `_approx_t` or `_approx_p` is used.
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
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> volumemelt(pres=1e7)
    -9.42178903288e-5
    >>> volumemelt(temp=270.)
    -1.04052121182e-4
    """
    temp, pres, dliq = eq_tp(temp=temp,pres=pres,dliq=dliq,chkvals=chkvals,
        chktol=chktol,temp0=temp0,pres0=pres0,dliq0=dliq0,chkbnd=chkbnd,
        mathargs=mathargs)
    vi = _ice_g(0,1,temp,pres,chkbnd=chkbnd)
    vl = dliq**(-1)
    vmelt = vl - vi
    return vmelt

