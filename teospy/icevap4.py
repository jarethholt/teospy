"""Ice-water vapour equilibrium functions.

This module provides thermodynamic properties of ice and water vapour in
equilibrium, e.g. the enthalpy of sublimation.

:Examples:

>>> temperature(pres=100.)
252.817910215
>>> densityvap(pres=100.)
8.57185487853e-4
>>> volumesubl(pres=100.)
1166.60755699
>>> entropysubl(pres=100.)
11225.8717816
>>> enthalpysubl(pres=100.)
2838101.44416
>>> pressure(temp=270.)
470.059067981
>>> densityvap(temp=270.)
3.77406140772e-3
>>> volumesubl(temp=270.)
264.965451558
>>> entropysubl(temp=270.)
10500.6135349
>>> enthalpysubl(temp=270.)
2835165.65442

:Functions:

* :func:`eq_tp`: Calculate ice-water vapour equilibrium properties at
  either temperature or pressure.
* :func:`temperature`: Temperature at ice-water vapour equilibrium.
* :func:`pressure`: Pressure at ice-water vapour equilibrium.
* :func:`densityvap`: Water vapour density at ice-water vapour
  equilibrium.
* :func:`chempot`: Chemical potential at ice-water vapour equilibrium.
* :func:`densityice`: Ice density at ice-water vapour equilibrium.
* :func:`enthalpyice`: Ice enthalpy at ice-water vapour equilibrium.
* :func:`enthalpyvap`: Water vapour enthalpy at ice-water vapour
  equilibrium.
* :func:`entropyice`: Ice entropy at ice-water vapour equilibrium.
* :func:`entropyvap`: Water vapour entropy at ice-water vapour
  equilibrium.
* :func:`volumesubl`: Specific volume of sublimation.
* :func:`entropysubl`: Specific entropy of sublimation.
* :func:`enthalpysubl`: Specific enthalpy of sublimation.

"""

__all__ = ['eq_tp','temperature','pressure','densityvap','chempot','densityice',
    'enthalpyice','enthalpyvap','entropyice','entropyvap','volumesubl',
    'entropysubl','enthalpysubl']

import warnings
import numpy
from teospy import constants0
from teospy import ice1
from teospy import flu2
from teospy import ice2
from teospy import maths3
from teospy import maths4

_CHKTOL = constants0.CHKTOL
_RWAT = constants0.RWAT
_TTP = constants0.TTP
_PTPE = constants0.PTPE
_LLVTP = constants0.LLVTP
_LILTP = constants0.LILTP
_CICE = constants0.CICE
_CVAP = constants0.CVAP
_chkflubnds = constants0.chkflubnds
_chkicebnds = constants0.chkicebnds
_ice_g = ice1.ice_g
_eq_chempot = flu2.eq_chempot
_eq_pressure = flu2.eq_pressure
_newton = maths3.newton
_AVI = (_LLVTP+_LILTP)/(_RWAT*_TTP)
_BVI = (_CICE-_CVAP)/_RWAT
_RAB = _AVI/_BVI


## Equilibrium functions
def _approx_t(temp):
    """Approximate PDv at T.
    
    Approximate the pressure and water vapour density of ice and water
    vapour in equilibrium at the given temperature. This approximation
    is based on constant heat capacities.
    
    :arg float temp: Temperature in K.
    :returns: Pressure in Pa and water vapour density in kg/m3.
    """
    earg = _AVI * (1 - _TTP/temp)
    earg += _BVI * (1 - _TTP/temp - numpy.log(temp/_TTP))
    pres = _PTPE * numpy.exp(earg)
    dvap = pres / (_RWAT * temp)
    return pres, dvap

def _approx_p(pres):
    """Approximate TDv at P.
    
    Approximate the temperature and water vapour density of ice and
    water vapour in equilibrium at the given pressure. This
    approximation is based on constant heat capacities.
    
    :arg float pres: Pressure in Pa.
    :returns: Temperature in K and water vapour density in kg/m3.
    """
    v = numpy.log(pres/_PTPE)/_BVI
    x = maths4.lamb2(v,_RAB)
    temp = _TTP/x
    dvap = pres / (_RWAT * temp)
    return temp, dvap

def _diff_t(p,dv,temp):
    """Calculate ice-vapour disequilibrium at T.
    
    Calculate both sides of the equations
    
        given pressure = pressure of water vapour
        chemical potential of ice = potential of water vapour
    
    and their Jacobians with respect to pressure and water vapour
    density. Solving these equations gives the pressure and water vapour
    density at the given temperature.
    
    :arg float p: Pressure in Pa.
    :arg float dv: Water vapour density in kg/m3.
    :arg float temp: Temperature in K.
    :returns: Left-hand side of the equation, right-hand side,
        Jacobian of LHS, and Jacobian of RHS.
    :rtype: tuple(array(float))
    """
    pv = _eq_pressure(0,0,temp,dv)
    gv = _eq_chempot(0,0,temp,dv)
    gi = _ice_g(0,0,temp,p)
    lhs = numpy.array([p, gi])
    rhs = numpy.array([pv, gv])
    
    pv_d = _eq_pressure(0,1,temp,dv)
    gi_p = _ice_g(0,1,temp,p)
    gv_d = _eq_chempot(0,1,temp,dv)
    dlhs = numpy.array([[1.,0.], [gi_p,0.]])
    drhs = numpy.array([[0.,pv_d], [0.,gv_d]])
    return lhs, rhs, dlhs, drhs

def _diff_p(t,dv,pres):
    """Calculate ice-vapour disequilibrium at P.
    
    Calculate both sides of the equations
    
        given pressure = pressure of water vapour
        chemical potential of ice = potential of water vapour
    
    and their Jacobians with respect to temperature and water vapour
    density. Solving these equations gives the temperature and water
    vapour density at the given pressure.
    
    :arg float t: Temperature in K.
    :arg float dv: Water vapour density in kg/m3.
    :arg float pres: Pressure in Pa.
    :returns: Left-hand side of the equation, right-hand side,
        Jacobian of LHS, and Jacobian of RHS.
    :rtype: tuple(array(float))
    """
    pv = _eq_pressure(0,0,t,dv)
    gv = _eq_chempot(0,0,t,dv)
    gi = _ice_g(0,0,t,pres)
    lhs = numpy.array([pres, gi])
    rhs = numpy.array([pv, gv])
    
    pv_t = _eq_pressure(1,0,t,dv)
    pv_d = _eq_pressure(0,1,t,dv)
    gv_t = _eq_chempot(1,0,t,dv)
    gi_t = _ice_g(1,0,t,pres)
    gv_d = _eq_chempot(0,1,t,dv)
    dlhs = numpy.array([[0.,0.], [gi_t,0.]])
    drhs = numpy.array([[pv_t,pv_d], [gv_t,gv_d]])
    return lhs, rhs, dlhs, drhs

def eq_tp(temp=None,pres=None,dvap=None,chkvals=False,chktol=_CHKTOL,
    temp0=None,pres0=None,dvap0=None,chkbnd=False,mathargs=None):
    """Get primary ice-vapour variables at T or P.
    
    Get the values of all primary variables for ice and water vapour in
    equilibrium at either of a given temperature or pressure.
    
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
    :arg dvap0: Initial guess for the water vapour density in kg/m3. If
        None (default) then `_approx_t` or `_approx_p` is used.
    :type dvap0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Temperature, pressure, and water vapour density (all in SI
        units).
    :raises ValueError: If neither of temp or pres is provided.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    """
    if temp is None and pres is None:
        errmsg = 'One of temp or pres must be provided'
        raise ValueError(errmsg)
    if temp is not None:
        if any(val is None for val in (pres,dvap)):
            x0 = (pres0,dvap0)
            fargs = (temp,)
            if mathargs is None:
                mathargs = dict()
            x1 = _newton(_diff_t,x0,_approx_t,fargs=fargs,**mathargs)
            pres, dvap = x1
    else:
        x0 = (temp0,dvap0)
        fargs = (pres,)
        if mathargs is None:
            mathargs = dict()
        x1 = _newton(_diff_p,x0,_approx_p,fargs=fargs,**mathargs)
        temp, dvap = x1
    
    _chkflubnds(temp,dvap,chkbnd=chkbnd)
    _chkicebnds(temp,pres,chkbnd=chkbnd)
    if not chkvals:
        return temp, pres, dvap
    
    lhs, rhs, __, __ = _diff_p(temp,dvap,pres)
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
    return temp, pres, dvap


### Thermodynamic properties
def temperature(temp=None,pres=None,dvap=None,chkvals=False,
    chktol=_CHKTOL,temp0=None,pres0=None,dvap0=None,chkbnd=False,
    mathargs=None):
    """Calculate ice-vapour temperature.
    
    Calculate the temperature of ice and water vapour in equilibrium.
    
    :arg temp: Temperature in K.
    :type temp: float or None
    :arg pres: Pressure in Pa.
    :type pres: float or None
    :arg dvap: Water vapour density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dvap: float or None
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
    :arg dvap0: Initial guess for the water vapour density in kg/m3. If
        None (default) then `_approx_t` or `_approx_p` is used.
    :type dvap0: float or None
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
    
    >>> temperature(pres=100.)
    252.817910215
    """
    temp, pres, dvap = eq_tp(temp=temp,pres=pres,dvap=dvap,chkvals=chkvals,
        chktol=chktol,temp0=temp0,pres0=pres0,dvap0=dvap0,chkbnd=chkbnd,
        mathargs=mathargs)
    return temp

def pressure(temp=None,pres=None,dvap=None,chkvals=False,chktol=_CHKTOL,
    temp0=None,pres0=None,dvap0=None,chkbnd=False,mathargs=None):
    """Calculate ice-vapour pressure.
    
    Calculate the pressure of ice and water vapour in equilibrium.
    
    :arg temp: Temperature in K.
    :type temp: float or None
    :arg pres: Pressure in Pa.
    :type pres: float or None
    :arg dvap: Water vapour density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dvap: float or None
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
    :arg dvap0: Initial guess for the water vapour density in kg/m3. If
        None (default) then `_approx_t` or `_approx_p` is used.
    :type dvap0: float or None
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
    470.059067981
    """
    temp, pres, dvap = eq_tp(temp=temp,pres=pres,dvap=dvap,chkvals=chkvals,
        chktol=chktol,temp0=temp0,pres0=pres0,dvap0=dvap0,chkbnd=chkbnd,
        mathargs=mathargs)
    return pres

def densityvap(temp=None,pres=None,dvap=None,chkvals=False,
    chktol=_CHKTOL,temp0=None,pres0=None,dvap0=None,chkbnd=False,
    mathargs=None):
    """Calculate ice-vapour vapour density.
    
    Calculate the density of water vapour for ice and water vapour in
    equilibrium.
    
    :arg temp: Temperature in K.
    :type temp: float or None
    :arg pres: Pressure in Pa.
    :type pres: float or None
    :arg dvap: Water vapour density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dvap: float or None
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
    :arg dvap0: Initial guess for the water vapour density in kg/m3. If
        None (default) then `_approx_t` or `_approx_p` is used.
    :type dvap0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Water vapour density in kg/m3.
    :raises ValueError: If neither of temp or pres is provided.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> densityvap(temp=270.)
    3.77406140772e-3
    >>> densityvap(pres=100.)
    8.57185487853e-4
    """
    temp, pres, dvap = eq_tp(temp=temp,pres=pres,dvap=dvap,chkvals=chkvals,
        chktol=chktol,temp0=temp0,pres0=pres0,dvap0=dvap0,chkbnd=chkbnd,
        mathargs=mathargs)
    return dvap

def chempot(temp=None,pres=None,dvap=None,chkvals=False,chktol=_CHKTOL,
    temp0=None,pres0=None,dvap0=None,chkbnd=False,mathargs=None):
    """Calculate ice-vapour chemical potential.
    
    Calculate the chemical potential of ice and water vapour in
    equilibrium.
    
    :arg temp: Temperature in K.
    :type temp: float or None
    :arg pres: Pressure in Pa.
    :type pres: float or None
    :arg dvap: Water vapour density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dvap: float or None
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
    :arg dvap0: Initial guess for the water vapour density in kg/m3. If
        None (default) then `_approx_t` or `_approx_p` is used.
    :type dvap0: float or None
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
    
    >>> chempot(temp=270.)
    -3895.26747392
    >>> chempot(pres=100.)
    -26421.2820403
    """
    temp, pres, dvap = eq_tp(temp=temp,pres=pres,dvap=dvap,chkvals=chkvals,
        chktol=chktol,temp0=temp0,pres0=pres0,dvap0=dvap0,chkbnd=chkbnd,
        mathargs=mathargs)
    g = _ice_g(0,0,temp,pres)
    return g

def densityice(temp=None,pres=None,dvap=None,chkvals=False,
    chktol=_CHKTOL,temp0=None,pres0=None,dvap0=None,chkbnd=False,
    mathargs=None):
    """Calculate ice-vapour ice density.
    
    Calculate the density of ice for ice and water vapour in
    equilibrium.
    
    :arg temp: Temperature in K.
    :type temp: float or None
    :arg pres: Pressure in Pa.
    :type pres: float or None
    :arg dvap: Water vapour density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dvap: float or None
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
    :arg dvap0: Initial guess for the water vapour density in kg/m3. If
        None (default) then `_approx_t` or `_approx_p` is used.
    :type dvap0: float or None
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
    
    >>> densityice(temp=270.)
    917.170465733
    >>> densityice(pres=100.)
    919.600269745
    """
    temp, pres, dvap = eq_tp(temp=temp,pres=pres,dvap=dvap,chkvals=chkvals,
        chktol=chktol,temp0=temp0,pres0=pres0,dvap0=dvap0,chkbnd=chkbnd,
        mathargs=mathargs)
    dice = ice2.density(temp,pres)
    return dice

def enthalpyice(temp=None,pres=None,dvap=None,chkvals=False,
    chktol=_CHKTOL,temp0=None,pres0=None,dvap0=None,chkbnd=False,
    mathargs=None):
    """Calculate ice-vapour ice enthalpy.
    
    Calculate the specific enthalpy of ice for ice and water vapour in
    equilibrium.
    
    :arg temp: Temperature in K.
    :type temp: float or None
    :arg pres: Pressure in Pa.
    :type pres: float or None
    :arg dvap: Water vapour density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dvap: float or None
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
    :arg dvap0: Initial guess for the water vapour density in kg/m3. If
        None (default) then `_approx_t` or `_approx_p` is used.
    :type dvap0: float or None
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
    
    >>> enthalpyice(temp=270.)
    -340033.434649
    >>> enthalpyice(pres=100.)
    -374576.247867
    """
    temp, pres, dvap = eq_tp(temp=temp,pres=pres,dvap=dvap,chkvals=chkvals,
        chktol=chktol,temp0=temp0,pres0=pres0,dvap0=dvap0,chkbnd=chkbnd,
        mathargs=mathargs)
    hi = ice2.enthalpy(temp,pres)
    return hi

def enthalpyvap(temp=None,pres=None,dvap=None,chkvals=False,
    chktol=_CHKTOL,temp0=None,pres0=None,dvap0=None,chkbnd=False,
    mathargs=None):
    """Calculate ice-vapour vapour enthalpy.
    
    Calculate the specific enthalpy of water vapour for ice and water
    vapour in equilibrium.
    
    :arg temp: Temperature in K.
    :type temp: float or None
    :arg pres: Pressure in Pa.
    :type pres: float or None
    :arg dvap: Water vapour density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dvap: float or None
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
    :arg dvap0: Initial guess for the water vapour density in kg/m3. If
        None (default) then `_approx_t` or `_approx_p` is used.
    :type dvap0: float or None
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
    
    >>> enthalpyvap(temp=270.)
    2495132.21977
    >>> enthalpyvap(pres=100.)
    2463525.19629
    """
    temp, pres, dvap = eq_tp(temp=temp,pres=pres,dvap=dvap,chkvals=chkvals,
        chktol=chktol,temp0=temp0,pres0=pres0,dvap0=dvap0,chkbnd=chkbnd,
        mathargs=mathargs)
    hv = flu2.enthalpy(temp,dvap)
    return hv

def entropyice(temp=None,pres=None,dvap=None,chkvals=False,
    chktol=_CHKTOL,temp0=None,pres0=None,dvap0=None,chkbnd=False,
    mathargs=None):
    """Calculate ice-vapour ice entropy.
    
    Calculate the specific entropy of ice for ice and water vapour in
    equilibrium.
    
    :arg temp: Temperature in K.
    :type temp: float or None
    :arg pres: Pressure in Pa.
    :type pres: float or None
    :arg dvap: Water vapour density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dvap: float or None
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
    :arg dvap0: Initial guess for the water vapour density in kg/m3. If
        None (default) then `_approx_t` or `_approx_p` is used.
    :type dvap0: float or None
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
    
    >>> entropyice(temp=270.)
    -1244.95617472
    >>> entropyice(pres=100.)
    -1377.09771247
    """
    temp, pres, dvap = eq_tp(temp=temp,pres=pres,dvap=dvap,chkvals=chkvals,
        chktol=chktol,temp0=temp0,pres0=pres0,dvap0=dvap0,chkbnd=chkbnd,
        mathargs=mathargs)
    si = ice2.entropy(temp,pres)
    return si

def entropyvap(temp=None,pres=None,dvap=None,chkvals=False,
    chktol=_CHKTOL,temp0=None,pres0=None,dvap0=None,chkbnd=False,
    mathargs=None):
    """Calculate ice-vapour vapour entropy.
    
    Calculate the specific entropy of water vapour for ice and water
    vapour in equilibrium.
    
    :arg temp: Temperature in K.
    :type temp: float or None
    :arg pres: Pressure in Pa.
    :type pres: float or None
    :arg dvap: Water vapour density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dvap: float or None
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
    :arg dvap0: Initial guess for the water vapour density in kg/m3. If
        None (default) then `_approx_t` or `_approx_p` is used.
    :type dvap0: float or None
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
    
    >>> entropyvap(temp=270.)
    9255.65736018
    >>> entropyvap(pres=100.)
    9848.77406912
    """
    temp, pres, dvap = eq_tp(temp=temp,pres=pres,dvap=dvap,chkvals=chkvals,
        chktol=chktol,temp0=temp0,pres0=pres0,dvap0=dvap0,chkbnd=chkbnd,
        mathargs=mathargs)
    sv = flu2.entropy(temp,dvap)
    return sv

def volumesubl(temp=None,pres=None,dvap=None,chkvals=False,
    chktol=_CHKTOL,temp0=None,pres0=None,dvap0=None,chkbnd=False,
    mathargs=None):
    """Calculate volume of sublimation.
    
    Calculate the specific volume of sublimation.
    
    :arg temp: Temperature in K.
    :type temp: float or None
    :arg pres: Pressure in Pa.
    :type pres: float or None
    :arg dvap: Water vapour density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dvap: float or None
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
    :arg dvap0: Initial guess for the water vapour density in kg/m3. If
        None (default) then `_approx_t` or `_approx_p` is used.
    :type dvap0: float or None
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
    
    >>> volumesubl(temp=270.)
    264.965451558
    >>> volumesubl(pres=100.)
    1166.60755699
    """
    temp, pres, dvap = eq_tp(temp=temp,pres=pres,dvap=dvap,chkvals=chkvals,
        chktol=chktol,temp0=temp0,pres0=pres0,dvap0=dvap0,chkbnd=chkbnd,
        mathargs=mathargs)
    vv = dvap**(-1)
    vi = _ice_g(0,1,temp,pres)
    vsubl = vv - vi
    return vsubl

def entropysubl(temp=None,pres=None,dvap=None,chkvals=False,
    chktol=_CHKTOL,temp0=None,pres0=None,dvap0=None,chkbnd=False,
    mathargs=None):
    """Calculate entropy of sublimation.
    
    Calculate the specific entropy of sublimation.
    
    :arg temp: Temperature in K.
    :type temp: float or None
    :arg pres: Pressure in Pa.
    :type pres: float or None
    :arg dvap: Water vapour density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dvap: float or None
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
    :arg dvap0: Initial guess for the water vapour density in kg/m3. If
        None (default) then `_approx_t` or `_approx_p` is used.
    :type dvap0: float or None
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
    
    >>> entropysubl(temp=270.)
    10500.6135349
    >>> entropysubl(pres=100.)
    11225.8717816
    """
    temp, pres, dvap = eq_tp(temp=temp,pres=pres,dvap=dvap,chkvals=chkvals,
        chktol=chktol,temp0=temp0,pres0=pres0,dvap0=dvap0,chkbnd=chkbnd,
        mathargs=mathargs)
    sv = flu2.entropy(temp,dvap)
    si = ice2.entropy(temp,pres)
    ssubl = sv - si
    return ssubl

def enthalpysubl(temp=None,pres=None,dvap=None,chkvals=False,
    chktol=_CHKTOL,temp0=None,pres0=None,dvap0=None,chkbnd=False,
    mathargs=None):
    """Calculate enthalpy of sublimation.
    
    Calculate the specific enthalpy of sublimation.
    
    :arg temp: Temperature in K.
    :type temp: float or None
    :arg pres: Pressure in Pa.
    :type pres: float or None
    :arg dvap: Water vapour density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dvap: float or None
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
    :arg dvap0: Initial guess for the water vapour density in kg/m3. If
        None (default) then `_approx_t` or `_approx_p` is used.
    :type dvap0: float or None
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
    
    >>> enthalpysubl(temp=270.)
    2835165.65442
    >>> enthalpysubl(pres=100.)
    2838101.44416
    """
    temp, pres, dvap = eq_tp(temp=temp,pres=pres,dvap=dvap,chkvals=chkvals,
        chktol=chktol,temp0=temp0,pres0=pres0,dvap0=dvap0,chkbnd=chkbnd,
        mathargs=mathargs)
    hv = flu2.enthalpy(temp,dvap)
    hi = ice2.enthalpy(temp,pres)
    hsubl = hv - hi
    return hsubl

