"""Seawater salinity function.

This module provides the salinity of seawater (liquid water and salt) as
a function of temperature, pressure, and density.

:Examples:

>>> salinity(273.15,101325.,1028.)
3.50315257709e-2

:Functions:

* :func:`salinity`: Seawater salinity.

"""

__all__ = ['salinity']

import warnings
import numpy
import constants0
import flu2
import sal2
import maths3
import flu3a

_CHKTOL = constants0.CHKTOL
_SAL0 = constants0.SAL0
_chkflubnds = constants0.chkflubnds
_chksalbnds = constants0.chksalbnds
_flu_eq_pressure = flu2.eq_pressure
_sal_g = sal2.sal_g
_newton = maths3.newton


## Salinity functions
def _approx_tpd(temp,pres,dsea):
    """Approximate SDl at TPD.
    
    Approximate the salinity and liquid water density of seawater for
    the given temperature, pressure, and seawater density. Uses
    properties at the reference salinity _SAL0.
    
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg float dsea: Seawater density in kg/m3.
    :returns: Temperature and liquid water density (in SI units).
    """
    dliq = flu3a._dliq_default(temp,pres)
    lhs = dsea**(-1) - dliq**(-1)
    lhs -= _sal_g(0,0,1,_SAL0,temp,pres)
    rhs = _sal_g(1,0,1,_SAL0,temp,pres)
    salt = _SAL0 + lhs/rhs
    salt = max(min(salt,1.),0.)  # Only 0 <= salt <= 1 is valid
    return salt, dliq

def _diff_tpd(s,dl,temp,pres,dsea,useext=False):
    """Calculate seawater disequilibrium at TPD.
    
    Calculate both sides of the equations
    
        given pressure = pressure in fluid water
        1/(given density) = specific volume of seawater
    
    and their Jacobians with respect to salinity and liquid water
    density. Solving these equations gives the salinity and liquid water
    density for the given temperature, pressure, and seawater density.
    
    :arg float s: Salinity in kg/kg.
    :arg float dl: Liquid water density in kg/m3.
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg float dsea: Seawater density in kg/m3.
    :arg bool useext: If False (default) then the salt contribution is
        calculated from _GSCOEFFS in sal1; if True, from _GSCOEFFS_EXT.
    :returns: Left-hand side of the equation, right-hand side,
        Jacobian of LHS, and Jacobian of RHS.
    :rtype: tuple(array(float))
    """
    pl = _flu_eq_pressure(0,0,temp,dl)
    vl = dl**(-1)
    vs = _sal_g(0,0,1,s,temp,pres,useext=useext)
    vsea = dsea**(-1)
    lhs = numpy.array([pres, vsea])
    rhs = numpy.array([pl, vl+vs])
    
    pl_d = _flu_eq_pressure(0,1,temp,dl)
    vl_d = -dl**(-2)
    vs_s = _sal_g(1,0,1,s,temp,pres,useext=useext)
    dlhs = numpy.array([[0.,0.], [0.,0.]])
    drhs = numpy.array([[0.,pl_d], [vs_s,vl_d]])
    return lhs, rhs, dlhs, drhs

def eq_tpd(temp,pres,dsea,salt=None,dliq=None,chkvals=False,
    chktol=_CHKTOL,salt0=None,dliq0=None,chkbnd=False,useext=False,
    mathargs=None):
    """Get primary seawater variables at TPD.
    
    Get the values of all primary seawater variables for the given
    temperature, pressure, and seawater density.
    
    If the calculation has already been done, the results can be passed
    to avoid unnecessary repeat calculations. If enough values are
    passed, they will be checked for consistency if chkvals is True.
    
    :arg float temp: Temperature in K.
    :arg float pres: In-situ pressure in Pa.
    :arg float dsea: Seawater density in kg/m3.
    :arg salt: Salinity in kg/kg. If unknown, pass None (default) and it
        will be calculated.
    :type salt: float or None
    :arg dliq: Liquid water density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dliq: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg salt0: Initial guess for the salinity in kg/kg. If None
        (default) then `_approx_tpd` is used.
    :type salt0: float or None
    :arg dliq0: Initial guess for the liquid water density in kg/m3. If
        None (default) then `_approx_tpd` is used.
    :type dliq0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Salinity and liquid water density (in SI units).
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    """
    if any(val is None for val in (salt,dliq)):
        x0 = (salt0,dliq0)
        fargs = (temp,pres,dsea)
        fkwargs = {'useext': useext}
        if mathargs is None:
            mathargs = dict()
        x1 = _newton(_diff_tpd,x0,_approx_tpd,fargs=fargs,fkwargs=fkwargs,
            **mathargs)
        salt, dliq = x1
    
    _chkflubnds(temp,dliq,chkbnd=chkbnd)
    _chksalbnds(salt,temp,pres,chkbnd=chkbnd)
    if not chkvals:
        return salt, dliq
    
    lhs, rhs, __, __ = _diff_tpd(salt,dliq,temp,pres,dsea,useext=useext)
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
    return salt, dliq

def salinity(temp,pres,dsea,salt=None,dliq=None,chkvals=False,chktol=_CHKTOL,salt0=None,dliq0=None,chkbnd=False,useext=False,mathargs=None):
    """Calculate seawater salinity.
    
    Calculate the absolute salinity of seawater from the temperature,
    pressure, and seawater density.
    
    :arg float temp: Temperature in K.
    :arg float pres: In-situ pressure in Pa.
    :arg float dsea: Seawater density in kg/m3.
    :arg salt: Salinity in kg/kg. If unknown, pass None (default) and it
        will be calculated.
    :type salt: float or None
    :arg dliq: Liquid water density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dliq: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg salt0: Initial guess for the salinity in kg/kg. If None
        (default) then `_approx_tpd` is used.
    :type salt0: float or None
    :arg dliq0: Initial guess for the liquid water density in kg/m3. If
        None (default) then `_approx_tpd` is used.
    :type dliq0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Salinity in kg/kg.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> salinity(273.15,101325.,1028.)
    3.50315257709e-2
    """
    salt, dliq = eq_tpd(temp,pres,dsea,salt=salt,dliq=dliq,chkvals=chkvals,
        chktol=chktol,salt0=salt0,dliq0=dliq0,chkbnd=chkbnd,useext=useext,
        mathargs=mathargs)
    return salt

