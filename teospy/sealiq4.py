"""Seawater-pure liquid water equilibrium functions.

This module provides the osmotic pressure, the additional pressure
required to bring seawater to equilibrium with pure liquid water, for
example to prevent the flow of pure water across a semipermeable
membrane.

:Examples:

>>> osmoticpressure(0.035,300.,1e5)
2594603.20968

:Functions:

* :func:`eq_stp`: Calculate pure water-seawater equilibrium properties
  at salinity, temperature, and pure water pressure.
* :func:`osmoticpressure`: Seawater osmotic pressure.

"""

__all__ = ['eq_stp','osmoticpressure']

import warnings
import numpy
from teospy import constants0
from teospy import flu2
from teospy import sal2
from teospy import maths3
from teospy import flu3a

_CHKTOL = constants0.CHKTOL
_MSAL = constants0.MSAL
_RUNIV = constants0.RUNIV
_DLTP = constants0.DLTP
_PCOEF = _RUNIV/_MSAL*_DLTP
_chkflubnds = constants0.chkflubnds
_chksalbnds = constants0.chksalbnds
_eq_chempot = flu2.eq_chempot
_eq_pressure = flu2.eq_pressure
_eq_liqpot = sal2.eq_liqpot
_newton = maths3.newton
_dliq_default = flu3a._dliq_default


## Equilibrium functions
def _approx_stp(salt,temp,ppur,dlpur):
    """Approximate PsDls at STPp.
    
    Approximate the pressure in seawater and density of liquid water in
    seawater for seawater at the given salinity and temperature in
    equilibrium with pure liquid water at the given temperature and
    pressure.
    
    :arg float salt: Salinity in kg/kg.
    :arg float temp: Temperature in K.
    :arg float ppur: Pure water pressure in Pa.
    :arg float dlpur: Pure water density in kg/m3 (unused).
    :returns: Seawater pressure and density of liquid water in seawater
        (both in SI units).
    """
    psea = ppur + _PCOEF*temp*salt
    dlsea = _dliq_default(temp,psea)
    return psea, dlsea

def _diff_stp(ps,dls,salt,temp,ppur,dlpur,useext=False):
    """Calculate seawater-pure water disequilibrium at STPp.
    
    Calculate both sides of the equations
    
        given seawater pressure = pressure in seawater
        potential of pure water = potential of water in seawater
    
    and their Jacobians with respect to seawater pressure and density.
    Solving these equations gives equilibrium values at the given
    seawater salinity, temperature, and pure water pressure.
    
    :arg float ps: Seawater pressure in Pa.
    :arg float dls: Density of liquid water in seawater in kg/m3.
    :arg float salt: Salinity in kg/kg.
    :arg float temp: Temperature in K.
    :arg float ppur: Pure water pressure in Pa.
    :arg float dlpur: Pure water density in kg/m3.
    :arg bool useext: If False (default) then the salt contribution is
        calculated from _GSCOEFFS; if True, from _GSCOEFFS_EXT.
    :returns: Left-hand side of the equation, right-hand side,
        Jacobian of LHS, and Jacobian of RHS.
    :rtype: tuple(array(float))
    """
    pls = _eq_pressure(0,0,temp,dls)
    muls = _eq_chempot(0,0,temp,dls)
    muls += _eq_liqpot(0,0,0,salt,temp,ps,useext=useext)
    mulp = _eq_chempot(0,0,temp,dlpur)
    lhs = numpy.array([ps, mulp])
    rhs = numpy.array([pls, muls])
    
    pls_d = _eq_pressure(0,1,temp,dls)
    muls_p = _eq_liqpot(0,0,1,salt,temp,ps,useext=useext)
    muls_d = _eq_chempot(0,1,temp,dls)
    dlhs = numpy.array([[1.,0.], [0.,0.]])
    drhs = numpy.array([[0.,pls_d], [muls_p,muls_d]])
    return lhs, rhs, dlhs, drhs

def eq_stp(salt,temp,ppur,dlpur=None,psea=None,dlsea=None,chkvals=False,
    chktol=_CHKTOL,dlpur0=None,psea0=None,dlsea0=None,chkbnd=False,
    useext=False,mathargs=None):
    """Get primary seawater-pure water variables at STPp.
    
    Get the values of all primary variables for seawater at the given
    salinity and temperature in equilibrium with pure liquid water at
    the given temperature and pressure.
    
    If the calculation has already been done, the results can be passed
    to avoid unnecessary repeat calculations. If enough values are
    passed, they will be checked for consistency if chkvals is True.
    
    :arg float salt: Salinity in kg/kg.
    :arg float temp: Temperature in K.
    :arg float ppur: Pure water pressure in Pa.
    :arg dlpur: Pure water density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dlpur: float or None
    :arg psea: Seawater pressure in Pa. If unknown, pass None (default)
        and it will be calculated.
    :type psea: float or None
    :arg dlsea: Density of liquid water in seawater in kg/m3. If
        unknown, pass None (default) and it will be calculated.
    :type dlsea: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg dlpur0: Initial guess for the pure water density in kg/m3. If
        None (default) then `flu3a._dliq_default` is used.
    :type dlpur0: float or None
    :arg psea0: Initial guess for the seawater pressure in Pa. If None
        (default) then `_approx_stp` is used.
    :type psea0: float or None
    :arg dlsea0: Initial guess for the density of liquid water in
        seawater in kg/m3. If None (default) then `_approx_stp` is used.
    :type dlsea0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg bool useext: If False (default) then the salt contribution is
        calculated from _GSCOEFFS; if True, from _GSCOEFFS_EXT.
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Pure water density, seawater pressure, and seawater liquid
        density (all in SI units).
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    """
    if mathargs is None:
        mathargs = dict()
    if dlpur is None:
        dlpur = flu3a.eq_tp_liq(temp,ppur,dliq0=dlpur0,mathargs=mathargs)
    
    if any(val is None for val in (psea,dlsea)):
        x0 = (psea0,dlsea0)
        fargs = (salt,temp,ppur,dlpur)
        fkwargs = {'useext': useext}
        x1 = _newton(_diff_stp,x0,_approx_stp,fargs=fargs,fkwargs=fkwargs,
            **mathargs)
        psea, dlsea = x1
    
    _chkflubnds(temp,dlpur)
    _chkflubnds(temp,dlsea)
    _chksalbnds(salt,temp,psea)
    if not chkvals:
        return dlpur, psea, dlsea
    
    lhs1, rhs1, __, __ = flu3a._diff_tp(dlpur,temp,ppur)
    lhs2, rhs2, __, __ = _diff_stp(psea,dlsea,salt,temp,ppur,dlpur)
    lhs = numpy.insert(lhs2,0,lhs1)
    rhs = numpy.insert(rhs2,0,rhs1)
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
    return dlpur, psea, dlsea


## Thermodynamic functions
def osmoticpressure(salt,temp,ppur,dlpur=None,psea=None,dlsea=None,
    chkvals=False,chktol=_CHKTOL,dlpur0=None,psea0=None,dlsea0=None,
    chkbnd=False,useext=False,mathargs=None):
    """Calculate seawater osmotic pressure.
    
    Calculate the osmotic pressure of seawater, the additional pressure
    required to bring seawater of the given salinity and temperature to
    equilibrium with pure liquid water at the given temperature and
    pressure. This is the additional pressure required to prevent the
    flow of pure water across a semipermeable membrane.
    
    :arg float salt: Salinity in kg/kg.
    :arg float temp: Temperature in K.
    :arg float ppur: Pure water pressure in Pa.
    :arg dlpur: Pure water density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dlpur: float or None
    :arg psea: Seawater pressure in Pa. If unknown, pass None (default)
        and it will be calculated.
    :type psea: float or None
    :arg dlsea: Density of liquid water in seawater in kg/m3. If
        unknown, pass None (default) and it will be calculated.
    :type dlsea: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg dlpur0: Initial guess for the pure water density in kg/m3. If
        None (default) then `flu3a._dliq_default` is used.
    :type dlpur0: float or None
    :arg psea0: Initial guess for the seawater pressure in Pa. If None
        (default) then `_approx_stp` is used.
    :type psea0: float or None
    :arg dlsea0: Initial guess for the density of liquid water in
        seawater in kg/m3. If None (default) then `_approx_stp` is used.
    :type dlsea0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg bool useext: If False (default) then the salt contribution is
        calculated from _GSCOEFFS; if True, from _GSCOEFFS_EXT.
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Osmotic pressure in Pa.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> osmoticpressure(0.035,300.,1e5)
    2594603.20968
    """
    dlpur, psea, dlsea = eq_stp(salt,temp,ppur,dlpur=dlpur,psea=psea,
        dlsea=dlsea,chkvals=chkvals,chktol=chktol,dlpur0=dlpur0,psea0=psea0,
        dlsea0=dlsea0,chkbnd=chkbnd,useext=useext,mathargs=mathargs)
    posm = psea - ppur
    return posm

