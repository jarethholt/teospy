"""Humid air enthalpy and related properties.
 
This module provides the enthalpy of humid air (dry air with water
vapour) with dry air mass fraction, entropy, and pressure as the primary
variables. It also provides functions `pot*` for calculating air
properties under (dry) adiabatic displacement.

:Examples:

>>> air_h(0,0,0,0.9,900.,1e5)
274592.611783
>>> air_h(0,0,1,0.9,900.,1e5)
0.903262695636
>>> air_h(0,2,0,0.9,900.,1e5)
0.223684689765
>>> pottemp(0.9,300.,5e4,1e5)
363.653905688
>>> potdensity(0.9,300.,5e4,1e5)
0.903509489711
>>> potenthalpy(0.9,300.,5e4,1e5)
348872.568665

:Functions:

* :func:`eq_aep`: Calculate equilibrium quantities at dry fraction,
  entropy, and pressure.
* :func:`air_h`: Humid air enthalpy with derivatives.
* :func:`eq_pot`: Calculate equilibrium quantities at dry fraction,
  in-situ temperature, in-situ pressure, and potential pressure.
* :func:`pottemp`: Calculate the humid air potential temperature.
* :func:`potdensity`: Calculate the humid air potential density.
* :func:`potenthalpy`: Calculate the humid air potential enthalpy.

"""

__all__ = ['eq_aep','air_h','eq_pot','pottemp','potdensity','potenthalpy']

import warnings
import numpy
from teospy import constants0
from teospy import convert0
from teospy import air2
from teospy import maths3
from teospy import air3a
from teospy import air3b

_CHKTOL = constants0.CHKTOL
_RUNIV = constants0.RUNIV
_MWAT = constants0.MWAT
_MDRY = constants0.MDRY
_RWAT = constants0.RWAT
_RDRY = constants0.RDRY
_PATM = constants0.PATM
_TCELS = constants0.TCELS
_TTP = constants0.TTP
_PTPE = constants0.PTPE
_LLVTP = constants0.LLVTP
_CVAP = constants0.CVAP
_CDRY = constants0.CDRY

_chkhumbnds = constants0.chkhumbnds
_air_f = air2.air_f
_eq_entropy = air2.eq_entropy
_eq_pressure = air2.eq_pressure
_newton = maths3.newton


## Auxiliary functions for AEP calculations
def _approx_aep(airf,entr,pres):
    """Approximate humid air TD at AEP.
    
    Approximate the temperature and density of humid air from the dry
    air mass fraction, specific entropy, and pressure. This
    approximation is based on constant heat capacities and latent heats
    for dry air and water vapour.
    
    :arg float airf: Dry air mass fraction in kg/kg.
    :arg float entr: Entropy in J/kg/K.
    :arg float pres: Pressure in Pa.
    :returns: Humid air temperature in K and density in kg/m3.
    """
    airx = convert0.air_molfractiondry(airf)
    pdry = pres*airx
    pvap = pres*(1-airx)
    s0 = (airf*_CDRY*numpy.log(_TTP/_TCELS) + (1-airf)*_LLVTP/_TTP
        - airf*_RDRY*numpy.log(pdry/_PATM)
        - (1-airf)*_RWAT*numpy.log(pvap/_PTPE))
    ceff = airf*_CDRY + (1-airf)*_CVAP
    
    # Invert for temperature, calculate density from ideal gas law
    temp = _TTP * numpy.exp((entr-s0)/ceff)
    raw = _RUNIV * (airf/_MDRY + (1-airf)/_MWAT)
    dhum = pres / (raw*temp)
    return temp, dhum

def _diff_aep(t,d,airf,entr,pres):
    """Calculate humid air disequilibrium at AEP.
    
    Calculate both sides of the equations
    
        given entropy = entropy of humid air
        given pressure = pressure of humid air
    
    and their Jacobians with respect to temperature and density. Solving
    these equations gives the equilibrium humid air temperature and
    density for the given dry air mass fraction, specific entropy, and
    pressure.
    
    :arg float t: Temperature in K.
    :arg float d: Humid air density in kg/m3.
    :arg float airf: Dry air mass fraction in kg/kg.
    :arg float entr: Entropy in J/kg/K.
    :arg float pres: Pressure in Pa.
    :returns: Left-hand side of the equations, right-hand side,
        Jacobian of LHS, and Jacobian of RHS.
    :rtype: tuple(array(float))
    """
    shum = _eq_entropy(0,0,0,airf,t,d)
    phum = _eq_pressure(0,0,0,airf,t,d)
    lhs = numpy.array([entr, pres])
    rhs = numpy.array([shum, phum])
    
    shum_t = _eq_entropy(0,1,0,airf,t,d)
    shum_d = _eq_entropy(0,0,1,airf,t,d)
    phum_t = _eq_pressure(0,1,0,airf,t,d)
    phum_d = _eq_pressure(0,0,1,airf,t,d)
    dlhs = numpy.array([[0.,0.], [0.,0.]])
    drhs = numpy.array([[shum_t,shum_d], [phum_t,phum_d]])
    return lhs, rhs, dlhs, drhs

def eq_aep(airf,entr,pres,temp=None,dhum=None,chkvals=False,
    chktol=_CHKTOL,temp0=None,dhum0=None,chkbnd=False,mathargs=None):
    """Get primary variables at AEP.
    
    Get the values of the equilibrium humid air temperature and density
    for the given dry air mass fraction, specific entropy, and pressure.
    
    If the calculation has already been done, the results can be passed
    to avoid unnecessary repeat calculations. If enough values are
    passed, they will be checked for consistency if chkvals is True.
    
    :arg float airf: Dry air mass fraction in kg/kg.
    :arg float entr: Entropy in J/kg/K.
    :arg float pres: Pressure in Pa.
    :arg temp: Temperature in K. If unknown, pass None (default) and it
        will be calculated.
    :type temp: float or None
    :arg dhum: Humid air density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dhum: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg temp0: Initial guess for the temperature in K. If None
        (default) then `_approx_aep` is used.
    :type temp0: float or None
    :arg dhum0: Initial guess for the humid air density in kg/m3. If
        None (default) then `_approx_aep` is used.
    :type dhum0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Temperature and density (in SI units).
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    """
    if any(val is None for val in (temp,dhum)):
        x0 = (temp0,dhum0)
        fargs = (airf,entr,pres)
        if mathargs is None:
            mathargs = dict()
        x1 = _newton(_diff_aep,x0,_approx_aep,fargs=fargs,**mathargs)
        temp, dhum = x1
    
    _chkhumbnds(airf,temp,dhum,chkbnd=chkbnd)
    if not chkvals:
        return temp, dhum
    
    lhs, rhs, __, __ = _diff_aep(temp,dhum,airf,entr,pres)
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
    return temp, dhum


## Public AEP functions
def air_h(drva,drve,drvp,airf,entr,pres,temp=None,dhum=None,chkvals=False,
    chktol=_CHKTOL,temp0=None,dhum0=None,chkbnd=False,mathargs=None):
    """Calculate humid air enthalpy with derivatives.
    
    Calculate the specific enthalpy of humid air or its derivatives with
    respect to dry air mass fraction, entropy, and pressure.
    
    :arg float airf: Dry air mass fraction in kg/kg.
    :arg float entr: Entropy in J/kg/K.
    :arg float pres: Pressure in Pa.
    :arg temp: Temperature in K. If unknown, pass None (default) and it
        will be calculated.
    :type temp: float or None
    :arg dhum: Humid air density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dhum: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg temp0: Initial guess for the temperature in K. If None
        (default) then `_approx_aep` is used.
    :type temp0: float or None
    :arg dhum0: Initial guess for the humid air density in kg/m3. If
        None (default) then `_approx_aep` is used.
    :type dhum0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Enthalpy in units of
        (J/kg) / (kg/kg)^drva / (J/kg/K)^drve / Pa^drvp.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> air_h(0,0,0,0.9,900.,1e5)
    274592.611783
    >>> air_h(0,0,1,0.9,900.,1e5)
    0.903262695636
    >>> air_h(0,2,0,0.9,900.,1e5)
    0.223684689765
    """
    temp, dhum = eq_aep(airf,entr,pres,temp=temp,dhum=dhum,chkvals=chkvals,
        chktol=chktol,temp0=temp0,dhum0=dhum0,chkbnd=chkbnd,mathargs=mathargs)
    
    # Simple cases: First-order derivatives
    if (drva,drve,drvp) == (0,0,0):
        f = _air_f(0,0,0,airf,temp,dhum)
        f_t = _air_f(0,1,0,airf,temp,dhum)
        f_d = _air_f(0,0,1,airf,temp,dhum)
        h = f + dhum*f_d - temp*f_t
        return h
    elif (drva,drve,drvp) == (1,0,0):
        f_a = _air_f(1,0,0,airf,temp,dhum)
        h_a = f_a
        return h_a
    elif (drva,drve,drvp) == (0,1,0):
        h_s = temp
        return h_s
    elif (drva,drve,drvp) == (0,0,1):
        h_p = dhum**(-1)
        return h_p
    
    # Second-order derivatives require this denominator
    __, __, __, sp_td = _diff_aep(temp,dhum,airf,entr,pres)
    if (drva,drve,drvp) == (2,0,0):
        f_aa = _air_f(2,0,0,airf,temp,dhum)
        f_at = _air_f(1,1,0,airf,temp,dhum)
        f_ad = _air_f(1,0,1,airf,temp,dhum)
        s_a = _eq_entropy(1,0,0,airf,temp,dhum)
        p_a = _eq_pressure(1,0,0,airf,temp,dhum)
        sp_a = numpy.array([s_a,p_a])
        f_ax = numpy.array([f_at,f_ad])
        td_a = numpy.linalg.solve(sp_td,-sp_a)
        h_aa = f_aa + f_ax.dot(td_a)
        return h_aa
    elif (drva,drve,drvp) == (1,1,0):
        f_at = _air_f(1,1,0,airf,temp,dhum)
        f_ad = _air_f(1,0,1,airf,temp,dhum)
        sp_s = numpy.array([1.,0.])
        td_s = numpy.linalg.solve(sp_td,sp_s)
        f_ax = numpy.array([f_at,f_ad])
        h_as = f_ax.dot(td_s)
        return h_as
    elif (drva,drve,drvp) == (1,0,1):
        f_at = _air_f(1,1,0,airf,temp,dhum)
        f_ad = _air_f(1,0,1,airf,temp,dhum)
        sp_p = numpy.array([0.,1.])
        td_p = numpy.linalg.solve(sp_td,sp_p)
        f_ax = numpy.array([f_at,f_ad])
        h_ap = f_ax.dot(td_p)
        return h_ap
    elif (drva,drve,drvp) == (0,2,0):
        sp_s = numpy.array([1.,0.])
        td_s = numpy.linalg.solve(sp_td,sp_s)
        h_ss = td_s[0]
        return h_ss
    elif (drva,drve,drvp) == (0,1,1):
        sp_p = numpy.array([0.,1.])
        td_p = numpy.linalg.solve(sp_td,sp_p)
        h_sp = td_p[0]
        return h_sp
    elif (drva,drve,drvp) == (0,0,2):
        sp_p = numpy.array([0.,1.])
        td_p = numpy.linalg.solve(sp_td,sp_p)
        h_pp = -td_p[1] / dhum**2
        return h_pp
    
    # Should not have made it this far!
    errmsg = 'Derivatives {0} not recognized'.format((drva,drve,drvp))
    raise ValueError(errmsg)


## Adiabatic ascent/descent calculations
def eq_pot(airf,temp,pres,ppot,dhum=None,tpot=None,dpot=None,
    chkvals=False,chktol=_CHKTOL,dhum0=None,tpot0=None,dpot0=None,
    chkbnd=False,mathargs=None):
    """Get primary variables at ATP1P2.
    
    Get the values of the equilibrium in-situ density, potential
    temperature, and potential density for the given dry air mass
    fraction, in-situ temperature, in-situ pressure, and potential
    pressure.
    
    If the calculation has already been done, the results can be passed
    to avoid unnecessary repeat calculations. If enough values are
    passed, they will be checked for consistency if chkvals is True.
    
    :arg float airf: Dry air mass fraction in kg/kg.
    :arg float temp: In-situ temperature in K.
    :arg float pres: In-situ pressure in Pa.
    :arg float ppot: Potential pressure in Pa.
    :arg dhum: Humid air in-situ density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dhum: float or None
    :arg tpot: Potential temperature in K. If unknown, pass None
        (default) and it will be calculated.
    :type tpot: float or None
    :arg dpot: Potential density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dpot: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg dhum0: Initial guess for the in-situ density in kg/m3. If None
        (default) then `air3a._approx_atp` is used.
    :type dhum0: float or None
    :arg tpot0: Initial guess for the potential temperature in K. If
        None (default) then `_approx_aep` is used.
    :type tpot0: float or None
    :arg dpot0: Initial guess for the potential density in kg/m3. If
        None (default) then `_approx_aep` is used.
    :type dpot0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: In-situ density, potential temperature, and potential
        density (in SI units).
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    """
    if dhum is None:
        dhum = air3a.eq_atp(airf,temp,pres,dhum0=dhum0,chkbnd=chkbnd,
            mathargs=mathargs)
    _chkhumbnds(airf,temp,dhum,chkbnd=chkbnd)
    
    entr = air3b.entropy(airf,temp,pres,dhum=dhum,chkvals=False)
    if any(val is None for val in (tpot,dpot)):
        x0 = (tpot0,dpot0)
        fargs = (airf,entr,ppot)
        if mathargs is None:
            mathargs = dict()
        x1 = _newton(_diff_aep,x0,_approx_aep,fargs=fargs,**mathargs)
        tpot, dpot = [val2 if val1 is None else val1
            for (val1,val2) in zip((tpot,dpot),x1)]
    _chkhumbnds(airf,tpot,dpot,chkbnd=chkbnd)
    if not chkvals:
        return dhum, tpot, dpot
    
    lhs1, rhs1, __, __ = air3a._diff_atp(dhum,airf,temp,pres)
    lhs2, rhs2, __, __ = _diff_aep(tpot,dpot,airf,entr,ppot)
    lhs = [lhs1] + lhs2.tolist()
    rhs = [rhs1] + rhs2.tolist()
    errs = list()
    for (l,r) in zip(lhs,rhs):
        if abs(r) > chktol:
            errs.append(abs(l/r - 1))
        else:
            errs.append(abs(l-r))
    if max(errs) > chktol:
        warnmsg = ('Given values {0} and solutions {1} disagree to more than '
            'the tolerance {2}').format(lhs,rhs,chktol)
        warnings.warn(warnmsg,RuntimeWarning)
    return dhum, tpot, dpot

def pottemp(airf,temp,pres,ppot,dhum=None,tpot=None,dpot=None,
    chkvals=False,chktol=_CHKTOL,dhum0=None,tpot0=None,dpot0=None,
    chkbnd=False,mathargs=None):
    """Calculate humid air potential temperature.
    
    Calculate the potential temperature of humid air from the dry air
    mass fraction, in-situ temperature, in-situ pressure, and potential
    pressure.
    
    :arg float airf: Dry air mass fraction in kg/kg.
    :arg float temp: In-situ temperature in K.
    :arg float pres: In-situ pressure in Pa.
    :arg float ppot: Potential pressure in Pa.
    :arg dhum: Humid air in-situ density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dhum: float or None
    :arg tpot: Potential temperature in K. If unknown, pass None
        (default) and it will be calculated.
    :type tpot: float or None
    :arg dpot: Potential density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dpot: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg dhum0: Initial guess for the in-situ density in kg/m3. If None
        (default) then `air3a._approx_atp` is used.
    :type dhum0: float or None
    :arg tpot0: Initial guess for the potential temperature in K. If
        None (default) then `_approx_aep` is used.
    :type tpot0: float or None
    :arg dpot0: Initial guess for the potential density in kg/m3. If
        None (default) then `_approx_aep` is used.
    :type dpot0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Potential temperature in K.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> pottemp(0.9,300.,5e4,1e5)
    363.653905688
    """
    dhum, tpot, dpot = eq_pot(airf,temp,pres,ppot,dhum=dhum,tpot=tpot,
        dpot=dpot,chkvals=chkvals,chktol=chktol,dhum0=dhum0,tpot0=tpot0,
        dpot0=dpot0,chkbnd=chkbnd,mathargs=mathargs)
    return tpot

def potdensity(airf,temp,pres,ppot,dhum=None,tpot=None,dpot=None,
    chkvals=False,chktol=_CHKTOL,dhum0=None,tpot0=None,dpot0=None,chkbnd=False,
    mathargs=None):
    """Calculate humid air potential density.
    
    Calculate the potential density of humid air from the dry air mass
    fraction, in-situ temperature, in-situ pressure, and potential
    pressure.
    
    :arg float airf: Dry air mass fraction in kg/kg.
    :arg float temp: In-situ temperature in K.
    :arg float pres: In-situ pressure in Pa.
    :arg float ppot: Potential pressure in Pa.
    :arg dhum: Humid air in-situ density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dhum: float or None
    :arg tpot: Potential temperature in K. If unknown, pass None
        (default) and it will be calculated.
    :type tpot: float or None
    :arg dpot: Potential density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dpot: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg dhum0: Initial guess for the in-situ density in kg/m3. If None
        (default) then `air3a._approx_atp` is used.
    :type dhum0: float or None
    :arg tpot0: Initial guess for the potential temperature in K. If
        None (default) then `_approx_aep` is used.
    :type tpot0: float or None
    :arg dpot0: Initial guess for the potential density in kg/m3. If
        None (default) then `_approx_aep` is used.
    :type dpot0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Potential density in kg/m3.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> potdensity(0.9,300.,5e4,1e5)
    0.903509489711
    """
    dhum, tpot, dpot = eq_pot(airf,temp,pres,ppot,dhum=dhum,tpot=tpot,
        dpot=dpot,chkvals=chkvals,chktol=chktol,dhum0=dhum0,tpot0=tpot0,
        dpot0=dpot0,chkbnd=chkbnd,mathargs=mathargs)
    return dpot

def potenthalpy(airf,temp,pres,ppot,dhum=None,tpot=None,dpot=None,
    chkvals=False,chktol=_CHKTOL,dhum0=None,tpot0=None,dpot0=None,
    chkbnd=False,mathargs=None):
    """Calculate humid air potential enthalpy.
    
    Calculate the potential enthalpy of humid air from the dry air mass
    fraction, in-situ temperature, in-situ pressure, and potential
    pressure.
    
    :arg float airf: Dry air mass fraction in kg/kg.
    :arg float temp: In-situ temperature in K.
    :arg float pres: In-situ pressure in Pa.
    :arg float ppot: Potential pressure in Pa.
    :arg dhum: Humid air in-situ density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dhum: float or None
    :arg tpot: Potential temperature in K. If unknown, pass None
        (default) and it will be calculated.
    :type tpot: float or None
    :arg dpot: Potential density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dpot: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg dhum0: Initial guess for the in-situ density in kg/m3. If None
        (default) then `air3a._approx_atp` is used.
    :type dhum0: float or None
    :arg tpot0: Initial guess for the potential temperature in K. If
        None (default) then `_approx_aep` is used.
    :type tpot0: float or None
    :arg dpot0: Initial guess for the potential density in kg/m3. If
        None (default) then `_approx_aep` is used.
    :type dpot0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Potential enthalpy in J/kg.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> potenthalpy(0.9,300.,5e4,1e5)
    348872.568665
    """
    dhum, tpot, dpot = eq_pot(airf,temp,pres,ppot,dhum=dhum,tpot=tpot,
        dpot=dpot,chkvals=chkvals,chktol=chktol,dhum0=dhum0,tpot0=tpot0,
        dpot0=dpot0,chkbnd=chkbnd,mathargs=mathargs)
    hpot = air2.enthalpy(airf,tpot,dpot)
    return hpot

