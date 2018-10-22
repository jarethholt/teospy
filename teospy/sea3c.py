"""Seawater entropy and related properties.

This module provides the entropy of seawater (liquid water and salt) as a function of salinity, enthalpy, and pressure. It also provides potential temperature functions, and expansion and contraction coefficients with respect to in-situ temperature, potential temperature (entropy), and potential enthalpy.

The potential temperature, enthalpy, etc. functions are also provided in `sea3a` and `sea3b`. The primary difference in this module is that the functions accept any of the in-situ temperature, in-situ enthalpy, potential temperature, or potential enthalpy as the thermal argument.

:Examples:

>>> temperature(0.035,1e7,enth=1e5)
295.985682129
>>> temperature(0.035,1e7,ppot=1e5,hpot=1e5)
298.413424848
>>> pottemp(0.035,1e7,1e5,enth=1e5)
295.782199115
>>> pottemp(0.035,1e7,1e5,hpot=1e5)
298.194955182
>>> contraction_t(0.035,1e7,enth=1e5)
0.728755239644
>>> contraction_t(0.035,1e7,ppot=1e5,hpot=1e5)
0.726349317428
>>> contraction_t(0.035,1e7,temp=300.)
0.724913833446
>>> contraction_t(0.035,1e7,ppot=15,tpot=300.)
0.724714253918
>>> contraction_h(0.035,1e7,1e5,enth=1e5)
0.718452125957
>>> contraction_h(0.035,1e7,1e5,hpot=1e5)
0.714531922616
>>> contraction_h(0.035,1e7,1e5,temp=300.)
0.712069013013
>>> contraction_h(0.035,1e7,1e5,tpot=300.)
0.711718411190
>>> contraction_theta(0.035,1e7,1e5,enth=1e5)
0.728499505688
>>> contraction_theta(0.035,1e7,1e5,hpot=1e5)
0.726099732703
>>> contraction_theta(0.035,1e7,1e5,temp=300.)
0.724667977117
>>> contraction_theta(0.035,1e7,1e5,tpot=300.)
0.724468894946

:Functions:

* :func:`eq_shp`: Calculate equilibrium quantities at salinity,
  enthalpy, and pressure.
* :func:`eq_pot`: Calculate equilibrium quantities at salinity, in-situ
  pressure, potential pressure, and in-situ/potential
  enthalpy/temperature.
* :func:`temperature`: Seawater temperature.
* :func:`density`: Seawater density.
* :func:`entropy`: Seawater entropy.
* :func:`contraction_t`: Seawater haline contraction coefficient at
  constant in-situ temperature.
* :func:`expansion_t`: Seawater thermal expansion coefficient with
  respect to in-situ temperature.
* :func:`pottemp`: Seawater potential temperature.
* :func:`potdensity`: Seawater potential density.
* :func:`contraction_h`: Seawater haline contraction coefficient at
  constant potential enthalpy.
* :func:`contraction_theta`: Seawater haline contraction coefficient at
  constant potential temperature.
* :func:`expansion_h`: Seawater thermal expansion coefficient with
  respect to potential enthalpy.
* :func:`expansion_theta`: Seawater thermal expansion coefficient with
  respect to potential temperature.

"""

__all__ = ['eq_shp','eq_pot','temperature','density','entropy','contraction_t',
    'expansion_t','pottemp','potdensity','contraction_h','contraction_theta',
    'expansion_h','expansion_theta']

import warnings
import numpy
import constants0
import flu2
import sal2
import maths3
import flu3a
import sea3a
import sea3b

_TTP = constants0.TTP
_PTP = constants0.PTPE
_DLTP = constants0.DLTP
_CLIQ = constants0.CLIQ
_CHKTOL = constants0.CHKTOL
_chkflubnds = constants0.chkflubnds
_chksalbnds = constants0.chksalbnds
_flu_eq_pressure = flu2.eq_pressure
_flu_eq_entropy = flu2.eq_entropy
_flu_eq_enthalpy = flu2.eq_enthalpy
_sal_eq_entropy = sal2.eq_entropy
_sal_eq_enthalpy = sal2.eq_enthalpy
_newton = maths3.newton


## Equilibrium functions
def _approx_shp(salt,enth,pres):
    """Approximate TDl at SHP.
    
    Approximate the temperature and liquid water density of seawater
    for the given salinity, enthalpy, and pressure.
    
    :arg float salt: Salinity in kg/kg.
    :arg float enth: Enthalpy in J/kg.
    :arg float pres: Pressure in Pa.
    :returns: Temperature and liquid water density (in SI units).
    """
    H2TP = 6815.26049527339
    hsal = H2TP * salt
    hvol = (pres-_PTP)/_DLTP
    temp = _TTP + (enth - hsal - hvol) / _CLIQ
    dliq = flu3a._dliq_default(temp,pres)
    return temp, dliq

def _diff_shp(t,d,salt,enth,pres,useext=False):
    """Calculate seawater disequilibrium at SHP.
    
    Calculate both sides of the equations
    
        given enthalpy = enthalpy of seawater
        given pressure = pressure in liquid water
    
    and their Jacobians with respect to temperature and liquid water
    density. Solving these equations gives the temperature and liquid
    water density for the given salinity, enthalpy, and pressure.
    
    :arg float t: Temperature in K.
    :arg float d: Liquid water density in kg/m3.
    :arg float salt: Salinity in kg/kg.
    :arg float enth: Enthalpy in J/kg.
    :arg float pres: Pressure in Pa.
    :arg bool useext: If False (default) then the salt contribution is
        calculated from _GSCOEFFS in sal1; if True, from _GSCOEFFS_EXT.
    :returns: Left-hand side of the equation, right-hand side,
        Jacobian of LHS, and Jacobian of RHS.
    :rtype: tuple(array(float))
    """
    pl = _flu_eq_pressure(0,0,t,d)
    hl = _flu_eq_enthalpy(0,0,t,d)
    hs = _sal_eq_enthalpy(0,0,0,salt,t,pres,useext=useext)
    h = hl + hs
    lhs = numpy.array([pres, enth])
    rhs = numpy.array([pl, h])
    
    pl_t = _flu_eq_pressure(1,0,t,d)
    pl_d = _flu_eq_pressure(0,1,t,d)
    hl_t = _flu_eq_enthalpy(1,0,t,d)
    hl_d = _flu_eq_enthalpy(0,1,t,d)
    hs_t = _sal_eq_enthalpy(0,1,0,salt,t,pres,useext=useext)
    h_t = hl_t + hs_t
    h_d = hl_d
    dlhs = numpy.array([[0.,0.], [0.,0.]])
    drhs = numpy.array([[pl_t,pl_d], [h_t,h_d]])
    return lhs, rhs, dlhs, drhs

def eq_shp(salt,pres,enth=None,temp=None,dliq=None,chkvals=False,
    chktol=_CHKTOL,temp0=None,dliq0=None,chkbnd=False,useext=False,
    mathargs=None):
    """Get primary seawater variables at SP and H/T.
    
    Get the values of all equilibrium seawater variables for the given
    salinity, pressure, and either temperature or enthalpy.
    
    If the calculation has already been done, the results can be passed
    to avoid unnecessary repeat calculations. If enough values are
    passed, they will be checked for consistency if chkvals is True.
    
    :arg float salt: Salinity in kg/kg.
    :arg float pres: In-situ pressure in Pa.
    :arg enth: In-situ enthalpy in J/kg. If unknown, pass None (default)
        and it will be calculated.
    :type enth: float or None
    :arg temp: In-situ temperature in K. If unknown, pass None (default)
        and it will be calculated.
    :type temp: float or None
    :arg dliq: In-situ liquid water density in kg/m3. If unknown, pass
        None (default) and it will be calculated.
    :type dliq: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg temp0: Initial guess for the in-situ temperature in K. If None
        (default) then `_approx_shp` is used.
    :type temp0: float or None
    :arg dliq0: Initial guess for the in-situ liquid water density in
        kg/m3. If None (default) then `_approx_shp` is used.
    :type dliq0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Enthalpy, temperature, and liquid water density (in SI units).
    :raises ValueError: If both enth and temp are None.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    """
    if temp is None:
        if enth is None:
            errmsg = 'One of enth or temp must be provided'
            raise ValueError(errmsg)
        x0 = (temp0,dliq0)
        fargs = (salt,enth,pres)
        fkwargs = {'useext': useext}
        if mathargs is None:
            mathargs = dict()
        x1 = _newton(_diff_shp,x0,_approx_shp,fargs=fargs,fkwargs=fkwargs,
            **mathargs)
        temp, dliq = x1
    else:
        if dliq is None:
            dliq = flu3a.eq_tp_liq(temp,pres,dliq0=dliq0,mathargs=mathargs)
        if enth is None:
            enth = sea3a.enthalpy(salt,temp,pres,dliq=dliq,useext=useext)
    
    _chkflubnds(temp,dliq,chkbnd=chkbnd)
    _chksalbnds(salt,temp,pres,chkbnd=chkbnd)
    if not chkvals:
        return enth, temp, dliq
    
    lhs, rhs, __, __ = _diff_shp(temp,dliq,salt,enth,pres,useext=useext)
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
    return enth, temp, dliq

def eq_pot(salt,pres,ppot,enth=None,temp=None,dliq=None,hpot=None,
    tpot=None,dlpot=None,chkvals=False,chktol=_CHKTOL,temp0=None,
    dliq0=None,tpot0=None,dlpot0=None,chkbnd=False,useext=False,
    mathargs=None):
    """Get primary seawater variables at SP1P2 and H1/T1/H2/T2.
    
    Get the values of all equilibrium seawater variables for adiabatic
    ascent (potential) calculations of a variety of inputs. Salinity,
    in-situ pressure, and potential pressure must always be given. Then
    any one of in-situ enthalpy, in-situ temperature, potential
    enthalpy, or potential temperature must be given. Any missing
    quantities will be calculated and all of these quantities will
    returned.
    
    If the calculation has already been done, the results can be passed
    to avoid unnecessary repeat calculations. If enough values are
    passed, they will be checked for consistency if chkvals is True.
    
    :arg float salt: Salinity in kg/kg.
    :arg float pres: In-situ pressure in Pa.
    :arg float ppot: Potential pressure in Pa.
    :arg enth: In-situ enthalpy in J/kg. If unknown, pass None (default)
        and it will be calculated.
    :type enth: float or None
    :arg temp: In-situ temperature in K. If unknown, pass None (default)
        and it will be calculated.
    :type temp: float or None
    :arg dliq: In-situ liquid water density in kg/m3. If unknown, pass
        None (default) and it will be calculated.
    :type dliq: float or None
    :arg hpot: Potential enthalpy in J/kg. If unknown, pass None
        (default) and it will be calculated.
    :type hpot: float or None
    :arg tpot: Potential temperature in K. If unknown, pass None
        (default) and it will be calculated.
    :type tpot: float or None
    :arg dlpot: Potential liquid water density in kg/m3. If unknown,
        pass None (default) and it will be calculated.
    :type dlpot: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg temp0: Initial guess for the in-situ temperature in K. If None
        (default) then `_approx_shp` is used.
    :type temp0: float or None
    :arg dliq0: Initial guess for the in-situ liquid water density in
        kg/m3. If None (default) then `_approx_shp` is used.
    :type dliq0: float or None
    :arg tpot0: Initial guess for the potential temperature in K. If
        None (default) then `sea3b._approx_sep` is used.
    :type tpot0: float or None
    :arg dlpot0: Initial guess for the potential liquid water density in
        kg/m3. If None (default) then `sea3b._approx_sep` is used.
    :type dlpot0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: In-situ enthalpy, in-situ temperature, in-situ liquid
        water density, potential enthalpy, potential temperature, and
        potential liquid water density (all in SI units).
    :raises ValueError: If insufficient or incompatible values are
        given.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    """
    if all(val is None for val in (enth,temp,hpot,tpot)):
        errmsg = 'One of enth, temp, hpot, or tpot must be provided'
        raise ValueError(errmsg)
    if temp is None and enth is None:
        if tpot is None:
            hpot, tpot, dlpot = eq_shp(salt,ppot,enth=hpot,dliq=dlpot,
                temp0=tpot0,dliq0=dlpot0,useext=useext,mathargs=mathargs)
        else:
            if dlpot is None:
                dlpot = flu3a.eq_tp_liq(tpot,ppot,dliq0=dlpot0,
                    mathargs=mathargs)
            if hpot is None:
                hpot = sea3a.enthalpy(salt,tpot,ppot,dliq=dlpot,useext=useext)
        
        __, __, __, temp, dliq = sea3b.eq_pot(salt,ppot,pres,temp=tpot,
            dliq=dlpot,tpot0=temp0,dlpot0=dliq0,useext=useext,mathargs=mathargs)
        enth = sea3a.enthalpy(salt,temp,pres,dliq=dliq,useext=useext)
    else:
        if temp is None:
            temp, dliq = eq_shp(salt,enth,pres,temp0=temp0,dliq0=dliq0,
                useext=useext,mathargs=mathargs)
        else:
            if dliq is None:
                dliq = flu3a.eq_tp_liq(temp,pres,dliq0=dliq0,mathargs=mathargs)
            if enth is None:
                enth = sea3a.enthalpy(salt,temp,pres,dliq=dliq,useext=useext)
        
        if tpot is None:
            __, __, __, tpot, dlpot = sea3b.eq_pot(salt,pres,ppot,temp=temp,
                dliq=dliq,tpot0=tpot0,dlpot0=dlpot0,useext=useext,
                mathargs=mathargs)
        if dlpot is None:
            dlpot = flu3a.eq_tp_liq(tpot,ppot,dliq0=dlpot0,mathargs=mathargs)
        if hpot is None:
            hpot = sea3a.enthalpy(salt,tpot,ppot,dliq=dlpot,useext=useext)
    
    _chkflubnds(temp,dliq,chkbnd=chkbnd)
    _chksalbnds(salt,temp,pres,chkbnd=chkbnd)
    _chkflubnds(tpot,dlpot,chkbnd=chkbnd)
    _chksalbnds(salt,tpot,ppot,chkbnd=chkbnd)
    if not chkvals:
        return enth, temp, dliq, hpot, tpot, dlpot
    
    # Check for consistency
    lhs1, rhs1, __, __ = _diff_shp(temp,dliq,salt,enth,pres,useext=useext)
    lhs2, rhs2, __, __ = _diff_shp(tpot,dlpot,salt,hpot,ppot,useext=useext)
    lhs = numpy.concatenate((lhs1,lhs2))
    rhs = numpy.concatenate((rhs1,rhs2))
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
    return enth, temp, dliq, hpot, tpot, dlpot


## In-situ thermodynamic properties
def temperature(salt,pres,ppot=None,enth=None,temp=None,dliq=None,
    hpot=None,tpot=None,dlpot=None,chkvals=False,chktol=_CHKTOL,
    temp0=None,dliq0=None,tpot0=None,dlpot0=None,chkbnd=False,
    useext=False,mathargs=None):
    """Calculate seawater temperature.
    
    Calculate the in-situ temperature of seawater under the given
    conditions, using any of in-situ/potential enthalpy/temperature as
    the thermal variable.
    
    :arg float salt: Salinity in kg/kg.
    :arg float pres: In-situ pressure in Pa.
    :arg ppot: Potential pressure in Pa. Use None (default) for purely
        in-situ calculations.
    :type ppot: float or None
    :arg enth: In-situ enthalpy in J/kg. If unknown, pass None (default)
        and it will be calculated.
    :type enth: float or None
    :arg temp: In-situ temperature in K. If unknown, pass None (default)
        and it will be calculated.
    :type temp: float or None
    :arg dliq: In-situ liquid water density in kg/m3. If unknown, pass
        None (default) and it will be calculated.
    :type dliq: float or None
    :arg hpot: Potential enthalpy in J/kg. If unknown, pass None
        (default) and it will be calculated.
    :type hpot: float or None
    :arg tpot: Potential temperature in K. If unknown, pass None
        (default) and it will be calculated.
    :type tpot: float or None
    :arg dlpot: Potential liquid water density in kg/m3. If unknown,
        pass None (default) and it will be calculated.
    :type dlpot: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg temp0: Initial guess for the in-situ temperature in K. If None
        (default) then `_approx_shp` is used.
    :type temp0: float or None
    :arg dliq0: Initial guess for the in-situ liquid water density in
        kg/m3. If None (default) then `_approx_shp` is used.
    :type dliq0: float or None
    :arg tpot0: Initial guess for the potential temperature in K. If
        None (default) then `sea3b._approx_sep` is used.
    :type tpot0: float or None
    :arg dlpot0: Initial guess for the potential liquid water density in
        kg/m3. If None (default) then `sea3b._approx_sep` is used.
    :type dlpot0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg bool useext: If False (default) then the salt contribution is
        calculated from _GSCOEFFS; if True, from _GSCOEFFS_EXT.
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: In-situ temperature in K.
    :raises ValueError: If insufficient or incompatible values are
        given.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> temperature(0.035,1e7,enth=1e5)
    295.985682129
    >>> temperature(0.035,1e7,ppot=1e5,hpot=1e5)
    298.413424848
    """
    if any(val is not None for val in (enth,temp)):
        enth, temp, dliq = eq_shp(salt,pres,enth=enth,temp=temp,dliq=dliq,
            chkvals=chkvals,chktol=chktol,temp0=temp0,dliq0=dliq0,useext=useext,
            mathargs=mathargs)
    else:
        if ppot is None:
            errmsg = 'ppot must be provided if both enth and temp are None'
            raise ValueError(errmsg)
        res = eq_pot(salt,pres,ppot,hpot=hpot,tpot=tpot,dlpot=dlpot,
            chkvals=chkvals,chktol=chktol,temp0=temp0,dliq0=dliq0,tpot0=tpot0,
            dlpot0=dlpot0,useext=useext,mathargs=mathargs)
        enth, temp, dliq = res[:3]
    _chkflubnds(temp,dliq,chkbnd=chkbnd)
    _chksalbnds(salt,temp,pres,chkbnd=chkbnd)
    return temp

def density(salt,pres,ppot=None,enth=None,temp=None,dliq=None,hpot=None,
    tpot=None,dlpot=None,chkvals=False,chktol=_CHKTOL,temp0=None,
    dliq0=None,tpot0=None,dlpot0=None,chkbnd=False,useext=False,
    mathargs=None):
    """Calculate seawater density.
    
    Calculate the in-situ density of seawater under the given
    conditions, using any of in-situ/potential enthalpy/temperature as
    the thermal variable.
    
    :arg float salt: Salinity in kg/kg.
    :arg float pres: In-situ pressure in Pa.
    :arg ppot: Potential pressure in Pa. Use None (default) for purely
        in-situ calculations.
    :type ppot: float or None
    :arg enth: In-situ enthalpy in J/kg. If unknown, pass None (default)
        and it will be calculated.
    :type enth: float or None
    :arg temp: In-situ temperature in K. If unknown, pass None (default)
        and it will be calculated.
    :type temp: float or None
    :arg dliq: In-situ liquid water density in kg/m3. If unknown, pass
        None (default) and it will be calculated.
    :type dliq: float or None
    :arg hpot: Potential enthalpy in J/kg. If unknown, pass None
        (default) and it will be calculated.
    :type hpot: float or None
    :arg tpot: Potential temperature in K. If unknown, pass None
        (default) and it will be calculated.
    :type tpot: float or None
    :arg dlpot: Potential liquid water density in kg/m3. If unknown,
        pass None (default) and it will be calculated.
    :type dlpot: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg temp0: Initial guess for the in-situ temperature in K. If None
        (default) then `_approx_shp` is used.
    :type temp0: float or None
    :arg dliq0: Initial guess for the in-situ liquid water density in
        kg/m3. If None (default) then `_approx_shp` is used.
    :type dliq0: float or None
    :arg tpot0: Initial guess for the potential temperature in K. If
        None (default) then `sea3b._approx_sep` is used.
    :type tpot0: float or None
    :arg dlpot0: Initial guess for the potential liquid water density in
        kg/m3. If None (default) then `sea3b._approx_sep` is used.
    :type dlpot0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg bool useext: If False (default) then the salt contribution is
        calculated from _GSCOEFFS; if True, from _GSCOEFFS_EXT.
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: In-situ density in kg/m3.
    :raises ValueError: If insufficient or incompatible values are
        given.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> density(0.035,1e7,enth=1e5)
    1028.10986556
    >>> density(0.035,1e7,ppot=1e5,hpot=1e5)
    1027.36529798
    """
    if any(val is not None for val in (enth,temp)):
        enth, temp, dliq = eq_shp(salt,pres,enth=enth,temp=temp,dliq=dliq,
            chkvals=chkvals,chktol=chktol,temp0=temp0,dliq0=dliq0,useext=useext,
            mathargs=mathargs)
    else:
        if ppot is None:
            errmsg = 'ppot must be provided if both enth and temp are None'
            raise ValueError(errmsg)
        res = eq_pot(salt,pres,ppot,hpot=hpot,tpot=tpot,dlpot=dlpot,
            chkvals=chkvals,chktol=chktol,temp0=temp0,dliq0=dliq0,tpot0=tpot0,
            dlpot0=dlpot0,useext=useext,mathargs=mathargs)
        enth, temp, dliq = res[:3]
    _chkflubnds(temp,dliq,chkbnd=chkbnd)
    _chksalbnds(salt,temp,pres,chkbnd=chkbnd)
    dsea = sea3a.density(salt,temp,pres,dliq=dliq,useext=useext)
    return dsea

def entropy(salt,pres,ppot=None,enth=None,temp=None,dliq=None,hpot=None,
    tpot=None,dlpot=None,chkvals=False,chktol=_CHKTOL,temp0=None,
    dliq0=None,tpot0=None,dlpot0=None,chkbnd=False,useext=False,
    mathargs=None):
    """Calculate seawater entropy.
    
    Calculate the entropy of seawater under the given conditions, using
    any of in-situ/potential enthalpy/temperature as the thermal
    variable.
    
    :arg float salt: Salinity in kg/kg.
    :arg float pres: In-situ pressure in Pa.
    :arg ppot: Potential pressure in Pa. Use None (default) for purely
        in-situ calculations.
    :type ppot: float or None
    :arg enth: In-situ enthalpy in J/kg. If unknown, pass None (default)
        and it will be calculated.
    :type enth: float or None
    :arg temp: In-situ temperature in K. If unknown, pass None (default)
        and it will be calculated.
    :type temp: float or None
    :arg dliq: In-situ liquid water density in kg/m3. If unknown, pass
        None (default) and it will be calculated.
    :type dliq: float or None
    :arg hpot: Potential enthalpy in J/kg. If unknown, pass None
        (default) and it will be calculated.
    :type hpot: float or None
    :arg tpot: Potential temperature in K. If unknown, pass None
        (default) and it will be calculated.
    :type tpot: float or None
    :arg dlpot: Potential liquid water density in kg/m3. If unknown,
        pass None (default) and it will be calculated.
    :type dlpot: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg temp0: Initial guess for the in-situ temperature in K. If None
        (default) then `_approx_shp` is used.
    :type temp0: float or None
    :arg dliq0: Initial guess for the in-situ liquid water density in
        kg/m3. If None (default) then `_approx_shp` is used.
    :type dliq0: float or None
    :arg tpot0: Initial guess for the potential temperature in K. If
        None (default) then `sea3b._approx_sep` is used.
    :type tpot0: float or None
    :arg dlpot0: Initial guess for the potential liquid water density in
        kg/m3. If None (default) then `sea3b._approx_sep` is used.
    :type dlpot0: float or None
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
    :raises ValueError: If insufficient or incompatible values are
        given.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> entropy(0.035,1e5,enth=1e5)
    350.310622663
    """
    if any(val is not None for val in (enth,temp)):
        enth, temp, dliq = eq_shp(salt,pres,enth=enth,temp=temp,dliq=dliq,
            chkvals=chkvals,chktol=chktol,temp0=temp0,dliq0=dliq0,useext=useext,
            mathargs=mathargs)
    else:
        if ppot is None:
            errmsg = 'ppot must be provided if both enth and temp are None'
            raise ValueError(errmsg)
        res = eq_pot(salt,pres,ppot,hpot=hpot,tpot=tpot,dlpot=dlpot,
            chkvals=chkvals,chktol=chktol,temp0=temp0,dliq0=dliq0,tpot0=tpot0,
            dlpot0=dlpot0,useext=useext,mathargs=mathargs)
        enth, temp, dliq = res[:3]
    _chkflubnds(temp,dliq,chkbnd=chkbnd)
    _chksalbnds(salt,temp,pres,chkbnd=chkbnd)
    entr = sea3a.entropy(salt,temp,pres,dliq=dliq,useext=useext)
    return entr

def contraction_t(salt,pres,ppot=None,enth=None,temp=None,dliq=None,
    hpot=None,tpot=None,dlpot=None,chkvals=False,chktol=_CHKTOL,
    temp0=None,dliq0=None,tpot0=None,dlpot0=None,chkbnd=False,
    useext=False,mathargs=None):
    """Calculate isothermal haline contraction coefficient.
    
    Calculate the haline contraction coefficient of seawater at constant
    in-situ temperature under the given conditions, using any of
    in-situ/potential enthalpy/temperature as the thermal variable.
    
    :arg float salt: Salinity in kg/kg.
    :arg float pres: In-situ pressure in Pa.
    :arg ppot: Potential pressure in Pa. Use None (default) for purely
        in-situ calculations.
    :type ppot: float or None
    :arg enth: In-situ enthalpy in J/kg. If unknown, pass None (default)
        and it will be calculated.
    :type enth: float or None
    :arg temp: In-situ temperature in K. If unknown, pass None (default)
        and it will be calculated.
    :type temp: float or None
    :arg dliq: In-situ liquid water density in kg/m3. If unknown, pass
        None (default) and it will be calculated.
    :type dliq: float or None
    :arg hpot: Potential enthalpy in J/kg. If unknown, pass None
        (default) and it will be calculated.
    :type hpot: float or None
    :arg tpot: Potential temperature in K. If unknown, pass None
        (default) and it will be calculated.
    :type tpot: float or None
    :arg dlpot: Potential liquid water density in kg/m3. If unknown,
        pass None (default) and it will be calculated.
    :type dlpot: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg temp0: Initial guess for the in-situ temperature in K. If None
        (default) then `_approx_shp` is used.
    :type temp0: float or None
    :arg dliq0: Initial guess for the in-situ liquid water density in
        kg/m3. If None (default) then `_approx_shp` is used.
    :type dliq0: float or None
    :arg tpot0: Initial guess for the potential temperature in K. If
        None (default) then `sea3b._approx_sep` is used.
    :type tpot0: float or None
    :arg dlpot0: Initial guess for the potential liquid water density in
        kg/m3. If None (default) then `sea3b._approx_sep` is used.
    :type dlpot0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg bool useext: If False (default) then the salt contribution is
        calculated from _GSCOEFFS; if True, from _GSCOEFFS_EXT.
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Contraction coefficient in 1/(kg/kg).
    :raises ValueError: If insufficient or incompatible values are
        given.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> contraction_t(0.035,1e7,enth=1e5)
    0.728755239644
    >>> contraction_t(0.035,1e7,ppot=1e5,hpot=1e5)
    0.726349317428
    >>> contraction_t(0.035,1e7,temp=300.)
    0.724913833446
    >>> contraction_t(0.035,1e7,ppot=15,tpot=300.)
    0.724714253918
    """
    if any(val is not None for val in (enth,temp)):
        enth, temp, dliq = eq_shp(salt,pres,enth=enth,temp=temp,dliq=dliq,
            chkvals=chkvals,chktol=chktol,temp0=temp0,dliq0=dliq0,useext=useext,
            mathargs=mathargs)
    else:
        if ppot is None:
            errmsg = 'ppot must be provided if both enth and temp are None'
            raise ValueError(errmsg)
        res = eq_pot(salt,pres,ppot,hpot=hpot,tpot=tpot,dlpot=dlpot,
            chkvals=chkvals,chktol=chktol,temp0=temp0,dliq0=dliq0,tpot0=tpot0,
            dlpot0=dlpot0,useext=useext,mathargs=mathargs)
        enth, temp, dliq = res[:3]
    _chkflubnds(temp,dliq,chkbnd=chkbnd)
    _chksalbnds(salt,temp,pres,chkbnd=chkbnd)
    beta = sea3a.contraction_t(salt,temp,pres,dliq=dliq,useext=useext)
    return beta

def expansion_t(salt,pres,ppot=None,enth=None,temp=None,dliq=None,
    hpot=None,tpot=None,dlpot=None,chkvals=False,chktol=_CHKTOL,
    temp0=None,dliq0=None,tpot0=None,dlpot0=None,chkbnd=False,
    useext=False,mathargs=None):
    """Calculate expansion coefficient wrt in-situ temperature.
    
    Calculate the expansion coefficient of seawater with respect to
    in-situ temperature under the given conditions, using any of
    in-situ/potential enthalpy/temperature as the thermal variable.
    
    :arg float salt: Salinity in kg/kg.
    :arg float pres: In-situ pressure in Pa.
    :arg ppot: Potential pressure in Pa. Use None (default) for purely
        in-situ calculations.
    :type ppot: float or None
    :arg enth: In-situ enthalpy in J/kg. If unknown, pass None (default)
        and it will be calculated.
    :type enth: float or None
    :arg temp: In-situ temperature in K. If unknown, pass None (default)
        and it will be calculated.
    :type temp: float or None
    :arg dliq: In-situ liquid water density in kg/m3. If unknown, pass
        None (default) and it will be calculated.
    :type dliq: float or None
    :arg hpot: Potential enthalpy in J/kg. If unknown, pass None
        (default) and it will be calculated.
    :type hpot: float or None
    :arg tpot: Potential temperature in K. If unknown, pass None
        (default) and it will be calculated.
    :type tpot: float or None
    :arg dlpot: Potential liquid water density in kg/m3. If unknown,
        pass None (default) and it will be calculated.
    :type dlpot: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg temp0: Initial guess for the in-situ temperature in K. If None
        (default) then `_approx_shp` is used.
    :type temp0: float or None
    :arg dliq0: Initial guess for the in-situ liquid water density in
        kg/m3. If None (default) then `_approx_shp` is used.
    :type dliq0: float or None
    :arg tpot0: Initial guess for the potential temperature in K. If
        None (default) then `sea3b._approx_sep` is used.
    :type tpot0: float or None
    :arg dlpot0: Initial guess for the potential liquid water density in
        kg/m3. If None (default) then `sea3b._approx_sep` is used.
    :type dlpot0: float or None
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
    :raises ValueError: If insufficient or incompatible values are
        given.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> expansion_t(0.035,1e7,enth=1e5)
    2.89480851145e-4
    >>> expansion_t(0.035,1e7,ppot=1e5,hpot=1e5)
    3.07242256461e-4
    >>> expansion_t(0.035,1e7,temp=300.)
    3.18513471410e-4
    >>> expansion_t(0.035,1e7,ppot=15,tpot=300.)
    3.20122324740e-4
    """
    if any(val is not None for val in (enth,temp)):
        enth, temp, dliq = eq_shp(salt,pres,enth=enth,temp=temp,dliq=dliq,
            chkvals=chkvals,chktol=chktol,temp0=temp0,dliq0=dliq0,useext=useext,
            mathargs=mathargs)
    else:
        if ppot is None:
            errmsg = 'ppot must be provided if both enth and temp are None'
            raise ValueError(errmsg)
        res = eq_pot(salt,pres,ppot,hpot=hpot,tpot=tpot,dlpot=dlpot,
            chkvals=chkvals,chktol=chktol,temp0=temp0,dliq0=dliq0,tpot0=tpot0,
            dlpot0=dlpot0,useext=useext,mathargs=mathargs)
        enth, temp, dliq = res[:3]
    _chkflubnds(temp,dliq,chkbnd=chkbnd)
    _chksalbnds(salt,temp,pres,chkbnd=chkbnd)
    alpha = sea3a.expansion_t(salt,temp,pres,dliq=dliq,useext=useext)
    return alpha


## Potential thermodynamic properties
def pottemp(salt,pres,ppot,enth=None,temp=None,dliq=None,hpot=None,
    tpot=None,dlpot=None,chkvals=False,chktol=_CHKTOL,temp0=None,
    dliq0=None,tpot0=None,dlpot0=None,chkbnd=False,useext=False,
    mathargs=None):
    """Calculate seawater potential temperature.
    
    Calculate the potential temperature of seawater under the given
    conditions, using any of in-situ/potential enthalpy/temperature as
    the thermal variable.
    
    :arg float salt: Salinity in kg/kg.
    :arg float pres: In-situ pressure in Pa.
    :arg float ppot: Potential pressure in Pa.
    :arg enth: In-situ enthalpy in J/kg. If unknown, pass None (default)
        and it will be calculated.
    :type enth: float or None
    :arg temp: In-situ temperature in K. If unknown, pass None (default)
        and it will be calculated.
    :type temp: float or None
    :arg dliq: In-situ liquid water density in kg/m3. If unknown, pass
        None (default) and it will be calculated.
    :type dliq: float or None
    :arg hpot: Potential enthalpy in J/kg. If unknown, pass None
        (default) and it will be calculated.
    :type hpot: float or None
    :arg tpot: Potential temperature in K. If unknown, pass None
        (default) and it will be calculated.
    :type tpot: float or None
    :arg dlpot: Potential liquid water density in kg/m3. If unknown,
        pass None (default) and it will be calculated.
    :type dlpot: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg temp0: Initial guess for the in-situ temperature in K. If None
        (default) then `_approx_shp` is used.
    :type temp0: float or None
    :arg dliq0: Initial guess for the in-situ liquid water density in
        kg/m3. If None (default) then `_approx_shp` is used.
    :type dliq0: float or None
    :arg tpot0: Initial guess for the potential temperature in K. If
        None (default) then `sea3b._approx_sep` is used.
    :type tpot0: float or None
    :arg dlpot0: Initial guess for the potential liquid water density in
        kg/m3. If None (default) then `sea3b._approx_sep` is used.
    :type dlpot0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg bool useext: If False (default) then the salt contribution is
        calculated from _GSCOEFFS; if True, from _GSCOEFFS_EXT.
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Potential temperature in K.
    :raises ValueError: If insufficient or incompatible values are
        given.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> pottemp(0.035,1e7,1e5,enth=1e5)
    295.782199115
    >>> pottemp(0.035,1e7,1e5,hpot=1e5)
    298.194955182
    """
    enth, temp, dliq, hpot, tpot, dlpot = eq_pot(salt,pres,ppot,enth=enth,
        temp=temp,dliq=dliq,hpot=hpot,tpot=tpot,dlpot=dlpot,chkvals=chkvals,
        chktol=chktol,temp0=temp0,dliq0=dliq0,tpot0=tpot0,dlpot0=dlpot0,
        chkbnd=chkbnd,useext=useext,mathargs=mathargs)
    return tpot

def potdensity(salt,pres,ppot,enth=None,temp=None,dliq=None,hpot=None,
    tpot=None,dlpot=None,chkvals=False,chktol=_CHKTOL,temp0=None,
    dliq0=None,tpot0=None,dlpot0=None,chkbnd=False,useext=False,
    mathargs=None):
    """Calculate seawater potential density.
    
    Calculate the potential density under the given conditions, using
    any of in-situ/potential enthalpy/temperature as the thermal
    variable.
    
    :arg float salt: Salinity in kg/kg.
    :arg float pres: In-situ pressure in Pa.
    :arg float ppot: Potential pressure in Pa.
    :arg enth: In-situ enthalpy in J/kg. If unknown, pass None (default)
        and it will be calculated.
    :type enth: float or None
    :arg temp: In-situ temperature in K. If unknown, pass None (default)
        and it will be calculated.
    :type temp: float or None
    :arg dliq: In-situ liquid water density in kg/m3. If unknown, pass
        None (default) and it will be calculated.
    :type dliq: float or None
    :arg hpot: Potential enthalpy in J/kg. If unknown, pass None
        (default) and it will be calculated.
    :type hpot: float or None
    :arg tpot: Potential temperature in K. If unknown, pass None
        (default) and it will be calculated.
    :type tpot: float or None
    :arg dlpot: Potential liquid water density in kg/m3. If unknown,
        pass None (default) and it will be calculated.
    :type dlpot: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg temp0: Initial guess for the in-situ temperature in K. If None
        (default) then `_approx_shp` is used.
    :type temp0: float or None
    :arg dliq0: Initial guess for the in-situ liquid water density in
        kg/m3. If None (default) then `_approx_shp` is used.
    :type dliq0: float or None
    :arg tpot0: Initial guess for the potential temperature in K. If
        None (default) then `sea3b._approx_sep` is used.
    :type tpot0: float or None
    :arg dlpot0: Initial guess for the potential liquid water density in
        kg/m3. If None (default) then `sea3b._approx_sep` is used.
    :type dlpot0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg bool useext: If False (default) then the salt contribution is
        calculated from _GSCOEFFS; if True, from _GSCOEFFS_EXT.
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Potential density in kg/m3.
    :raises ValueError: If insufficient or incompatible values are
        given.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> potdensity(0.035,1e7,1e5,enth=1e5)
    1023.91737474
    >>> potdensity(0.035,1e7,1e5,hpot=1e5)
    1023.20527737
    """
    enth, temp, dliq, hpot, tpot, dlpot = eq_pot(salt,pres,ppot,enth=enth,
        temp=temp,dliq=dliq,hpot=hpot,tpot=tpot,dlpot=dlpot,chkvals=chkvals,
        chktol=chktol,temp0=temp0,dliq0=dliq0,tpot0=tpot0,dlpot0=dlpot0,
        chkbnd=chkbnd,useext=useext,mathargs=mathargs)
    dpot = sea3a.density(salt,tpot,ppot,dliq=dlpot,useext=useext)
    return dpot

def contraction_h(salt,pres,ppot,enth=None,temp=None,dliq=None,
    hpot=None,tpot=None,dlpot=None,chkvals=False,chktol=_CHKTOL,
    temp0=None,dliq0=None,tpot0=None,dlpot0=None,chkbnd=False,
    useext=False,mathargs=None):
    """Calculate seawater isenthalpic haline contraction coefficient.
    
    Calculate the haline contraction coefficient of seawater at constant
    potential enthalpy under the given conditions, using any of
    in-situ/potential enthalpy/temperature as the thermal variable.
    
    :arg float salt: Salinity in kg/kg.
    :arg float pres: In-situ pressure in Pa.
    :arg float ppot: Potential pressure in Pa.
    :arg enth: In-situ enthalpy in J/kg. If unknown, pass None (default)
        and it will be calculated.
    :type enth: float or None
    :arg temp: In-situ temperature in K. If unknown, pass None (default)
        and it will be calculated.
    :type temp: float or None
    :arg dliq: In-situ liquid water density in kg/m3. If unknown, pass
        None (default) and it will be calculated.
    :type dliq: float or None
    :arg hpot: Potential enthalpy in J/kg. If unknown, pass None
        (default) and it will be calculated.
    :type hpot: float or None
    :arg tpot: Potential temperature in K. If unknown, pass None
        (default) and it will be calculated.
    :type tpot: float or None
    :arg dlpot: Potential liquid water density in kg/m3. If unknown,
        pass None (default) and it will be calculated.
    :type dlpot: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg temp0: Initial guess for the in-situ temperature in K. If None
        (default) then `_approx_shp` is used.
    :type temp0: float or None
    :arg dliq0: Initial guess for the in-situ liquid water density in
        kg/m3. If None (default) then `_approx_shp` is used.
    :type dliq0: float or None
    :arg tpot0: Initial guess for the potential temperature in K. If
        None (default) then `sea3b._approx_sep` is used.
    :type tpot0: float or None
    :arg dlpot0: Initial guess for the potential liquid water density in
        kg/m3. If None (default) then `sea3b._approx_sep` is used.
    :type dlpot0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg bool useext: If False (default) then the salt contribution is
        calculated from _GSCOEFFS; if True, from _GSCOEFFS_EXT.
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Contraction coefficient in 1/(kg/kg).
    :raises ValueError: If insufficient or incompatible values are
        given.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> contraction_h(0.035,1e7,1e5,enth=1e5)
    0.718452125957
    >>> contraction_h(0.035,1e7,1e5,hpot=1e5)
    0.714531922616
    >>> contraction_h(0.035,1e7,1e5,temp=300.)
    0.712069013013
    >>> contraction_h(0.035,1e7,1e5,tpot=300.)
    0.711718411190
    """
    enth, temp, dliq, hpot, tpot, dlpot = eq_pot(salt,pres,ppot,enth=enth,
        temp=temp,dliq=dliq,hpot=hpot,tpot=tpot,dlpot=dlpot,chkvals=chkvals,
        chktol=chktol,temp0=temp0,dliq0=dliq0,tpot0=tpot0,dlpot0=dlpot0,
        chkbnd=chkbnd,useext=useext,mathargs=mathargs)
    beta = sea3b.contraction_h(salt,pres,ppot,temp=temp,dliq=dliq,tpot=tpot,
        dlpot=dlpot,useext=useext)
    return beta

def contraction_theta(salt,pres,ppot,enth=None,temp=None,dliq=None,
    hpot=None,tpot=None,dlpot=None,chkvals=False,chktol=_CHKTOL,
    temp0=None,dliq0=None,tpot0=None,dlpot0=None,chkbnd=False,
    useext=False,mathargs=None):
    """Calculate seawater isentropic haline contraction coefficient.
    
    Calculate the haline contraction coefficient of seawater at constant
    potential temperature under the given conditions, using any of
    in-situ/potential enthalpy/temperature as the thermal variable.
    
    :arg float salt: Salinity in kg/kg.
    :arg float pres: In-situ pressure in Pa.
    :arg float ppot: Potential pressure in Pa.
    :arg enth: In-situ enthalpy in J/kg. If unknown, pass None (default)
        and it will be calculated.
    :type enth: float or None
    :arg temp: In-situ temperature in K. If unknown, pass None (default)
        and it will be calculated.
    :type temp: float or None
    :arg dliq: In-situ liquid water density in kg/m3. If unknown, pass
        None (default) and it will be calculated.
    :type dliq: float or None
    :arg hpot: Potential enthalpy in J/kg. If unknown, pass None
        (default) and it will be calculated.
    :type hpot: float or None
    :arg tpot: Potential temperature in K. If unknown, pass None
        (default) and it will be calculated.
    :type tpot: float or None
    :arg dlpot: Potential liquid water density in kg/m3. If unknown,
        pass None (default) and it will be calculated.
    :type dlpot: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg temp0: Initial guess for the in-situ temperature in K. If None
        (default) then `_approx_shp` is used.
    :type temp0: float or None
    :arg dliq0: Initial guess for the in-situ liquid water density in
        kg/m3. If None (default) then `_approx_shp` is used.
    :type dliq0: float or None
    :arg tpot0: Initial guess for the potential temperature in K. If
        None (default) then `sea3b._approx_sep` is used.
    :type tpot0: float or None
    :arg dlpot0: Initial guess for the potential liquid water density in
        kg/m3. If None (default) then `sea3b._approx_sep` is used.
    :type dlpot0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg bool useext: If False (default) then the salt contribution is
        calculated from _GSCOEFFS; if True, from _GSCOEFFS_EXT.
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Contraction coefficient in 1/(kg/kg).
    :raises ValueError: If insufficient or incompatible values are
        given.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> contraction_theta(0.035,1e7,1e5,enth=1e5)
    0.728499505688
    >>> contraction_theta(0.035,1e7,1e5,hpot=1e5)
    0.726099732703
    >>> contraction_theta(0.035,1e7,1e5,temp=300.)
    0.724667977117
    >>> contraction_theta(0.035,1e7,1e5,tpot=300.)
    0.724468894946
    """
    enth, temp, dliq, hpot, tpot, dlpot = eq_pot(salt,pres,ppot,enth=enth,
        temp=temp,dliq=dliq,hpot=hpot,tpot=tpot,dlpot=dlpot,chkvals=chkvals,
        chktol=chktol,temp0=temp0,dliq0=dliq0,tpot0=tpot0,dlpot0=dlpot0,
        chkbnd=chkbnd,useext=useext,mathargs=mathargs)
    beta = sea3b.contraction_theta(salt,pres,ppot,temp=temp,dliq=dliq,tpot=tpot,
        dlpot=dlpot,useext=useext)
    return beta

def expansion_h(salt,pres,ppot,enth=None,temp=None,dliq=None,hpot=None,
    tpot=None,dlpot=None,chkvals=False,chktol=_CHKTOL,temp0=None,
    dliq0=None,tpot0=None,dlpot0=None,chkbnd=False,useext=False,
    mathargs=None):
    """Calculate seawater enthalpy expansion coefficient.
    
    Calculate the expansion coefficient of seawater with respect to
    potential enthalpy under the given conditions, using any of
    in-situ/potential enthalpy/temperature as the thermal variable.
    
    :arg float salt: Salinity in kg/kg.
    :arg float pres: In-situ pressure in Pa.
    :arg float ppot: Potential pressure in Pa.
    :arg enth: In-situ enthalpy in J/kg. If unknown, pass None (default)
        and it will be calculated.
    :type enth: float or None
    :arg temp: In-situ temperature in K. If unknown, pass None (default)
        and it will be calculated.
    :type temp: float or None
    :arg dliq: In-situ liquid water density in kg/m3. If unknown, pass
        None (default) and it will be calculated.
    :type dliq: float or None
    :arg hpot: Potential enthalpy in J/kg. If unknown, pass None
        (default) and it will be calculated.
    :type hpot: float or None
    :arg tpot: Potential temperature in K. If unknown, pass None
        (default) and it will be calculated.
    :type tpot: float or None
    :arg dlpot: Potential liquid water density in kg/m3. If unknown,
        pass None (default) and it will be calculated.
    :type dlpot: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg temp0: Initial guess for the in-situ temperature in K. If None
        (default) then `_approx_shp` is used.
    :type temp0: float or None
    :arg dliq0: Initial guess for the in-situ liquid water density in
        kg/m3. If None (default) then `_approx_shp` is used.
    :type dliq0: float or None
    :arg tpot0: Initial guess for the potential temperature in K. If
        None (default) then `sea3b._approx_sep` is used.
    :type tpot0: float or None
    :arg dlpot0: Initial guess for the potential liquid water density in
        kg/m3. If None (default) then `sea3b._approx_sep` is used.
    :type dlpot0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg bool useext: If False (default) then the salt contribution is
        calculated from _GSCOEFFS; if True, from _GSCOEFFS_EXT.
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Expansion coefficient in 1/(J/kg).
    :raises ValueError: If insufficient or incompatible values are
        given.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> expansion_h(0.035,1e7,1e5,enth=1e5)
    7.28514646021e-8
    >>> expansion_h(0.035,1e7,1e5,hpot=1e5)
    7.72876772245e-8
    >>> expansion_h(0.035,1e7,1e5,temp=300.)
    8.01009066333e-8
    >>> expansion_h(0.035,1e7,1e5,tpot=300.)
    8.05023387611e-8
    """
    enth, temp, dliq, hpot, tpot, dlpot = eq_pot(salt,pres,ppot,enth=enth,
        temp=temp,dliq=dliq,hpot=hpot,tpot=tpot,dlpot=dlpot,chkvals=chkvals,
        chktol=chktol,temp0=temp0,dliq0=dliq0,tpot0=tpot0,dlpot0=dlpot0,
        chkbnd=chkbnd,useext=useext,mathargs=mathargs)
    alpha = sea3b.expansion_h(salt,pres,ppot,temp=temp,dliq=dliq,tpot=tpot,
        dlpot=dlpot,useext=useext)
    return alpha

def expansion_theta(salt,pres,ppot,enth=None,temp=None,dliq=None,
    hpot=None,tpot=None,dlpot=None,chkvals=False,chktol=_CHKTOL,
    temp0=None,dliq0=None,tpot0=None,dlpot0=None,chkbnd=False,
    useext=False,mathargs=None):
    """Calculate seawater entropy expansion coefficient.
    
    Calculate the expansion coefficient of seawater with respect to
    potential temperature (entropy) under the given conditions, using
    any of in-situ/potential enthalpy/temperature as the thermal
    variable.
    
    :arg float salt: Salinity in kg/kg.
    :arg float pres: In-situ pressure in Pa.
    :arg float ppot: Potential pressure in Pa.
    :arg enth: In-situ enthalpy in J/kg. If unknown, pass None (default)
        and it will be calculated.
    :type enth: float or None
    :arg temp: In-situ temperature in K. If unknown, pass None (default)
        and it will be calculated.
    :type temp: float or None
    :arg dliq: In-situ liquid water density in kg/m3. If unknown, pass
        None (default) and it will be calculated.
    :type dliq: float or None
    :arg hpot: Potential enthalpy in J/kg. If unknown, pass None
        (default) and it will be calculated.
    :type hpot: float or None
    :arg tpot: Potential temperature in K. If unknown, pass None
        (default) and it will be calculated.
    :type tpot: float or None
    :arg dlpot: Potential liquid water density in kg/m3. If unknown,
        pass None (default) and it will be calculated.
    :type dlpot: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg temp0: Initial guess for the in-situ temperature in K. If None
        (default) then `_approx_shp` is used.
    :type temp0: float or None
    :arg dliq0: Initial guess for the in-situ liquid water density in
        kg/m3. If None (default) then `_approx_shp` is used.
    :type dliq0: float or None
    :arg tpot0: Initial guess for the potential temperature in K. If
        None (default) then `sea3b._approx_sep` is used.
    :type tpot0: float or None
    :arg dlpot0: Initial guess for the potential liquid water density in
        kg/m3. If None (default) then `sea3b._approx_sep` is used.
    :type dlpot0: float or None
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
    :raises ValueError: If insufficient or incompatible values are
        given.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> expansion_h(0.035,1e7,1e5,enth=1e5)
    2.91293979902e-4
    >>> expansion_h(0.035,1e7,1e5,hpot=1e5)
    3.09134848554e-4
    >>> expansion_h(0.035,1e7,1e5,temp=300.)
    3.20454167783e-4
    >>> expansion_h(0.035,1e7,1e5,tpot=300.)
    3.22069710839e-4
    """
    enth, temp, dliq, hpot, tpot, dlpot = eq_pot(salt,pres,ppot,enth=enth,
        temp=temp,dliq=dliq,hpot=hpot,tpot=tpot,dlpot=dlpot,chkvals=chkvals,
        chktol=chktol,temp0=temp0,dliq0=dliq0,tpot0=tpot0,dlpot0=dlpot0,
        chkbnd=chkbnd,useext=useext,mathargs=mathargs)
    alpha = sea3b.expansion_theta(salt,pres,ppot,temp=temp,dliq=dliq,tpot=tpot,
        dlpot=dlpot,useext=useext)
    return alpha

