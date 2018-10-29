"""Seawater-pure water vapour equilibrium functions.

This module provides thermodynamic properties for sea water in
equilibrium with pure water vapour (seawater-vapour), e.g. the enthalpy
of evaporation. It also provides a Gibbs free energy function for
seawater-vapour parcels, with primary variables being the total salinity
(mass of salt per mass of salt, liquid water and water vapour),
temperature, and pressure.

:Examples:

>>> vapourpressure(0.035,274.)
638.044692615
>>> boilingtemperature(0.035,640.)
274.042417608
>>> brinesalinity(274.,640.)
2.94396298289e-02
>>> enthalpyevap(salt=0.035,temp=274.)
2498395.40101
>>> volumeevap(salt=0.035,temp=274.)
198.075461154
>>> enthalpyevap(salt=0.035,pres=640.)
2498295.32717
>>> volumeevap(salt=0.035,pres=640.)
197.500648674
>>> enthalpyevap(temp=274.,pres=640.)
2498551.19875
>>> volumeevap(temp=274.,pres=640.)
197.469911653
>>> seavap_g(0,0,0,0.035,274.,610.)
-2748.82963245
>>> seavap_g(0,0,1,0.035,274.,610.)
137.534028399
>>> seavap_g(0,1,1,0.035,274.,610.)
63.10933482
>>> cp(0.035,274.,610.)
756270.4316
>>> expansion(0.035,274.,610.)
0.4588634213
>>> kappa_t(0.035,274.,610.)
1.19990585451e-2
>>> brinefraction(0.0035,274.,640.)
0.118887364425

:Functions:

* :func:`eq_stp`: Calculate primary variables for seawater-vapour at any
  two of the seawater salinity, temperature, and pressure.
* :func:`pressure`: Seawater-vapour pressure.
* :func:`salinity`: Seawater-vapour seawater salinity.
* :func:`temperature`: Seawater-vapour temperature.
* :func:`densitysea`: Seawater-vapour seawater density.
* :func:`densityvap`: Seawater-vapour vapour density.
* :func:`enthalpysea`: Seawater-vapour seawater enthalpy.
* :func:`enthalpyvap`: Seawater-vapour vapour enthalpy.
* :func:`entropysea`: Seawater-vapour seawater entropy.
* :func:`entropyvap`: Seawater-vapour vapour entropy.
* :func:`enthalpyevap`: Enthalpy of evaporation.
* :func:`volumeevap`: Specific volume of evaporation.
* :func:`boilingtemperature`: Seawater boiling temperature.
* :func:`brinesalinity`: Seawater salinity in equilibrium with vapour.
* :func:`vapourpressure`: Vapour pressure over seawater.
* :func:`eq_seavap`: Calculate primary variables for a seawater-vapour
  parcel at the given total salinity, temperature, and pressure.
* :func:`seavap_g`: Seawater-vapour Gibbs free energy with derivatives.
* :func:`cp`: Seawater-vapour isobaric heat capacity.
* :func:`density`: Seawater-vapour total density.
* :func:`enthalpy`: Seawater-vapour specific enthalpy.
* :func:`entropy`: Seawater-vapour specific entropy.
* :func:`expansion`: Seawater-vapour thermal expansion coefficient.
* :func:`kappa_t`: Seawater-vapour isothermal compressibility.
* :func:`brinefraction`: Seawater-vapour seawater salinity.

"""

__all__ = ['eq_stp','pressure','salinity','temperature','densitysea',
    'densityvap','enthalpysea','enthalpyvap','entropysea','entropyvap',
    'enthalpyevap','volumeevap',
    'boilingtemperature','brinesalinity','vapourpressure',
    'eq_seavap','seavap_g','cp','density','enthalpy','entropy','expansion',
    'kappa_t','brinefraction']

import warnings
import numpy
import constants0
import flu1
import flu2
import sal2
import maths3
import flu3a
import sea3a
import maths4

_CHKTOL = constants0.CHKTOL
_RWAT = constants0.RWAT
_MWAT = constants0.MWAT
_MSAL = constants0.MSAL
_TTP = constants0.TTP
_PTPE = constants0.PTPE
_LLVTP = constants0.LLVTP
_CLIQ = constants0.CLIQ
_CVAP = constants0.CVAP
_EPSS = _MWAT/_MSAL
_AVL = _LLVTP/(_RWAT*_TTP)
_BVL = (_CLIQ-_CVAP)/_RWAT
_chkflubnds = constants0.chkflubnds
_chksalbnds = constants0.chksalbnds
_flu_f = flu1.flu_f
_eq_chempot = flu2.eq_chempot
_eq_pressure = flu2.eq_pressure
_sal_g = sal2.sal_g
_eq_liqpot = sal2.eq_liqpot
_newton = maths3.newton
_dliq_default = flu3a._dliq_default
_dvap_default = flu3a._dvap_default


## Equilibrium functions
def _approx_st(salt,temp):
    """Approximate PDlDv at ST.
    
    Approximate the pressure, seawater liquid water density, and water
    vapour density for seawater in equilibrium with pure water vapour at
    the given salinity and temperature.
    
    :arg float salt: Salinity in kg/kg.
    :arg float temp: Temperature in K.
    :returns: Pressure, liquid water density, and water vapour density
        (all in SI units).
    """
    earg = _AVL*(1 - _TTP/temp)
    earg += _BVL*(1 - _TTP/temp - numpy.log(temp/_TTP))
    earg += -_EPSS*salt
    pres = _PTPE * numpy.exp(earg)
    dliq = _dliq_default(temp,pres)
    dvap = _dvap_default(temp,pres)
    return pres, dliq, dvap

def _approx_sp(salt,pres):
    """Approximate TDlDv at SP.
    
    Approximate the temperature, seawater liquid water density, and
    water vapour density for seawater in equilibrium with pure water
    vapour at the given salinity and pressure.
    
    :arg float salt: Salinity in kg/kg.
    :arg float pres: Pressure in Pa.
    :returns: Temperature, liquid water density, and water vapour
        density (all in SI units).
    """
    r = _AVL/_BVL
    v = (_EPSS*salt + numpy.log(pres/_PTPE))/_BVL
    x = maths4.lamb2(v,r)
    temp = _TTP/x
    dliq = _dliq_default(temp,pres)
    dvap = _dvap_default(temp,pres)
    return temp, dliq, dvap

def _approx_tp(temp,pres,dliq,dvap):
    """Approximate S at TP.
    
    Approximate the salinity, seawater liquid water density, and water
    vapour density for seawater in equilibrium with pure water vapour at
    the given temperature and pressure.
    
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg float dliq: Seawater liquid water density in kg/m3 (unused).
    :arg float dvap: Water vapour density in kg/m3 (unused).
    :returns: Salinity in kg/kg.
    """
    salt = _AVL*(1 - _TTP/temp)
    salt += _BVL*(1 - _TTP/temp - numpy.log(temp/_TTP))
    salt += -numpy.log(pres/_PTPE)
    salt /= _EPSS
    salt = max(salt,_CHKTOL)
    return salt

def _diff_st(p,dl,dv,salt,temp,useext=False):
    """Calculate seawater-vapour disequilibrium at ST.
    
    Calculate both sides of the equations
    
        given pressure = pressure in seawater
        given pressure = pressure in water vapour
        chemical potential of water in seawater = potential of vapour
    
    and their Jacobians with respect to pressure, seawater liquid water
    density, and water vapour density. Solving these equations gives
    equilibrium values at the given salinity and temperature.
    
    :arg float p: Pressure in Pa.
    :arg float dl: Seawater liquid water density in kg/m3.
    :arg float dv: Water vapour density in kg/m3.
    :arg float salt: Salinity in kg/kg.
    :arg float temp: Temperature in K.
    :arg bool useext: If False (default) then the salt contribution is
        calculated from _GSCOEFFS; if True, from _GSCOEFFS_EXT.
    :returns: Left-hand side of the equation, right-hand side,
        Jacobian of LHS, and Jacobian of RHS.
    :rtype: tuple(array(float))
    """
    pls = _eq_pressure(0,0,temp,dl)
    pv = _eq_pressure(0,0,temp,dv)
    muls = _eq_chempot(0,0,temp,dl)
    muls += _eq_liqpot(0,0,0,salt,temp,p,useext=useext)
    muv = _eq_chempot(0,0,temp,dv)
    lhs = numpy.array([p, p, muls])
    rhs = numpy.array([pls, pv, muv])
    
    pls_d = _eq_pressure(0,1,temp,dl)
    pv_d = _eq_pressure(0,1,temp,dv)
    muls_p = _eq_liqpot(0,0,1,salt,temp,p,useext=useext)
    muls_d = _eq_chempot(0,1,temp,dl)
    muv_d = _eq_chempot(0,1,temp,dv)
    dlhs = numpy.array([[1.,0.,0.], [1.,0.,0.], [muls_p,muls_d,0.]])
    drhs = numpy.array([[0.,pls_d,0.], [0.,0.,pv_d], [0.,0.,muv_d]])
    return lhs, rhs, dlhs, drhs

def _diff_sp(t,dl,dv,salt,pres,useext=False):
    """Calculate seawater-vapour disequilibrium at SP.
    
    Calculate both sides of the equations
    
        given pressure = pressure in seawater
        given pressure = pressure in water vapour
        chemical potential of water in seawater = potential of vapour
    
    and their Jacobians with respect to temperature, seawater liquid
    water density, and water vapour density. Solving these equations
    gives equilibrium values at the given salinity and pressure.

    :arg float t: Temperature in K.
    :arg float dl: Seawater liquid water density in kg/m3.
    :arg float dv: Water vapour density in kg/m3.
    :arg float salt: Salinity in kg/kg.
    :arg float pres: Pressure in Pa.
    :arg bool useext: If False (default) then the salt contribution is
        calculated from _GSCOEFFS; if True, from _GSCOEFFS_EXT.
    :returns: Left-hand side of the equation, right-hand side,
        Jacobian of LHS, and Jacobian of RHS.
    :rtype: tuple(array(float))
    """
    pls = _eq_pressure(0,0,t,dl)
    pv = _eq_pressure(0,0,t,dv)
    muls = _eq_chempot(0,0,t,dl)
    muls += _eq_liqpot(0,0,0,salt,t,pres,useext=useext)
    muv = _eq_chempot(0,0,t,dv)
    lhs = numpy.array([pres, pres, muls])
    rhs = numpy.array([pls, pv, muv])
    
    pls_t = _eq_pressure(1,0,t,dl)
    pls_d = _eq_pressure(0,1,t,dl)
    pv_t = _eq_pressure(1,0,t,dv)
    pv_d = _eq_pressure(0,1,t,dv)
    muls_t = _eq_chempot(1,0,t,dl)
    muls_t += _eq_liqpot(0,1,0,salt,t,pres,useext=useext)
    muls_d = _eq_chempot(0,1,t,dl)
    muv_t = _eq_chempot(1,0,t,dv)
    muv_d = _eq_chempot(0,1,t,dv)
    dlhs = numpy.array([[0.,0.,0.], [0.,0.,0.], [muls_t,muls_d,0.]])
    drhs = numpy.array([[pls_t,pls_d,0.], [pv_t,0.,pv_d], [muv_t,0.,muv_d]])
    return lhs, rhs, dlhs, drhs

def _diff_tp(s,temp,pres,dliq,dvap,useext=False):
    """Calculate seawater-vapour disequilibrium at TP.
    
    Calculate both sides of the equation
    
        chemical potential of water in seawater = potential of vapour
    
    and their derivatives with respect to salinity. Solving this
    equation gives the equilibrium salinity at the given temperature and
    pressure.

    :arg float s: Salinity in kg/kg.
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg float dliq: Seawater liquid water density in kg/m3.
    :arg float dvap: Water vapour density in kg/m3.
    :arg bool useext: If False (default) then the salt contribution is
        calculated from _GSCOEFFS; if True, from _GSCOEFFS_EXT.
    :returns: Left-hand side of the equation, right-hand side,
        derivative of LHS, and derivative of RHS.
    :rtype: tuple(float)
    """
    muls = _eq_chempot(0,0,temp,dliq)
    muls += _eq_liqpot(0,0,0,s,temp,pres,useext=useext)
    muv = _eq_chempot(0,0,temp,dvap)
    lhs = muls
    rhs = muv
    
    muls_s = _eq_liqpot(1,0,0,s,temp,pres,useext=useext)
    dlhs = muls_s
    drhs = 0.
    return lhs, rhs, dlhs, drhs

def eq_stp(salt=None,temp=None,pres=None,dliq=None,dvap=None,
    chkvals=False,chktol=_CHKTOL,salt0=None,temp0=None,pres0=None,
    dliq0=None,dvap0=None,chkbnd=False,useext=False,mathargs=None):
    """Get primary seawater-vapour variables at STP.
    
    Get the values of all primary variables for seawater and pure water
    vapour in equilibrium. Any two of the salinity, temperature, and
    pressure must be given.
    
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
        (default) then `_approx_tp` is used.
    :type salt0: float or None
    :arg temp0: Initial guess for the temperature in K. If None
        (default) then `_approx_sp` is used.
    :type temp0: float or None
    :arg pres0: Initial guess for the pressure in Pa. If None (default)
        then `_approx_st` is used.
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
        dvap = flu3a.eq_tp_vap(temp,pres,dvap=dvap,dvap0=dvap0,
            mathargs=mathargs)
        fargs = (temp,pres,dliq,dvap)
        salt = _newton(_diff_tp,salt0,_approx_tp,fargs=fargs,fkwargs=fkwargs,
            **mathargs)
    elif temp is None:
        x0 = (temp0,dliq0,dvap0)
        fargs = (salt,pres)
        x1 = _newton(_diff_sp,x0,_approx_sp,fargs=fargs,fkwargs=fkwargs,
            **mathargs)
        temp, dliq, dvap = x1
    elif pres is None:
        x0 = (pres0,dliq0,dvap0)
        fargs = (salt,temp)
        x1 = _newton(_diff_st,x0,_approx_st,fargs=fargs,fkwargs=fkwargs,
            **mathargs)
        pres, dliq, dvap = x1
    elif dliq is None:
        dliq = flu3a.eq_tp_liq(temp,pres,dliq=dliq,dliq0=dliq0,
            mathargs=mathargs)
    elif dvap is None:
        dvap = flu3a.eq_tp_vap(temp,pres,dvap=dvap,dvap0=dvap0,
            mathargs=mathargs)
    
    _chkflubnds(temp,dliq,chkbnd=chkbnd)
    _chkflubnds(temp,dvap,chkbnd=chkbnd)
    _chksalbnds(salt,temp,pres,chkbnd=chkbnd)
    if not chkvals:
        return salt, temp, pres, dliq, dvap
    
    lhs, rhs, __, __ = _diff_st(pres,dliq,dvap,salt,temp,useext=useext)
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
def pressure(salt=None,temp=None,pres=None,dliq=None,dvap=None,
    chkvals=False,chktol=_CHKTOL,salt0=None,temp0=None,pres0=None,
    dliq0=None,dvap0=None,chkbnd=False,useext=False,mathargs=None):
    """Calculate seawater-vapour pressure.
    
    Calculate the pressure of seawater and water vapour in equilibrium.
    
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
        (default) then `_approx_tp` is used.
    :type salt0: float or None
    :arg temp0: Initial guess for the temperature in K. If None
        (default) then `_approx_sp` is used.
    :type temp0: float or None
    :arg pres0: Initial guess for the pressure in Pa. If None (default)
        then `_approx_st` is used.
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
    :raises ValueError: If fewer than two of salt, temp, and pres are
        provided.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> pressure(salt=0.035,temp=274.)
    638.044692615
    """
    salt, temp, pres, dliq, dvap = eq_stp(salt=salt,temp=temp,pres=pres,
        dliq=dliq,dvap=dvap,chkvals=chkvals,chktol=chktol,salt0=salt0,
        temp0=temp0,pres0=pres0,dliq0=dliq0,dvap0=dvap0,chkbnd=chkbnd,
        useext=useext,mathargs=mathargs)
    return pres

def salinity(salt=None,temp=None,pres=None,dliq=None,dvap=None,
    chkvals=False,chktol=_CHKTOL,salt0=None,temp0=None,pres0=None,
    dliq0=None,dvap0=None,chkbnd=False,useext=False,mathargs=None):
    """Calculate seawater-vapour salinity.
    
    Calculate the salinity of seawater and water vapour in equilibrium.
    
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
        (default) then `_approx_tp` is used.
    :type salt0: float or None
    :arg temp0: Initial guess for the temperature in K. If None
        (default) then `_approx_sp` is used.
    :type temp0: float or None
    :arg pres0: Initial guess for the pressure in Pa. If None (default)
        then `_approx_st` is used.
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
    :raises ValueError: If fewer than two of salt, temp, and pres are
        provided.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> salinity(temp=274.,pres=640.)
    2.94396298289e-02
    """
    salt, temp, pres, dliq, dvap = eq_stp(salt=salt,temp=temp,pres=pres,
        dliq=dliq,dvap=dvap,chkvals=chkvals,chktol=chktol,salt0=salt0,
        temp0=temp0,pres0=pres0,dliq0=dliq0,dvap0=dvap0,chkbnd=chkbnd,
        useext=useext,mathargs=mathargs)
    return salt

def temperature(salt=None,temp=None,pres=None,dliq=None,dvap=None,
    chkvals=False,chktol=_CHKTOL,salt0=None,temp0=None,pres0=None,
    dliq0=None,dvap0=None,chkbnd=False,useext=False,mathargs=None):
    """Calculate seawater-vapour temperature.
    
    Calculate the temperature of seawater and water vapour in
    equilibrium.
    
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
        (default) then `_approx_tp` is used.
    :type salt0: float or None
    :arg temp0: Initial guess for the temperature in K. If None
        (default) then `_approx_sp` is used.
    :type temp0: float or None
    :arg pres0: Initial guess for the pressure in Pa. If None (default)
        then `_approx_st` is used.
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
    :raises ValueError: If fewer than two of salt, temp, and pres are
        provided.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> temperature(salt=0.035,pres=640.)
    274.042417608
    """
    salt, temp, pres, dliq, dvap = eq_stp(salt=salt,temp=temp,pres=pres,
        dliq=dliq,dvap=dvap,chkvals=chkvals,chktol=chktol,salt0=salt0,
        temp0=temp0,pres0=pres0,dliq0=dliq0,dvap0=dvap0,chkbnd=chkbnd,
        useext=useext,mathargs=mathargs)
    return temp

def densitysea(salt=None,temp=None,pres=None,dliq=None,dvap=None,
    chkvals=False,chktol=_CHKTOL,salt0=None,temp0=None,pres0=None,
    dliq0=None,dvap0=None,chkbnd=False,useext=False,mathargs=None):
    """Calculate seawater-vapour seawater density.
    
    Calculate the density of seawater in equilibrium with water vapour.
    
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
        (default) then `_approx_tp` is used.
    :type salt0: float or None
    :arg temp0: Initial guess for the temperature in K. If None
        (default) then `_approx_sp` is used.
    :type temp0: float or None
    :arg pres0: Initial guess for the pressure in Pa. If None (default)
        then `_approx_st` is used.
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
    :raises ValueError: If fewer than two of salt, temp, and pres are
        provided.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> densitysea(salt=0.035,pres=640.)
    1027.87349551
    >>> densitysea(salt=0.035,temp=274.)
    1027.87626132
    >>> densitysea(temp=274.,pres=640.)
    1023.42713047
    """
    salt, temp, pres, dliq, dvap = eq_stp(salt=salt,temp=temp,pres=pres,
        dliq=dliq,dvap=dvap,chkvals=chkvals,chktol=chktol,salt0=salt0,
        temp0=temp0,pres0=pres0,dliq0=dliq0,dvap0=dvap0,chkbnd=chkbnd,
        useext=useext,mathargs=mathargs)
    dsea = sea3a.density(salt,temp,pres,dliq=dliq,useext=useext)
    return dsea

def densityvap(salt=None,temp=None,pres=None,dliq=None,dvap=None,
    chkvals=False,chktol=_CHKTOL,salt0=None,temp0=None,pres0=None,
    dliq0=None,dvap0=None,chkbnd=False,useext=False,mathargs=None):
    """Calculate seawater-vapour vapour density.
    
    Calculate the density of water vapour in equilibrium with seawater.
    
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
        (default) then `_approx_tp` is used.
    :type salt0: float or None
    :arg temp0: Initial guess for the temperature in K. If None
        (default) then `_approx_sp` is used.
    :type temp0: float or None
    :arg pres0: Initial guess for the pressure in Pa. If None (default)
        then `_approx_st` is used.
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
    :raises ValueError: If fewer than two of salt, temp, and pres are
        provided.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> densityvap(salt=0.035,pres=640.)
    5.06324888818e-03
    >>> densityvap(salt=0.035,temp=274.)
    5.04855547811e-3
    >>> densityvap(temp=274.,pres=640.)
    5.06403699513e-3
    """
    salt, temp, pres, dliq, dvap = eq_stp(salt=salt,temp=temp,pres=pres,
        dliq=dliq,dvap=dvap,chkvals=chkvals,chktol=chktol,salt0=salt0,
        temp0=temp0,pres0=pres0,dliq0=dliq0,dvap0=dvap0,chkbnd=chkbnd,
        useext=useext,mathargs=mathargs)
    return dvap

def enthalpysea(salt=None,temp=None,pres=None,dliq=None,dvap=None,
    chkvals=False,chktol=_CHKTOL,salt0=None,temp0=None,pres0=None,
    dliq0=None,dvap0=None,chkbnd=False,useext=False,mathargs=None):
    """Calculate seawater-vapour seawater enthalpy.
    
    Calculate the specific enthalpy of seawater in equilibrium with
    water vapour.
    
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
        (default) then `_approx_tp` is used.
    :type salt0: float or None
    :arg temp0: Initial guess for the temperature in K. If None
        (default) then `_approx_sp` is used.
    :type temp0: float or None
    :arg pres0: Initial guess for the pressure in Pa. If None (default)
        then `_approx_st` is used.
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
    :returns: Seawater enthalpy in J/kg.
    :raises ValueError: If fewer than two of salt, temp, and pres are
        provided.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> enthalpysea(salt=0.035,pres=640.)
    3465.122066
    >>> enthalpysea(salt=0.035,temp=274.)
    3295.966293
    >>> enthalpysea(temp=274.,pres=640.)
    3405.933537
    """
    salt, temp, pres, dliq, dvap = eq_stp(salt=salt,temp=temp,pres=pres,
        dliq=dliq,dvap=dvap,chkvals=chkvals,chktol=chktol,salt0=salt0,
        temp0=temp0,pres0=pres0,dliq0=dliq0,dvap0=dvap0,chkbnd=chkbnd,
        useext=useext,mathargs=mathargs)
    hsea = sea3a.enthalpy(salt,temp,pres,dliq=dliq,useext=useext)
    return hsea

def enthalpyvap(salt=None,temp=None,pres=None,dliq=None,dvap=None,
    chkvals=False,chktol=_CHKTOL,salt0=None,temp0=None,pres0=None,
    dliq0=None,dvap0=None,chkbnd=False,useext=False,mathargs=None):
    """Calculate seawater-vapour vapour enthalpy.
    
    Calculate the specific enthalpy of water vapour in equilibrium with
    seawater.
    
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
        (default) then `_approx_tp` is used.
    :type salt0: float or None
    :arg temp0: Initial guess for the temperature in K. If None
        (default) then `_approx_sp` is used.
    :type temp0: float or None
    :arg pres0: Initial guess for the pressure in Pa. If None (default)
        then `_approx_st` is used.
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
    :returns: Water vapour enthalpy in J/kg.
    :raises ValueError: If fewer than two of salt, temp, and pres are
        provided.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> enthalpyvap(salt=0.035,pres=640.)
    2502546.89505
    >>> enthalpyvap(salt=0.035,temp=274.)
    2502469.07187
    >>> enthalpyvap(temp=274.,pres=640.)
    2502466.96633
    """
    salt, temp, pres, dliq, dvap = eq_stp(salt=salt,temp=temp,pres=pres,
        dliq=dliq,dvap=dvap,chkvals=chkvals,chktol=chktol,salt0=salt0,
        temp0=temp0,pres0=pres0,dliq0=dliq0,dvap0=dvap0,chkbnd=chkbnd,
        useext=useext,mathargs=mathargs)
    hvap = flu2.enthalpy(temp,dvap)
    return hvap

def entropysea(salt=None,temp=None,pres=None,dliq=None,dvap=None,
    chkvals=False,chktol=_CHKTOL,salt0=None,temp0=None,pres0=None,
    dliq0=None,dvap0=None,chkbnd=False,useext=False,mathargs=None):
    """Calculate seawater-vapour seawater entropy.
    
    Calculate the specific entropy of seawater in equilibrium with water
    vapour.
    
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
        (default) then `_approx_tp` is used.
    :type salt0: float or None
    :arg temp0: Initial guess for the temperature in K. If None
        (default) then `_approx_sp` is used.
    :type temp0: float or None
    :arg pres0: Initial guess for the pressure in Pa. If None (default)
        then `_approx_st` is used.
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
    :returns: Seawater entropy in J/kg/K.
    :raises ValueError: If fewer than two of salt, temp, and pres are
        provided.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> entropysea(salt=0.035,pres=640.)
    13.06170045
    >>> entropysea(salt=0.035,temp=274.)
    12.4443983378
    >>> entropysea(temp=274.,pres=640.)
    14.0256815113
    """
    salt, temp, pres, dliq, dvap = eq_stp(salt=salt,temp=temp,pres=pres,
        dliq=dliq,dvap=dvap,chkvals=chkvals,chktol=chktol,salt0=salt0,
        temp0=temp0,pres0=pres0,dliq0=dliq0,dvap0=dvap0,chkbnd=chkbnd,
        useext=useext,mathargs=mathargs)
    ssea = sea3a.entropy(salt,temp,pres,dliq=dliq,useext=useext)
    return ssea

def entropyvap(salt=None,temp=None,pres=None,dliq=None,dvap=None,
    chkvals=False,chktol=_CHKTOL,salt0=None,temp0=None,pres0=None,
    dliq0=None,dvap0=None,chkbnd=False,useext=False,mathargs=None):
    """Calculate seawater-vapour vapour entropy.
    
    Calculate the specific entropy of water vapour in equilibrium with
    seawater.
    
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
        (default) then `_approx_tp` is used.
    :type salt0: float or None
    :arg temp0: Initial guess for the temperature in K. If None
        (default) then `_approx_sp` is used.
    :type temp0: float or None
    :arg pres0: Initial guess for the pressure in Pa. If None (default)
        then `_approx_st` is used.
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
    :returns: Water vapour entropy in J/kg/K.
    :raises ValueError: If fewer than two of salt, temp, and pres are
        provided.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> entropyvap(salt=0.035,pres=640.)
    9140.56256600
    >>> entropyvap(salt=0.035,temp=274.)
    9141.68990452
    >>> entropyvap(temp=274.,pres=640.)
    9140.27087793
    """
    salt, temp, pres, dliq, dvap = eq_stp(salt=salt,temp=temp,pres=pres,
        dliq=dliq,dvap=dvap,chkvals=chkvals,chktol=chktol,salt0=salt0,
        temp0=temp0,pres0=pres0,dliq0=dliq0,dvap0=dvap0,chkbnd=chkbnd,
        useext=useext,mathargs=mathargs)
    svap = flu2.entropy(temp,dvap)
    return svap

def enthalpyevap(salt=None,temp=None,pres=None,dliq=None,dvap=None,
    chkvals=False,chktol=_CHKTOL,salt0=None,temp0=None,pres0=None,
    dliq0=None,dvap0=None,chkbnd=False,useext=False,mathargs=None):
    """Calculate seawater-vapour enthalpy of evaporation.
    
    Calculate the enthalpy of evaporation of seawater.
    
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
        (default) then `_approx_tp` is used.
    :type salt0: float or None
    :arg temp0: Initial guess for the temperature in K. If None
        (default) then `_approx_sp` is used.
    :type temp0: float or None
    :arg pres0: Initial guess for the pressure in Pa. If None (default)
        then `_approx_st` is used.
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
    :returns: Enthalpy in J/kg.
    :raises ValueError: If fewer than two of salt, temp, and pres are
        provided.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> enthalpyevap(salt=0.035,pres=640.)
    2498295.32717
    >>> enthalpyevap(salt=0.035,temp=274.)
    2498395.40101
    >>> enthalpyevap(temp=274.,pres=640.)
    2498551.19875
    """
    salt, temp, pres, dliq, dvap = eq_stp(salt=salt,temp=temp,pres=pres,
        dliq=dliq,dvap=dvap,chkvals=chkvals,chktol=chktol,salt0=salt0,
        temp0=temp0,pres0=pres0,dliq0=dliq0,dvap0=dvap0,chkbnd=chkbnd,
        useext=useext,mathargs=mathargs)
    
    hl = flu2.enthalpy(temp,dliq)
    hs = salt*sal2.saltenthalpy(salt,temp,pres,useext=useext)
    hs_s = _sal_g(1,0,0,salt,temp,pres,useext=useext)
    hs_s -= temp*_sal_g(1,1,0,salt,temp,pres,useext=useext)
    hv = flu2.enthalpy(temp,dvap)
    hevap = hv - (hl + hs - salt*hs_s)
    return hevap

def volumeevap(salt=None,temp=None,pres=None,dliq=None,dvap=None,
    chkvals=False,chktol=_CHKTOL,salt0=None,temp0=None,pres0=None,
    dliq0=None,dvap0=None,chkbnd=False,useext=False,mathargs=None):
    """Calculate seawater-vapour volume of evaporation.
    
    Calculate the specific volume of evaporation of seawater.
    
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
        (default) then `_approx_tp` is used.
    :type salt0: float or None
    :arg temp0: Initial guess for the temperature in K. If None
        (default) then `_approx_sp` is used.
    :type temp0: float or None
    :arg pres0: Initial guess for the pressure in Pa. If None (default)
        then `_approx_st` is used.
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
    :returns: Specific volume in m3/kg.
    :raises ValueError: If fewer than two of salt, temp, and pres are
        provided.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> volumeevap(salt=0.035,pres=640.)
    197.500648674
    >>> volumeevap(salt=0.035,temp=274.)
    198.075461154
    >>> volumeevap(temp=274.,pres=640.)
    197.469911653
    """
    salt, temp, pres, dliq, dvap = eq_stp(salt=salt,temp=temp,pres=pres,
        dliq=dliq,dvap=dvap,chkvals=chkvals,chktol=chktol,salt0=salt0,
        temp0=temp0,pres0=pres0,dliq0=dliq0,dvap0=dvap0,chkbnd=chkbnd,
        useext=useext,mathargs=mathargs)
    
    vl = dliq**(-1)
    vv = dvap**(-1)
    vs = _sal_g(0,0,1,salt,temp,pres,useext=useext)
    vs_s = _sal_g(1,0,1,salt,temp,pres,useext=useext)
    vevap = vv - (vl + vs - salt*vs_s)
    return vevap


## Functions of 2 variables
def boilingtemperature(salt,pres,temp=None,dliq=None,dvap=None,
    chkvals=False,chktol=_CHKTOL,temp0=None,dliq0=None,dvap0=None,
    chkbnd=False,useext=False,mathargs=None):
    """Calculate seawater boiling temperature.
    
    Calculate the temperature at which seawater of the given salinity
    and pressure will boil to produce pure water vapour.
    
    :arg float salt: Salinity in kg/kg.
    :arg float pres: Pressure in Pa.
    :arg temp: Temperature in K. If unknown, pass None (default) and it
        will be calculated.
    :type temp: float or None
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
    :arg temp0: Initial guess for the temperature in K. If None
        (default) then `_approx_sp` is used.
    :type temp0: float or None
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
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> boilingtemperature(0.035,640.)
    274.042417608
    """
    __, temp, __, dliq, dvap = eq_stp(salt=salt,pres=pres,temp=temp,
        dliq=dliq,dvap=dvap,chkvals=chkvals,chktol=chktol,temp0=temp0,
        dliq0=dliq0,dvap0=dvap0,chkbnd=chkbnd,useext=useext,mathargs=mathargs)
    return temp

def brinesalinity(temp,pres,salt=None,dliq=None,dvap=None,
    chkvals=False,chktol=_CHKTOL,salt0=None,dliq0=None,dvap0=None,
    chkbnd=False,useext=False,mathargs=None):
    """Calculate seawater-vapour salinity.
    
    Calculate the salinity of seawater in equilibrium with water vapour
    at the given temperature and pressure.
    
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg salt: Salinity in kg/kg. If unknown, pass None (default) and it
        will be calculated.
    :type salt: float or None
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
        (default) then `_approx_tp` is used.
    :type salt0: float or None
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
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> brinesalinity(274.,640.)
    2.94396298289e-02
    """
    salt, __, __, dliq, dvap = eq_stp(temp=temp,pres=pres,salt=salt,
        dliq=dliq,dvap=dvap,chkvals=chkvals,chktol=chktol,salt0=salt0,
        dliq0=dliq0,dvap0=dvap0,chkbnd=chkbnd,useext=useext,mathargs=mathargs)
    return salt

def vapourpressure(salt,temp,pres=None,dliq=None,dvap=None,
    chkvals=False,chktol=_CHKTOL,pres0=None,dliq0=None,dvap0=None,
    chkbnd=False,useext=False,mathargs=None):
    """Calculate seawater-vapour vapour pressure.
    
    Calculate the pressure of pure water vapour over seawater of the
    given salinity and temperature.
    
    :arg float salt: Salinity in kg/kg.
    :arg float temp: Temperature in K.
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
    :arg pres0: Initial guess for the pressure in Pa. If None (default)
        then `_approx_st` is used.
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
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> vapourpressure(0.035,274.)
    638.044692615
    """
    __, __, pres, dliq, dvap = eq_stp(salt=salt,temp=temp,pres=pres,
        dliq=dliq,dvap=dvap,chkvals=chkvals,chktol=chktol,pres0=pres0,
        dliq0=dliq0,dvap0=dvap0,chkbnd=chkbnd,useext=useext,mathargs=mathargs)
    return pres


## Seawater-vapour parcel properties
def eq_seavap(svsal,temp,pres,salt=None,dliq=None,dvap=None,
    chkvals=False,chktol=_CHKTOL,salt0=None,dliq0=None,dvap0=None,
    chkbnd=False,useext=False,mathargs=None):
    """Get primary seawater-vapour variables at SsvTP.
    
    Get the values of all primary variables for a seawater-water vapour
    parcel at the given total salinity, temperature, and pressure. Total
    salinity here is the ratio of the mass of salt to the total parcel
    mass (salt + liquid water + water vapour).
    
    If the calculation has already been done, the results can be passed
    to avoid unnecessary repeat calculations. If enough values are
    passed, they will be checked for consistency if chkvals is True.
    
    :arg float svsal: Total sea-vapour salinity in kg/kg.
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg salt: Seawater salinity in kg/kg. If unknown, pass None
        (default) and it will be calculated.
    :type salt: float or None
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
    :arg salt0: Initial guess for the seawater salinity in kg/kg. If
        None (default) then `_approx_tp` is used.
    :type salt0: float or None
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
    :returns: Seawater salinity, liquid water density, and water vapour
        density (all in SI units).
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    :raises RuntimeWarning: If the equilibrium seawater salinity is
        lower than the total parcel salinity.
    """
    if any(val is None for val in (salt,dliq,dvap)):
        salt, __, __, dliq, dvap = eq_stp(temp=temp,pres=pres,salt=salt,
            dliq=dliq,dvap=dvap,chkvals=chkvals,chktol=chktol,salt0=salt0,
            dliq0=dliq0,dvap0=dvap0,chkbnd=chkbnd,useext=useext,
            mathargs=mathargs)
    if salt < svsal:
        warnmsg = ('Equilibrium salinity {0} is lower than the total parcel '
            'salinity {1}').format(salt,svsal)
        warnings.warn(warnmsg,RuntimeWarning)
        salt = svsal
    return salt, dliq, dvap

def seavap_g(drvs,drvt,drvp,svsal,temp,pres,salt=None,dliq=None,
    dvap=None,chkvals=False,chktol=_CHKTOL,salt0=None,dliq0=None,
    dvap0=None,chkbnd=False,useext=False,mathargs=None):
    """Calculate seawater-vapour Gibbs energy with derivatives.
    
    Calculate the specific Gibbs free energy of a seawater-water vapour
    parcel or its derivatives with respect to total salinity,
    temperature, and pressure.
    
    :arg int drvs: Number of total salinity derivatives.
    :arg int drvt: Number of temperature derivatives.
    :arg int drvp: Number of pressure derivatives.
    :arg float svsal: Total sea-vapour salinity in kg/kg.
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg salt: Seawater salinity in kg/kg. If unknown, pass None
        (default) and it will be calculated.
    :type salt: float or None
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
    :arg salt0: Initial guess for the seawater salinity in kg/kg. If
        None (default) then `_approx_tp` is used.
    :type salt0: float or None
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
    :returns: Gibbs free energy in units of
        (J/kg) / (kg/kg)^drvs / K^drvt / Pa^drvp.
    :raises ValueError: If any of (drvs,drvt,drvp) are negative or if
        (drvs+drvt+drvp) > 2.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    :raises RuntimeWarning: If the equilibrium seawater salinity is
        lower than the total parcel salinity.
    
    :Examples:
    
    >>> seavap_g(0,0,0,0.035,274.,610.)
    -2748.82963245
    >>> seavap_g(1,0,0,0.035,274.,610.)
    151028.257424
    >>> seavap_g(0,1,0,0.035,274.,610.)
    -6072.50817709
    >>> seavap_g(0,0,1,0.035,274.,610.)
    137.534028399
    >>> seavap_g(2,0,0,0.035,274.,610.)
    0.
    >>> seavap_g(1,1,0,0.035,274.,610.)
    88286.3861825
    >>> seavap_g(1,0,1,0.035,274.,610.)
    -1990.13848555
    >>> seavap_g(0,2,0,0.035,274.,610.)
    -2760.11106421
    >>> seavap_g(0,1,1,0.035,274.,610.)
    63.10933482
    >>> seavap_g(0,0,2,0.035,274.,610.)
    -1.65027885871
    """
    drvtup = (drvs,drvt,drvp)
    if any(drv < 0 for drv in drvtup) or sum(drvtup) > 2:
        errmsg = 'Derivatives {0} not recognized'.format(drvtup)
        raise ValueError(errmsg)
    salt, dliq, dvap = eq_seavap(svsal,temp,pres,salt=salt,dliq=dliq,
        dvap=dvap,chkvals=chkvals,chktol=chktol,salt0=salt0,dliq0=dliq0,
        dvap0=dvap0,chkbnd=chkbnd,useext=useext,mathargs=mathargs)
    seaf = svsal/salt
    
    # Straightforward derivatives
    if (drvs,drvt,drvp) == (0,0,0):
        gl = _eq_chempot(0,0,temp,dliq)
        gs = _sal_g(0,0,0,salt,temp,pres,useext=useext)
        gv = _eq_chempot(0,0,temp,dvap)
        g = seaf*(gl + gs) + (1-seaf)*gv
        return g
    elif (drvs,drvt,drvp) == (1,0,0):
        gs_s = _sal_g(1,0,0,salt,temp,pres,useext=useext)
        g_s = gs_s
        return g_s
    elif (drvs,drvt,drvp) == (0,1,0):
        fl_t = _flu_f(1,0,temp,dliq)
        gs_t = _sal_g(0,1,0,salt,temp,pres,useext=useext)
        fv_t = _flu_f(1,0,temp,dvap)
        g_t = seaf*(fl_t + gs_t) + (1-seaf)*fv_t
        return g_t
    elif (drvs,drvt,drvp) == (0,0,1):
        gs_p = _sal_g(0,0,1,salt,temp,pres,useext=useext)
        g_p = seaf*(dliq**(-1) + gs_p) + (1-seaf)/dvap
        return g_p
    elif (drvs,drvt,drvp) == (2,0,0):
        return 0.0
    elif (drvs,drvt,drvp) == (1,1,0):
        fl_t = _flu_f(1,0,temp,dliq)
        gs_t = _sal_g(0,1,0,salt,temp,pres,useext=useext)
        fv_t = _flu_f(1,0,temp,dvap)
        g_st = (fl_t + gs_t - fv_t)/salt
        return g_st
    elif (drvs,drvt,drvp) == (1,0,1):
        gs_p = _sal_g(0,0,1,salt,temp,pres,useext=useext)
        g_sp = (gs_p + dliq**(-1) - dvap**(-1))/salt
        return g_sp
    
    # Other derivatives require inversion
    cl = _eq_pressure(0,1,temp,dliq)
    cv = _eq_pressure(0,1,temp,dvap)
    gs_ss = _sal_g(2,0,0,salt,temp,pres,useext=useext)
    if drvt > 0:
        fl_t = _flu_f(1,0,temp,dliq)
        gs_t = _sal_g(0,1,0,salt,temp,pres,useext=useext)
        gs_st = _sal_g(1,1,0,salt,temp,pres,useext=useext)
        fv_t = _flu_f(1,0,temp,dvap)
        dentr = fl_t + gs_t - salt*gs_st - fv_t
    if drvp > 0:
        gs_p = _sal_g(0,0,1,salt,temp,pres,useext=useext)
        gs_sp = _sal_g(1,0,1,salt,temp,pres,useext=useext)
        dvol = dliq**(-1) + gs_p - salt*gs_sp - dvap**(-1)
        s_p = dvol / (salt*gs_ss)
        dl_p = cl**(-1)
        dv_p = cv**(-1)
    
    if (drvs,drvt,drvp) == (0,2,0):
        fl_tt = _flu_f(2,0,temp,dliq)
        fl_td = _flu_f(1,1,temp,dliq)
        gs_tt = _sal_g(0,2,0,salt,temp,pres,useext=useext)
        fv_tt = _flu_f(2,0,temp,dvap)
        fv_td = _flu_f(1,1,temp,dvap)
        
        s_t = dentr / (salt*gs_ss)
        dl_t = -dliq**2*fl_td/cl
        dv_t = -dvap**2*fv_td/cv
        gb_tt = fl_tt + fl_td*dl_t + gs_tt
        gv_tt = fv_tt + fv_td*dv_t
        g_tt = -seaf/salt*dentr*s_t + seaf*gb_tt + (1-seaf)*gv_tt
        return g_tt
    elif (drvs,drvt,drvp) == (0,1,1):
        fl_td = _flu_f(1,1,temp,dliq)
        gs_tp = _sal_g(0,1,1,salt,temp,pres,useext=useext)
        fv_td = _flu_f(1,1,temp,dvap)
        gb_tp = fl_td*dl_p + gs_tp
        gv_tp = fv_td*dv_p
        g_tp = -seaf/salt*dentr*s_p + seaf*gb_tp + (1-seaf)*gv_tp
        return g_tp
    elif (drvs,drvt,drvp) == (0,0,2):
        gs_pp = _sal_g(0,0,2,salt,temp,pres,useext=useext)
        gb_pp = -dl_p/dliq**2 + gs_pp
        gv_pp = -dv_p/dvap**2
        g_pp = -seaf/salt*dvol*s_p + seaf*gb_pp + (1-seaf)*gv_pp
        return g_pp
    
    # Should not have made it this far!
    errmsg = 'Derivatives {0} not recognized'.format(drvtup)
    raise ValueError(errmsg)

def cp(svsal,temp,pres,salt=None,dliq=None,dvap=None,chkvals=False,
    chktol=_CHKTOL,salt0=None,dliq0=None,dvap0=None,chkbnd=False,
    useext=False,mathargs=None):
    """Calculate seawater-vapour isobaric heat capacity.
    
    Calculate the isobaric heat capacity of a seawater-vapour parcel.
    
    :arg float svsal: Total sea-vapour salinity in kg/kg.
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg salt: Seawater salinity in kg/kg. If unknown, pass None
        (default) and it will be calculated.
    :type salt: float or None
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
    :arg salt0: Initial guess for the seawater salinity in kg/kg. If
        None (default) then `_approx_tp` is used.
    :type salt0: float or None
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
    :returns: Heat capacity in J/kg/K.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    :raises RuntimeWarning: If the equilibrium seawater salinity is
        lower than the total parcel salinity.
    
    :Examples:
    
    >>> cp(0.035,274.,610.)
    756270.4316
    """
    g_tt = seavap_g(0,2,0,svsal,temp,pres,salt=salt,dliq=dliq,dvap=dvap,
        chkvals=chkvals,chktol=chktol,salt0=salt0,dliq0=dliq0,dvap0=dvap0,
        chkbnd=chkbnd,useext=useext,mathargs=mathargs)
    cp = -temp*g_tt
    return cp

def density(svsal,temp,pres,salt=None,dliq=None,dvap=None,
    chkvals=False,chktol=_CHKTOL,salt0=None,dliq0=None,dvap0=None,
    chkbnd=False,useext=False,mathargs=None):
    """Calculate seawater-vapour density.
    
    Calculate the total density of a seawater-vapour parcel.
    
    :arg float svsal: Total sea-vapour salinity in kg/kg.
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg salt: Seawater salinity in kg/kg. If unknown, pass None
        (default) and it will be calculated.
    :type salt: float or None
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
    :arg salt0: Initial guess for the seawater salinity in kg/kg. If
        None (default) then `_approx_tp` is used.
    :type salt0: float or None
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
    :returns: Density in kg/m3.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    :raises RuntimeWarning: If the equilibrium seawater salinity is
        lower than the total parcel salinity.
    
    :Examples:
    
    >>> density(0.035,274.,610.)
    7.27092786882e-3
    """
    g_p = seavap_g(0,0,1,svsal,temp,pres,salt=salt,dliq=dliq,dvap=dvap,
        chkvals=chkvals,chktol=chktol,salt0=salt0,dliq0=dliq0,dvap0=dvap0,
        chkbnd=chkbnd,useext=useext,mathargs=mathargs)
    rho = g_p**(-1)
    return rho

def enthalpy(svsal,temp,pres,salt=None,dliq=None,dvap=None,
    chkvals=False,chktol=_CHKTOL,salt0=None,dliq0=None,dvap0=None,
    chkbnd=False,useext=False,mathargs=None):
    """Calculate seawater-vapour enthalpy.
    
    Calculate the specific enthalpy of a seawater-vapour parcel.
    
    :arg float svsal: Total sea-vapour salinity in kg/kg.
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg salt: Seawater salinity in kg/kg. If unknown, pass None
        (default) and it will be calculated.
    :type salt: float or None
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
    :arg salt0: Initial guess for the seawater salinity in kg/kg. If
        None (default) then `_approx_tp` is used.
    :type salt0: float or None
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
    :returns: Enthalpy in J/kg.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    :raises RuntimeWarning: If the equilibrium seawater salinity is
        lower than the total parcel salinity.
    
    :Examples:
    
    >>> enthalpy(0.035,274.,610.)
    1661118.41089
    """
    salt, dliq, dvap = eq_seavap(svsal,temp,pres,salt=salt,dliq=dliq,
        dvap=dvap,chkvals=chkvals,chktol=chktol,salt0=salt0,dliq0=dliq0,
        dvap0=dvap0,chkbnd=chkbnd,useext=useext,mathargs=mathargs)
    g = seavap_g(0,0,0,svsal,temp,pres,salt=salt,dliq=dliq,dvap=dvap,
        useext=useext)
    g_t = seavap_g(0,1,0,svsal,temp,pres,salt=salt,dliq=dliq,dvap=dvap,
        useext=useext)
    h = g - temp*g_t
    return h

def entropy(svsal,temp,pres,salt=None,dliq=None,dvap=None,
    chkvals=False,chktol=_CHKTOL,salt0=None,dliq0=None,dvap0=None,
    chkbnd=False,useext=False,mathargs=None):
    """Calculate seawater-vapour entropy.
    
    Calculate the specific entropy of a seawater-vapour parcel.
    
    :arg float svsal: Total sea-vapour salinity in kg/kg.
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg salt: Seawater salinity in kg/kg. If unknown, pass None
        (default) and it will be calculated.
    :type salt: float or None
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
    :arg salt0: Initial guess for the seawater salinity in kg/kg. If
        None (default) then `_approx_tp` is used.
    :type salt0: float or None
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
    :returns: Entropy in J/kg/K.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    :raises RuntimeWarning: If the equilibrium seawater salinity is
        lower than the total parcel salinity.
    
    :Examples:
    
    >>> entropy(0.035,274.,610.)
    6072.50817709
    """
    g_t = seavap_g(0,1,0,svsal,temp,pres,salt=salt,dliq=dliq,dvap=dvap,
        chkvals=chkvals,chktol=chktol,salt0=salt0,dliq0=dliq0,dvap0=dvap0,
        chkbnd=chkbnd,useext=useext,mathargs=mathargs)
    s = -g_t
    return s

def expansion(svsal,temp,pres,salt=None,dliq=None,dvap=None,
    chkvals=False,chktol=_CHKTOL,salt0=None,dliq0=None,dvap0=None,
    chkbnd=False,useext=False,mathargs=None):
    """Calculate seawater-vapour thermal expansion coefficient.
    
    Calculate the thermal expansion coefficient of a seawater-vapour
    parcel.
    
    :arg float svsal: Total sea-vapour salinity in kg/kg.
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg salt: Seawater salinity in kg/kg. If unknown, pass None
        (default) and it will be calculated.
    :type salt: float or None
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
    :arg salt0: Initial guess for the seawater salinity in kg/kg. If
        None (default) then `_approx_tp` is used.
    :type salt0: float or None
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
    :returns: Expansion coefficient in 1/K.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    :raises RuntimeWarning: If the equilibrium seawater salinity is
        lower than the total parcel salinity.
    
    :Examples:
    
    >>> expansion(0.035,274.,610.)
    0.4588634213
    """
    salt, dliq, dvap = eq_seavap(svsal,temp,pres,salt=salt,dliq=dliq,
        dvap=dvap,chkvals=chkvals,chktol=chktol,salt0=salt0,dliq0=dliq0,
        dvap0=dvap0,chkbnd=chkbnd,useext=useext,mathargs=mathargs)
    g_p = seavap_g(0,0,1,svsal,temp,pres,salt=salt,dliq=dliq,dvap=dvap,
        useext=useext)
    g_tp = seavap_g(0,1,1,svsal,temp,pres,salt=salt,dliq=dliq,dvap=dvap,
        useext=useext)
    alpha = g_tp / g_p
    return alpha

def kappa_t(svsal,temp,pres,salt=None,dliq=None,dvap=None,
    chkvals=False,chktol=_CHKTOL,salt0=None,dliq0=None,dvap0=None,
    chkbnd=False,useext=False,mathargs=None):
    """Calculate seawater-vapour isothermal compressibility.
    
    Calculate the isothermal compressibility of a seawater-vapour
    parcel.
    
    :arg float svsal: Total sea-vapour salinity in kg/kg.
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg salt: Seawater salinity in kg/kg. If unknown, pass None
        (default) and it will be calculated.
    :type salt: float or None
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
    :arg salt0: Initial guess for the seawater salinity in kg/kg. If
        None (default) then `_approx_tp` is used.
    :type salt0: float or None
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
    :returns: Compressibility in 1/Pa.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    :raises RuntimeWarning: If the equilibrium seawater salinity is
        lower than the total parcel salinity.
    
    :Examples:
    
    >>> kappa_t(0.035,274.,610.)
    1.19990585451e-2
    """
    salt, dliq, dvap = eq_seavap(svsal,temp,pres,salt=salt,dliq=dliq,
        dvap=dvap,chkvals=chkvals,chktol=chktol,salt0=salt0,dliq0=dliq0,
        dvap0=dvap0,chkbnd=chkbnd,useext=useext,mathargs=mathargs)
    g_p = seavap_g(0,0,1,svsal,temp,pres,salt=salt,dliq=dliq,dvap=dvap,
        useext=useext)
    g_pp = seavap_g(0,0,2,svsal,temp,pres,salt=salt,dliq=dliq,dvap=dvap,
        useext=useext)
    kappa = -g_pp / g_p
    return kappa

def brinefraction(svsal,temp,pres,salt=None,dliq=None,dvap=None,
    chkvals=False,chktol=_CHKTOL,salt0=None,dliq0=None,dvap0=None,
    chkbnd=False,useext=False,mathargs=None):
    """Calculate seawater-vapour brine fraction.
    
    Calculate the mass fraction of brine in a seawater-vapour parcel,
    the ratio of the mass of seawater (salt + liquid water) to the total
    mass (salt + liquid water + water vapour).
    
    :arg float svsal: Total sea-vapour salinity in kg/kg.
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg salt: Seawater salinity in kg/kg. If unknown, pass None
        (default) and it will be calculated.
    :type salt: float or None
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
    :arg salt0: Initial guess for the seawater salinity in kg/kg. If
        None (default) then `_approx_tp` is used.
    :type salt0: float or None
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
    :returns: Brine fraction in kg/kg.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    :raises RuntimeWarning: If the equilibrium seawater salinity is
        lower than the total parcel salinity.
    
    :Examples:
    
    >>> brinefraction(0.0035,274.,640.)
    0.118887364425
    """
    salt, dliq, dvap = eq_seavap(svsal,temp,pres,salt=salt,dliq=dliq,
        dvap=dvap,chkvals=chkvals,chktol=chktol,salt0=salt0,dliq0=dliq0,
        dvap0=dvap0,chkbnd=chkbnd,useext=useext,mathargs=mathargs)
    seaf = svsal/salt
    return seaf

