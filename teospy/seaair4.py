"""Seawater-humid air equilibrium functions.

This module provides thermodynamic properties for humid air in
equilibrium with seawater (sea-air).

:Examples:

>>> massfractionair(0.035,300.,1e5)
0.978029483888
>>> vapourpressure(0.035,300.,1e5)
3485.92986681
>>> condensetemp(0.035,0.99,1e5)
287.367451766
>>> densityair(salt=0.035,temp=300.,pres=1e5)
1.14642944448
>>> densityvap(salt=0.035,temp=300.,pres=1e5)
2.51876465812e-2
>>> entropyair(salt=0.035,temp=300.,pres=1e5)
293.150672957
>>> enthalpyevap(salt=0.035,temp=300.,pres=1e5)
2434570.55229
>>> densityair(salt=0.035,airf=0.99,pres=1e5)
1.20553990576
>>> densityvap(salt=0.035,airf=0.99,pres=1e5)
1.20553990576e-02
>>> enthalpyevap(salt=0.035,airf=0.99,pres=1e5)
2464765.77588
>>> chempotevap(0.9,0.035,300.,1e5)
1.45584069071

:Functions:

* :func:`eq_satp`: Calculate primary variables for sea-air at any three
  of the seawater salinity, humid air dry fraction, temperature, and
  pressure.
* :func:`densityair`: Sea-air humid air density.
* :func:`densityvap`: Sea-air water vapour density.
* :func:`entropyair`: Sea-air humid air entropy.
* :func:`enthalpyevap`: Enthalpy of evaporation for sea air.
* :func:`massfractionair`: Sea-air humid air dry fraction.
* :func:`vapourpressure`: Sea-air vapour pressure.
* :func:`condensetemp`: Sea-air condensation temperature.
* :func:`chempotevap`: Scaled chemical potential difference between
  seawater and water vapour in sea-air; relative Onsager force.

"""

__all__ = ['eq_satp','densityair','densityvap','entropyair','enthalpyevap',
    'massfractionair','vapourpressure','condensetemp','chempotevap']

import numpy
import warnings
import constants0
import convert0
import flu1
import air2
import flu2
import sal2
import maths3
import air3a
import air3b
import flu3a
import sea3a
import maths4

_CHKTOL = constants0.CHKTOL
_MWAT = constants0.MWAT
_MDRY = constants0.MDRY
_MSAL = constants0.MSAL
_RUNIV = constants0.RUNIV
_RWAT = constants0.RWAT
_RDRY = constants0.RDRY
_TTP = constants0.TTP
_PTPE = constants0.PTPE
_LLVTP = constants0.LLVTP
_CLIQ = constants0.CLIQ
_CVAP = constants0.CVAP
_RSAL = _RUNIV/_MSAL
_EPSW = _MWAT/_MDRY
_EPSS = _MWAT/_MSAL
_AVL = _LLVTP / (_RWAT*_TTP)
_BVL = (_CLIQ-_CVAP)/_RWAT
_chkflubnds = constants0.chkflubnds
_chkhumbnds = constants0.chkhumbnds
_chksalbnds = constants0.chksalbnds
_flu_f = flu1.flu_f
_air_f = air2.air_f
_air_eq_pressure = air2.eq_pressure
_eq_vappot = air2.eq_vappot
_eq_chempot = flu2.eq_chempot
_flu_eq_pressure = flu2.eq_pressure
_sal_g = sal2.sal_g
_eq_liqpot = sal2.eq_liqpot
_newton = maths3.newton
_dliq_default = flu3a._dliq_default


## Equilibrium functions
def _approx_stp(salt,temp,pres,dliq):
    """Approximate ADh at STP.
    
    Approximate the humid air dry fraction and density for sea-air at
    the given salinity, temperature, and pressure.
    
    :arg float salt: Salinity in kg/kg.
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg float dliq: Seawater liquid water density in kg/m3 (unused).
    :returns: Humid air dry fraction and density (both in SI units).
    """
    earg = (_AVL*(1-_TTP/temp) + _BVL*(1-_TTP/temp+numpy.log(_TTP/temp))
        - _EPSS*salt)
    pvap = _PTPE*numpy.exp(earg)
    airf = (pres-pvap) / (pres-pvap + _EPSW*pvap)
    dhum = pres/(_RDRY*temp) / (airf+(1-airf)/_EPSW)
    return airf, dhum

def _approx_sap(salt,airf,pres):
    """Approximate TDhDl at SAP.
    
    Approximate the temperature, humid air density, and seawater liquid
    water density for sea-air at the given salinity, humid air dry
    fraction, and pressure.
    
    :arg float salt: Salinity in kg/kg.
    :arg float airf: Humid air dry fraction in kg/kg.
    :arg float pres: Pressure in Pa.
    :returns: Temperature, humid air density, and liquid water density
        (all in SI units).
    """
    pvap = pres * (1-airf)/(_EPSW*airf + 1-airf)
    r = _AVL/_BVL
    v = (_EPSS*salt + numpy.log(pvap/_PTPE))/_BVL
    x = maths4.lamb2(v,r)
    temp = _TTP/x
    dhum = pres/(_RDRY*temp) / (airf + (1-airf)/_EPSW)
    dliq = _dliq_default(temp,pres)
    return temp, dhum, dliq

def _approx_sat(salt,airf,temp):
    """Approximate PDhDl at SAT.
    
    Approximate the pressure, humid air density, and liquid water
    density of sea-air with the given salinity, humid air dry fraction,
    and temperature.
    
    :arg float salt: Salinity in kg/kg.
    :arg float airf: Humid air dry fraction in kg/kg.
    :arg float temp: Temperature in K.
    :returns: Pressure, humid air density, and liquid water density (all
        in SI units).
    """
    earg = (_AVL*(1-_TTP/temp) + _BVL*(1-_TTP/temp+numpy.log(_TTP/temp))
        - _EPSS*salt)
    pvap = _PTPE*numpy.exp(earg)
    pres = pvap * (_EPSW*airf + 1-airf)/(1-airf)
    dhum = pres/(_RDRY*temp) / (airf + (1-airf)/_EPSW)
    dliq = _dliq_default(temp,pres)
    return pres, dhum, dliq

def _approx_atp(airf,temp,pres,dhum,dliq):
    """Approximate S at ATP.
    
    Approximate the salinity of sea-air with the given humid air dry fraction, temperature, and pressure.
    
    :arg float airf: Humid air dry fraction in kg/kg.
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg float dhum: Humid air density in kg/m3 (unused).
    :arg float dliq: Seawater liquid water density in kg/m3 (unused).
    """
    pvap = pres * (1-airf)/(_EPSW*airf + 1-airf)
    dmu = ((_CLIQ-_CVAP)*(temp-_TTP-temp*numpy.log(temp/_TTP))
        - _RWAT*temp*numpy.log(pvap/_PTPE))
    salt = dmu/(_RSAL*temp)
    return salt

def _diff_stp(a,dh,salt,temp,pres,dliq,useext=False):
    """Calculate sea-air disequilibrium at STP.
    
    Calculate both sides of the equations
    
        given pressure = pressure in humid air
        chemical potential of liquid water = potential of water vapour
    
    and their Jacobians with respect to humid air dry fraction and
    density. Solving these equations gives equilibrium values at the
    given salinity, temperature, and pressure.
    
    :arg float a: Humid air dry fraction in kg/kg.
    :arg float dh: Humid air density in kg/m3.
    :arg float salt: Salinity in kg/kg.
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg float dliq: Seawater liquid water density in kg/m3.
    :arg bool useext: If False (default) then the salt contribution is
        calculated from _GSCOEFFS; if True, from _GSCOEFFS_EXT.
    :returns: Left-hand side of the equation, right-hand side,
        Jacobian of LHS, and Jacobian of RHS.
    :rtype: tuple(array(float))
    """
    ph = _air_eq_pressure(0,0,0,a,temp,dh)
    gl = _eq_chempot(0,0,temp,dliq)
    gl += _eq_liqpot(0,0,0,salt,temp,pres,useext=useext)
    gv = _eq_vappot(0,0,0,a,temp,dh)
    lhs = numpy.array([pres, gl])
    rhs = numpy.array([ph, gv])
    
    ph_a = _air_eq_pressure(1,0,0,a,temp,dh)
    ph_d = _air_eq_pressure(0,0,1,a,temp,dh)
    gv_a = _eq_vappot(1,0,0,a,temp,dh)
    gv_d = _eq_vappot(0,0,1,a,temp,dh)
    dlhs = numpy.array([[0.,0.], [0.,0.]])
    drhs = numpy.array([[ph_a,ph_d], [gv_a,gv_d]])
    return lhs, rhs, dlhs, drhs

def _diff_sap(t,dh,dl,salt,airf,pres,useext=False):
    """Calculate sea-air disequilibrium at SAP.
    
    Calculate both sides of the equations
    
        given pressure = pressure in humid air
        given pressure = pressure of liquid water
        chemical potential of liquid water = potential of water vapour
    
    and their Jacobians with respect to temperature, humid air density,
    and liquid water density. Solving these equations gives equilibrium
    values at the given salinity, humid air dry fraction, and pressure.

    :arg float t: Temperature in K.
    :arg float dh: Humid air density in kg/m3.
    :arg float dl: Seawater liquid water density in kg/m3.
    :arg float salt: Salinity in kg/kg.
    :arg float airf: Humid air dry fraction in kg/kg.
    :arg float pres: Pressure in Pa.
    :arg bool useext: If False (default) then the salt contribution is
        calculated from _GSCOEFFS; if True, from _GSCOEFFS_EXT.
    :returns: Left-hand side of the equation, right-hand side,
        Jacobian of LHS, and Jacobian of RHS.
    :rtype: tuple(array(float))
    """
    ph = _air_eq_pressure(0,0,0,airf,t,dh)
    pl = _flu_eq_pressure(0,0,t,dl)
    gl = _eq_chempot(0,0,t,dl)
    gl += _eq_liqpot(0,0,0,salt,t,pres,useext=useext)
    gv = _eq_vappot(0,0,0,airf,t,dh)
    lhs = numpy.array([pres, pres, gl])
    rhs = numpy.array([ph, pl, gv])
    
    ph_t = _air_eq_pressure(0,1,0,airf,t,dh)
    ph_d = _air_eq_pressure(0,0,1,airf,t,dh)
    pl_t = _flu_eq_pressure(1,0,t,dl)
    pl_d = _flu_eq_pressure(0,1,t,dl)
    gl_t = _eq_chempot(1,0,t,dl)
    gl_t += _eq_liqpot(0,1,0,salt,t,pres,useext=useext)
    gl_d = _eq_chempot(0,1,t,dl)
    gv_t = _eq_vappot(0,1,0,airf,t,dh)
    gv_d = _eq_vappot(0,0,1,airf,t,dh)
    dlhs = numpy.array([[0.,0.,0.], [0.,0.,0.], [gl_t,0.,gl_d]])
    drhs = numpy.array([[ph_t,ph_d,0.], [pl_t,0.,pl_d], [gv_t,gv_d,0.]])
    return lhs, rhs, dlhs, drhs

def _diff_sat(p,dh,dl,salt,airf,temp,useext=False):
    """Calculate sea-air disequilibrium at SAT.
    
    Calculate both sides of the equations
    
        given pressure = pressure in humid air
        given pressure = pressure of liquid water
        chemical potential of liquid water = potential of water vapour
    
    and their Jacobians with respect to pressure, humid air density, and
    liquid water density. Solving these equations gives equilibrium
    values at the given salinity, humid air dry fraction, and
    temperature.

    :arg float p: Pressure in Pa.
    :arg float dh: Humid air density in kg/m3.
    :arg float dl: Seawater liquid water density in kg/m3.
    :arg float salt: Salinity in kg/kg.
    :arg float airf: Humid air dry fraction in kg/kg.
    :arg float temp: Temperature in K.
    :arg bool useext: If False (default) then the salt contribution is
        calculated from _GSCOEFFS; if True, from _GSCOEFFS_EXT.
    :returns: Left-hand side of the equation, right-hand side,
        Jacobian of LHS, and Jacobian of RHS.
    :rtype: tuple(array(float))
    """
    ph = _air_eq_pressure(0,0,0,airf,temp,dh)
    pl = _flu_eq_pressure(0,0,temp,dl)
    gl = _eq_chempot(0,0,temp,dl)
    gl += _eq_liqpot(0,0,0,salt,temp,p,useext=useext)
    gv = _eq_vappot(0,0,0,airf,temp,dh)
    lhs = numpy.array([p, p, gl])
    rhs = numpy.array([ph, pl, gv])
    
    ph_d = _air_eq_pressure(0,0,1,airf,temp,dh)
    pl_d = _flu_eq_pressure(0,1,temp,dl)
    gl_d = _eq_chempot(0,1,temp,dl)
    gl_p = _eq_liqpot(0,0,1,salt,temp,p,useext=useext)
    gv_d = _eq_vappot(0,0,1,airf,temp,dh)
    dlhs = numpy.array([[1.,0.,0.], [1.,0.,0.], [gl_p,0.,gl_d]])
    drhs = numpy.array([[0.,ph_d,0.], [0.,0.,pl_d], [0.,gv_d,0.]])
    return lhs, rhs, dlhs, drhs

def _diff_atp(s,airf,temp,pres,dhum,dliq,useext=False):
    """Calculate sea-air disequilibrium at ATP.
    
    Calculate both sides of the equation
    
        chemical potential of liquid water = potential of water vapour
    
    and its derivative with respect to salinity. Solving this equation
    gives the equilibrium salinity at the given humid air dry fraction,
    temperature, and pressure.

    :arg float s: Salinity in kg/kg.
    :arg float airf: Humid air dry fraction in kg/kg.
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg float dhum: Humid air density in kg/m3.
    :arg float dliq: Seawater liquid water density in kg/m3.
    :arg bool useext: If False (default) then the salt contribution is
        calculated from _GSCOEFFS; if True, from _GSCOEFFS_EXT.
    :returns: Left-hand side of the equation, right-hand side,
        derivative of LHS, and derivative of RHS.
    :rtype: tuple(float)
    """
    gl = _eq_chempot(0,0,temp,dliq)
    gl += _eq_liqpot(0,0,0,s,temp,pres,useext=useext)
    gv = _eq_vappot(0,0,0,airf,temp,dhum)
    lhs = gl
    rhs = gv
    
    gl_s = _eq_liqpot(1,0,0,s,temp,pres,useext=useext)
    dlhs = gl_s
    drhs = 0.
    return lhs, rhs, dlhs, drhs

def eq_satp(salt=None,airf=None,temp=None,pres=None,dhum=None,dliq=None,
    chkvals=False,chktol=_CHKTOL,salt0=None,airf0=None,temp0=None,
    pres0=None,dhum0=None,dliq0=None,chkbnd=False,useext=False,
    mathargs=None):
    """Get primary sea-ice-vapour variables at SATP.
    
    Get the values of all primary variables for sea-ice-vapour. Any
    three of the seawater salinity, humid air dry fraction, temperature,
    and pressure can be provided.
    
    If the calculation has already been done, the results can be passed
    to avoid unnecessary repeat calculations. If enough values are
    passed, they will be checked for consistency if chkvals is True.
    
    :arg salt: Salinity in kg/kg. If unknown, pass None (default) and it
        will be calculated.
    :type salt: float or None
    :arg airf: Humid air dry fraction in kg/kg. If unknown, pass None
        (default) and it will be calculated.
    :type airf: float or None
    :arg temp: Temperature in K. If unknown, pass None (default) and it
        will be calculated.
    :type temp: float or None
    :arg pres: Pressure in Pa. If unknown, pass None (default) and it
        will be calculated.
    :type pres: float or None
    :arg dliq: Seawater liquid water density in kg/m3. If unknown, pass
        None (default) and it will be calculated.
    :type dliq: float or None
    :arg dhum: Humid air density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dhum: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg salt0: Initial guess for the salinity in kg/kg. If None
        (default) then the appropriate `_approx_*` is used.
    :type salt0: float or None
    :arg airf0: Initial guess for the humid air dry fraction in kg/kg.
        If None (default) then the appropriate `_approx_*` is used.
    :type airf0: float or None
    :arg temp0: Initial guess for the temperature in K. If None
        (default) then the appropriate `_approx_*` is used.
    :type temp0: float or None
    :arg pres0: Initial guess for the pressure in Pa. If None (default)
        then the appropriate `_approx_*` is used.
    :type pres0: float or None
    :arg dliq0: Initial guess for the liquid water density in kg/m3. If
        None (default) then `flu3a._dliq_default` is used.
    :type dliq0: float or None
    :arg dhum0: Initial guess for the liquid water density in kg/m3. If
        None (default) then the appropriate `_approx_*` is used.
    :type dhum0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg bool useext: If False (default) then the salt contribution is
        calculated from _GSCOEFFS; if True, from _GSCOEFFS_EXT.
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Salinity, humid air dry fraction, temperature, pressure,
        seawater liquid density, and humid air density (all in SI
        units).
    :raises ValueError: If fewer than three of salt, airf, temp, and
        pres are provided.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    """
    if sum(val is None for val in (salt,airf,temp,pres)) > 1:
        errmsg = 'Must provide at least three of (salt,airf,temp,pres)'
        raise ValueError(errmsg)
    if mathargs is None:
        mathargs = dict()
    fkwargs = {'useext': useext}
    
    if salt is None:
        dhum = air3a.eq_atp(airf,temp,pres,dhum=dhum,dhum0=dhum0,
            mathargs=mathargs)
        dliq = flu3a.eq_tp_liq(temp,pres,dliq=dliq,dliq0=dliq0,
            mathargs=mathargs)
        fargs = (airf,temp,pres,dhum,dliq)
        salt = _newton(_diff_atp,salt0,_approx_atp,fargs=fargs,fkwargs=fkwargs,
            **mathargs)
    elif airf is None:
        dliq = flu3a.eq_tp_liq(temp,pres,dliq=dliq,dliq0=dliq0,
            mathargs=mathargs)
        x0 = (airf0,dhum0)
        fargs = (salt,temp,pres,dliq)
        x1 = _newton(_diff_stp,x0,_approx_stp,fargs=fargs,fkwargs=fkwargs,
            **mathargs)
        airf, dhum = x1
    elif temp is None:
        x0 = (temp0,dhum0,dliq0)
        fargs = (salt,airf,pres)
        x1 = _newton(_diff_sap,x0,_approx_sap,fargs=fargs,fkwargs=fkwargs,
            **mathargs)
        temp, dhum, dliq = x1
    elif pres is None:
        x0 = (pres0,dhum0,dliq0)
        fargs = (salt,airf,temp)
        x1 = _newton(_diff_sat,x0,_approx_sat,fargs=fargs,fkwargs=fkwargs,
            **mathargs)
        pres, dhum, dliq = x1
    if dhum is None:
        dhum = air3a.eq_atp(airf,temp,pres,dhum0=dhum0,mathargs=mathargs)
    if dliq is None:
        dliq = flu3a.eq_tp_liq(temp,pres,dliq0=dliq0,mathargs=mathargs)
    
    _chkflubnds(temp,dliq,chkbnd=chkbnd)
    _chkhumbnds(airf,temp,dhum,chkbnd=chkbnd)
    _chksalbnds(salt,temp,pres,chkbnd=chkbnd)
    if not chkvals:
        return salt, airf, temp, pres, dhum, dliq
    
    lhs, rhs, __, __ = _diff_stp(airf,dhum,dliq,salt,temp,pres,useext=useext)
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
    return salt, airf, temp, pres, dhum, dliq


## Thermodynamic properties
def densityair(salt=None,airf=None,temp=None,pres=None,dhum=None,
    dliq=None,chkvals=False,chktol=_CHKTOL,salt0=None,airf0=None,
    temp0=None,pres0=None,dhum0=None,dliq0=None,chkbnd=False,
    useext=False,mathargs=None):
    """Calculate sea-air humid air density.
    
    Calculate the density of humid air in sea-air.
    
    :arg salt: Salinity in kg/kg. If unknown, pass None (default) and it
        will be calculated.
    :type salt: float or None
    :arg airf: Humid air dry fraction in kg/kg. If unknown, pass None
        (default) and it will be calculated.
    :type airf: float or None
    :arg temp: Temperature in K. If unknown, pass None (default) and it
        will be calculated.
    :type temp: float or None
    :arg pres: Pressure in Pa. If unknown, pass None (default) and it
        will be calculated.
    :type pres: float or None
    :arg dliq: Seawater liquid water density in kg/m3. If unknown, pass
        None (default) and it will be calculated.
    :type dliq: float or None
    :arg dhum: Humid air density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dhum: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg salt0: Initial guess for the salinity in kg/kg. If None
        (default) then the appropriate `_approx_*` is used.
    :type salt0: float or None
    :arg airf0: Initial guess for the humid air dry fraction in kg/kg.
        If None (default) then the appropriate `_approx_*` is used.
    :type airf0: float or None
    :arg temp0: Initial guess for the temperature in K. If None
        (default) then the appropriate `_approx_*` is used.
    :type temp0: float or None
    :arg pres0: Initial guess for the pressure in Pa. If None (default)
        then the appropriate `_approx_*` is used.
    :type pres0: float or None
    :arg dliq0: Initial guess for the liquid water density in kg/m3. If
        None (default) then `flu3a._dliq_default` is used.
    :type dliq0: float or None
    :arg dhum0: Initial guess for the liquid water density in kg/m3. If
        None (default) then the appropriate `_approx_*` is used.
    :type dhum0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg bool useext: If False (default) then the salt contribution is
        calculated from _GSCOEFFS; if True, from _GSCOEFFS_EXT.
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Humid air density in kg/m3.
    :raises ValueError: If fewer than three of salt, airf, temp, and
        pres are provided.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> densityair(salt=0.035,airf=0.99,pres=1e5)
    1.20553990576
    >>> densityair(salt=0.035,temp=300.,pres=1e5)
    1.14642944448
    """
    salt, airf, temp, pres, dhum, dliq = eq_satp(salt=salt,airf=airf,temp=temp,
        pres=pres,dhum=dhum,dliq=dliq,chkvals=chkvals,chktol=chktol,salt0=salt0,
        airf0=airf0,temp0=temp0,pres0=pres0,dhum0=dhum0,dliq0=dliq0,
        chkbnd=chkbnd,useext=useext,mathargs=mathargs)
    return dhum

def densityvap(salt=None,airf=None,temp=None,pres=None,dhum=None,
    dliq=None,chkvals=False,chktol=_CHKTOL,salt0=None,airf0=None,
    temp0=None,pres0=None,dhum0=None,dliq0=None,chkbnd=False,
    useext=False,mathargs=None):
    """Calculate sea-air water vapour density.
    
    Calculate the partial density of water vapour in sea-air.
    
    :arg salt: Salinity in kg/kg. If unknown, pass None (default) and it
        will be calculated.
    :type salt: float or None
    :arg airf: Humid air dry fraction in kg/kg. If unknown, pass None
        (default) and it will be calculated.
    :type airf: float or None
    :arg temp: Temperature in K. If unknown, pass None (default) and it
        will be calculated.
    :type temp: float or None
    :arg pres: Pressure in Pa. If unknown, pass None (default) and it
        will be calculated.
    :type pres: float or None
    :arg dliq: Seawater liquid water density in kg/m3. If unknown, pass
        None (default) and it will be calculated.
    :type dliq: float or None
    :arg dhum: Humid air density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dhum: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg salt0: Initial guess for the salinity in kg/kg. If None
        (default) then the appropriate `_approx_*` is used.
    :type salt0: float or None
    :arg airf0: Initial guess for the humid air dry fraction in kg/kg.
        If None (default) then the appropriate `_approx_*` is used.
    :type airf0: float or None
    :arg temp0: Initial guess for the temperature in K. If None
        (default) then the appropriate `_approx_*` is used.
    :type temp0: float or None
    :arg pres0: Initial guess for the pressure in Pa. If None (default)
        then the appropriate `_approx_*` is used.
    :type pres0: float or None
    :arg dliq0: Initial guess for the liquid water density in kg/m3. If
        None (default) then `flu3a._dliq_default` is used.
    :type dliq0: float or None
    :arg dhum0: Initial guess for the liquid water density in kg/m3. If
        None (default) then the appropriate `_approx_*` is used.
    :type dhum0: float or None
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
    :raises ValueError: If fewer than three of salt, airf, temp, and
        pres are provided.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> densityvap(salt=0.035,airf=0.99,pres=1e5)
    1.20553990576e-02
    >>> densityvap(salt=0.035,temp=300.,pres=1e5)
    2.51876465812e-02
    """
    salt, airf, temp, pres, dhum, dliq = eq_satp(salt=salt,airf=airf,temp=temp,
        pres=pres,dhum=dhum,dliq=dliq,chkvals=chkvals,chktol=chktol,salt0=salt0,
        airf0=airf0,temp0=temp0,pres0=pres0,dhum0=dhum0,dliq0=dliq0,
        chkbnd=chkbnd,useext=useext,mathargs=mathargs)
    dvap = dhum * (1 - airf)
    return dvap

def entropyair(salt=None,airf=None,temp=None,pres=None,dhum=None,
    dliq=None,chkvals=False,chktol=_CHKTOL,salt0=None,airf0=None,
    temp0=None,pres0=None,dhum0=None,dliq0=None,chkbnd=False,
    useext=False,mathargs=None):
    """Calculate sea-air humid air entropy.
    
    Calculate the specific entropy of humid air in sea-air.
    
    :arg salt: Salinity in kg/kg. If unknown, pass None (default) and it
        will be calculated.
    :type salt: float or None
    :arg airf: Humid air dry fraction in kg/kg. If unknown, pass None
        (default) and it will be calculated.
    :type airf: float or None
    :arg temp: Temperature in K. If unknown, pass None (default) and it
        will be calculated.
    :type temp: float or None
    :arg pres: Pressure in Pa. If unknown, pass None (default) and it
        will be calculated.
    :type pres: float or None
    :arg dliq: Seawater liquid water density in kg/m3. If unknown, pass
        None (default) and it will be calculated.
    :type dliq: float or None
    :arg dhum: Humid air density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dhum: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg salt0: Initial guess for the salinity in kg/kg. If None
        (default) then the appropriate `_approx_*` is used.
    :type salt0: float or None
    :arg airf0: Initial guess for the humid air dry fraction in kg/kg.
        If None (default) then the appropriate `_approx_*` is used.
    :type airf0: float or None
    :arg temp0: Initial guess for the temperature in K. If None
        (default) then the appropriate `_approx_*` is used.
    :type temp0: float or None
    :arg pres0: Initial guess for the pressure in Pa. If None (default)
        then the appropriate `_approx_*` is used.
    :type pres0: float or None
    :arg dliq0: Initial guess for the liquid water density in kg/m3. If
        None (default) then `flu3a._dliq_default` is used.
    :type dliq0: float or None
    :arg dhum0: Initial guess for the liquid water density in kg/m3. If
        None (default) then the appropriate `_approx_*` is used.
    :type dhum0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg bool useext: If False (default) then the salt contribution is
        calculated from _GSCOEFFS; if True, from _GSCOEFFS_EXT.
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Humid air entropy in J/kg/K.
    :raises ValueError: If fewer than three of salt, airf, temp, and
        pres are provided.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> entropyair(salt=0.035,temp=300.,pres=1e5)
    293.150672957
    """
    salt, airf, temp, pres, dhum, dliq = eq_satp(salt=salt,airf=airf,temp=temp,
        pres=pres,dhum=dhum,dliq=dliq,chkvals=chkvals,chktol=chktol,salt0=salt0,
        airf0=airf0,temp0=temp0,pres0=pres0,dhum0=dhum0,dliq0=dliq0,
        chkbnd=chkbnd,useext=useext,mathargs=mathargs)
    shum = -_air_f(0,1,0,airf,temp,dhum)
    return shum

def enthalpyevap(salt=None,airf=None,temp=None,pres=None,dhum=None,
    dliq=None,chkvals=False,chktol=_CHKTOL,salt0=None,airf0=None,
    temp0=None,pres0=None,dhum0=None,dliq0=None,chkbnd=False,
    useext=False,mathargs=None):
    """Calculate sea-air enthalpy of evaporation.
    
    Calculate the specific enthalpy of evaporation for sea-air.
    
    :arg salt: Salinity in kg/kg. If unknown, pass None (default) and it
        will be calculated.
    :type salt: float or None
    :arg airf: Humid air dry fraction in kg/kg. If unknown, pass None
        (default) and it will be calculated.
    :type airf: float or None
    :arg temp: Temperature in K. If unknown, pass None (default) and it
        will be calculated.
    :type temp: float or None
    :arg pres: Pressure in Pa. If unknown, pass None (default) and it
        will be calculated.
    :type pres: float or None
    :arg dliq: Seawater liquid water density in kg/m3. If unknown, pass
        None (default) and it will be calculated.
    :type dliq: float or None
    :arg dhum: Humid air density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dhum: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg salt0: Initial guess for the salinity in kg/kg. If None
        (default) then the appropriate `_approx_*` is used.
    :type salt0: float or None
    :arg airf0: Initial guess for the humid air dry fraction in kg/kg.
        If None (default) then the appropriate `_approx_*` is used.
    :type airf0: float or None
    :arg temp0: Initial guess for the temperature in K. If None
        (default) then the appropriate `_approx_*` is used.
    :type temp0: float or None
    :arg pres0: Initial guess for the pressure in Pa. If None (default)
        then the appropriate `_approx_*` is used.
    :type pres0: float or None
    :arg dliq0: Initial guess for the liquid water density in kg/m3. If
        None (default) then `flu3a._dliq_default` is used.
    :type dliq0: float or None
    :arg dhum0: Initial guess for the liquid water density in kg/m3. If
        None (default) then the appropriate `_approx_*` is used.
    :type dhum0: float or None
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
    :raises ValueError: If fewer than three of salt, airf, temp, and
        pres are provided.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> enthalpyevap(salt=0.035,airf=0.99,pres=1e5)
    2464765.77588
    >>> enthalpyevap(salt=0.035,temp=300.,pres=1e5)
    2434570.55229
    """
    salt, airf, temp, pres, dhum, dliq = eq_satp(salt=salt,airf=airf,temp=temp,
        pres=pres,dhum=dhum,dliq=dliq,chkvals=chkvals,chktol=chktol,salt0=salt0,
        airf0=airf0,temp0=temp0,pres0=pres0,dhum0=dhum0,dliq0=dliq0,
        chkbnd=chkbnd,useext=useext,mathargs=mathargs)
    
    # Calculate water vapour enthalpy
    fh = _air_f(0,0,0,airf,temp,dhum)
    fh_a = _air_f(1,0,0,airf,temp,dhum)
    fh_t = _air_f(0,1,0,airf,temp,dhum)
    fh_d = _air_f(0,0,1,airf,temp,dhum)
    fh_at = _air_f(1,1,0,airf,temp,dhum)
    fh_ad = _air_f(1,0,1,airf,temp,dhum)
    fh_td = _air_f(0,1,1,airf,temp,dhum)
    fh_dd = _air_f(0,0,2,airf,temp,dhum)
    ch = 2*fh_d + dhum*fh_dd
    hh = fh - temp*fh_t + dhum*fh_d
    hh_a = fh_a - temp*fh_at + temp*dhum/ch*fh_ad*fh_td
    hv = hh - airf*hh_a
    hl = flu2.enthalpy(temp,dliq)
    
    # Treat pure water case separately to avoid numerical errors
    if salt == 0:
        hevap = hv - hl
        return hevap
    
    gs = _sal_g(0,0,0,salt,temp,pres,useext=useext)
    gs_s = _sal_g(1,0,0,salt,temp,pres,useext=useext)
    gs_t = _sal_g(0,1,0,salt,temp,pres,useext=useext)
    gs_st = _sal_g(1,1,0,salt,temp,pres,useext=useext)
    hls = gs - temp*gs_t - salt*(gs_s - temp*gs_st)
    hevap = hv - hl - hls
    return hevap


## Thermodynamic functions of 3 variables
def massfractionair(salt,temp,pres,airf=None,dhum=None,dliq=None,
    chkvals=False,chktol=_CHKTOL,airf0=None,dhum0=None,dliq0=None,
    chkbnd=False,useext=False,mathargs=None):
    """Calculate sea-air dry fraction.
    
    Calculate the humid air dry fraction for sea-air.

    :arg float airf: Humid air dry fraction in kg/kg.
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg salt: Salinity in kg/kg. If unknown, pass None (default) and it
        will be calculated.
    :type salt: float or None
    :arg dliq: Seawater liquid water density in kg/m3. If unknown, pass
        None (default) and it will be calculated.
    :type dliq: float or None
    :arg dhum: Humid air density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dhum: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg airf0: Initial guess for the humid air dry fraction in kg/kg.
        If None (default) then the appropriate `_approx_*` is used.
    :type airf0: float or None
    :arg dliq0: Initial guess for the liquid water density in kg/m3. If
        None (default) then `flu3a._dliq_default` is used.
    :type dliq0: float or None
    :arg dhum0: Initial guess for the liquid water density in kg/m3. If
        None (default) then the appropriate `_approx_*` is used.
    :type dhum0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg bool useext: If False (default) then the salt contribution is
        calculated from _GSCOEFFS; if True, from _GSCOEFFS_EXT.
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Humid air dry fraction in kg/kg.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> massfractionair(0.035,300.,1e5)
    0.978029483888
    """
    __, airf, __, __, dhum, dliq = eq_satp(salt=salt,temp=temp,pres=pres,
        airf=airf,dhum=dhum,dliq=dliq,chkvals=chkvals,chktol=chktol,airf0=airf0,
        dhum0=dhum0,dliq0=dliq0,chkbnd=chkbnd,useext=useext,mathargs=mathargs)
    return airf

def vapourpressure(salt,temp,pres,airf=None,dhum=None,dliq=None,
    chkvals=False,chktol=_CHKTOL,airf0=None,dhum0=None,dliq0=None,
    chkbnd=False,useext=False,mathargs=None):
    """Calculate sea-air vapour pressure.
    
    Calculate the partial pressure of water vapour in sea-air.
    
    :arg float salt: Salinity in kg/kg.
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg airf: Humid air dry fraction in kg/kg. If unknown, pass None
        (default) and it will be calculated.
    :type airf: float or None
    :arg dliq: Seawater liquid water density in kg/m3. If unknown, pass
        None (default) and it will be calculated.
    :type dliq: float or None
    :arg dhum: Humid air density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dhum: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg airf0: Initial guess for the humid air dry fraction in kg/kg.
        If None (default) then the appropriate `_approx_*` is used.
    :type airf0: float or None
    :arg dliq0: Initial guess for the liquid water density in kg/m3. If
        None (default) then `flu3a._dliq_default` is used.
    :type dliq0: float or None
    :arg dhum0: Initial guess for the liquid water density in kg/m3. If
        None (default) then the appropriate `_approx_*` is used.
    :type dhum0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg bool useext: If False (default) then the salt contribution is
        calculated from _GSCOEFFS; if True, from _GSCOEFFS_EXT.
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Vapour pressure in Pa.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> vapourpressure(0.035,300.,1e5)
    3485.92986681
    """
    __, airf, __, __, dhum, dliq = eq_satp(salt=salt,temp=temp,pres=pres,
        airf=airf,dhum=dhum,dliq=dliq,chkvals=chkvals,chktol=chktol,airf0=airf0,
        dhum0=dhum0,dliq0=dliq0,chkbnd=chkbnd,useext=useext,mathargs=mathargs)
    xv = convert0.air_molfractionvap(airf)
    pvap = xv * pres
    return pvap

def condensetemp(salt,airf,pres,temp=None,dhum=None,dliq=None,
    chkvals=False,chktol=_CHKTOL,temp0=None,dhum0=None,dliq0=None,
    chkbnd=False,useext=False,mathargs=None):
    """Calculate sea-air condensation temperature.
    
    Calculate the condensation temperature of sea-air.
    
    :arg salt: Salinity in kg/kg. If unknown, pass None (default) and it
        will be calculated.
    :type salt: float or None
    :arg airf: Humid air dry fraction in kg/kg. If unknown, pass None
        (default) and it will be calculated.
    :type airf: float or None
    :arg temp: Temperature in K. If unknown, pass None (default) and it
        will be calculated.
    :type temp: float or None
    :arg pres: Pressure in Pa. If unknown, pass None (default) and it
        will be calculated.
    :type pres: float or None
    :arg dliq: Seawater liquid water density in kg/m3. If unknown, pass
        None (default) and it will be calculated.
    :type dliq: float or None
    :arg dhum: Humid air density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dhum: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg salt0: Initial guess for the salinity in kg/kg. If None
        (default) then the appropriate `_approx_*` is used.
    :type salt0: float or None
    :arg airf0: Initial guess for the humid air dry fraction in kg/kg.
        If None (default) then the appropriate `_approx_*` is used.
    :type airf0: float or None
    :arg temp0: Initial guess for the temperature in K. If None
        (default) then the appropriate `_approx_*` is used.
    :type temp0: float or None
    :arg pres0: Initial guess for the pressure in Pa. If None (default)
        then the appropriate `_approx_*` is used.
    :type pres0: float or None
    :arg dliq0: Initial guess for the liquid water density in kg/m3. If
        None (default) then `flu3a._dliq_default` is used.
    :type dliq0: float or None
    :arg dhum0: Initial guess for the liquid water density in kg/m3. If
        None (default) then the appropriate `_approx_*` is used.
    :type dhum0: float or None
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
    
    >>> condensetemp(0.035,0.99,1e5)
    287.367451766
    """
    __, __, temp, __, dhum, dliq = eq_satp(salt=salt,airf=airf,pres=pres,
        temp=temp,dhum=dhum,dliq=dliq,chkvals=chkvals,chktol=chktol,temp0=temp0,
        dhum0=dhum0,dliq0=dliq0,chkbnd=chkbnd,useext=useext,mathargs=mathargs)
    return temp


## Chemical potential difference
def chempotevap(airf,salt,temp,pres,dhum=None,dliq=None,chkvals=False,
    chktol=_CHKTOL,dhum0=None,dliq0=None,chkbnd=False,useext=False,
    mathargs=None):
    """Calculate seawater-air chemical potential difference.
    
    Calculate the chemical potential difference between liquid water in
    seawater (at the given salinity, temperature, and pressure) and
    water vapour in humid air (at the given dry air mass fraction,
    temperature, and pressure), scaled by the temperature. This is a
    measure of the sea-air disequilibrium, and the Onsager force in the
    entropy equation.

    :arg airf: Humid air dry fraction in kg/kg.
    :arg salt: Salinity in kg/kg.
    :arg temp: Temperature in K.
    :arg pres: Pressure in Pa.
    :arg dliq: Seawater liquid water density in kg/m3. If unknown, pass
        None (default) and it will be calculated.
    :type dliq: float or None
    :arg dhum: Humid air density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dhum: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg dliq0: Initial guess for the liquid water density in kg/m3. If
        None (default) then `flu3a._dliq_default` is used.
    :type dliq0: float or None
    :arg dhum0: Initial guess for the liquid water density in kg/m3. If
        None (default) then the appropriate `_approx_*` is used.
    :type dhum0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg bool useext: If False (default) then the salt contribution is
        calculated from _GSCOEFFS; if True, from _GSCOEFFS_EXT.
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Scaled potential difference (unitless).
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> chempotevap(0.9,0.035,300.,1e5)
    1.45584069071
    """
    muv = air3b.vappot(airf,temp,pres,dhum=dhum,chkvals=chkvals,chktol=chktol,
        dhum0=dhum0,chkbnd=chkbnd,mathargs=mathargs)
    mul = sea3a.liqpot(salt,temp,pres,dliq=dliq,chkvals=chkvals,chktol=chktol,
        dliq0=dliq0,chkbnd=chkbnd,useext=useext,mathargs=mathargs)
    dphi = (muv - mul) / (_RWAT*temp)
    return dphi

