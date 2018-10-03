"""Liquid water-humid air equilibrium functions.

This module provides functions to get the values of primary variables
for liquid water saturated (wet) air, the equilibrium of dry air, water
vapour, and liquid water.

:Examples:

>>> pressure(airf=.99,temp=300.)
223057.741750
>>> entropy(airf=.99,temp=300.)
-41.9991507402
>>> enthalpyevap(airf=.99,temp=300.)
2433111.29696
>>> temperature(airf=.99,pres=1e5)
287.078299795
>>> entropy(airf=.99,pres=1e5)
145.863545194
>>> enthalpyevap(airf=.99,pres=1e5)
2465656.38630
>>> enthalpyevap(temp=300.,pres=1e5)
2434585.53919
>>> temperature(airf=.99,entr=100.)
290.107386673
>>> pressure(airf=.99,entr=100.)
121546.373652
>>> enthalpyevap(airf=.99,entr=100.)
2458121.74961
>>> icl(0.99,300.,1e5)
82723.6047631
>>> airffromrh_wmo(0.8,300.,1e5)
0.982004037135
>>> rhfromairf_wmo(0.99,300.,1e5)
0.440953686019

:Functions:

* :func:`eq_atpe`: Calculate wet air equilibrium properties at any two
  of dry fraction, temperature, and pressure, or at dry fraction and
  entropy.
* :func:`temperature`: Temperature at liquid-air equilibrium.
* :func:`pressure`: Pressure at liquid-air equilibrium.
* :func:`densityair`: Humid air density at liquid-air equilibrium.
* :func:`densityliq`: Liquid water density at liquid-air equilibrium.
* :func:`entropy`: Entropy at liquid-air equilibrium.
* :func:`densityvap`: Water vapour density at liquid-air equilibrium.
* :func:`enthalpyevap`: Enthalpy of evaporation.
* :func:`condensationpressure`: Pressure at which liquid water will
  condense out of humid air.
* :func:`dewpoint`: Temperature at which liquid water will condense out
  of humid air.
* :func:`massfractionair`: Dry air mass fraction at liquid-air
  equilibrium.
* :func:`eq_icl`: Calculate wet air equilibrium properties at the
  isentropic condensation level.
* :func:`icl`: Isentropic condensation pressure.
* :func:`ict`: Isentropic condensation temperature.
* :func:`airffromrh_wmo`: Calculate dry air fraction from relative
  humidity (WMO definition).
* :func:`rhfromairf_wmo`: Calculate relative humidity (WMO definition)
  from dry air fraction.
* :func:`airffromrh_cct`: Calculate dry air fraction from relative
  humidity (IUPAC/CCT definition).
* :func:`rhfromairf_cct`: Calculate relative humidity (IUPAC/CCT
  definition) from dry air fraction.

"""

__all__ = ['eq_atpe',
    'temperature','pressure','densityair','densityliq','entropy','densityvap',
    'enthalpyevap',
    'condensationpressure','dewpoint','massfractionair',
    'eq_icl','icl','ict',
    'airffromrh_wmo','rhfromairf_wmo','airffromrh_cct','rhfromairf_cct']

import warnings
import numpy
import constants0
import convert0
import air2
import flu2
import maths3
import air3a
import flu3a
import maths4

_CHKTOL = constants0.CHKTOL
_RWAT = constants0.RWAT
_RDRY = constants0.RDRY
_PATM = constants0.PATM
_TCELS = constants0.TCELS
_TTP = constants0.TTP
_PTPE = constants0.PTPE
_LLVTP = constants0.LLVTP
_CDRY = constants0.CDRY
_CVAP = constants0.CVAP
_CLIQ = constants0.CLIQ
_chkhumbnds = constants0.chkhumbnds
_chkflubnds = constants0.chkflubnds
_air_f = air2.air_f
_air_eq_pressure = air2.eq_pressure
_air_eq_vappot = air2.eq_vappot
_air_eq_entropy = air2.eq_entropy
_flu_eq_chempot = flu2.eq_chempot
_flu_eq_pressure = flu2.eq_pressure
_newton = maths3.newton
_dliq_default = flu3a._dliq_default
_EPSW = _RDRY/_RWAT
_AVL = _LLVTP / (_RWAT*_TTP)
_BVL = (_CLIQ-_CVAP) / _RWAT


## Equilibrium functions
def _approx_tp(temp,pres):
    """Approximate ADhDl at TP.
    
    Approximate the humid air dry fraction, humid air density, and
    liquid water density of wet air at the given temperature and
    pressure.
    
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :returns: Humid air dry fraction, humid air density, and liquid
        water density (all in SI units).
    """
    earg = _AVL*(1-_TTP/temp) + _BVL*(1-_TTP/temp - numpy.log(temp/_TTP))
    pvap = _PTPE * numpy.exp(earg)
    airf = (pres-pvap)/(pres-pvap + _EPSW*pvap)
    dhum = (pres-pvap + _EPSW*pvap)/(_RDRY*temp)
    dliq = _dliq_default(temp,pres)
    return airf, dhum, dliq

def _approx_at(airf,temp):
    """Approximate DhDl at AT.
    
    Approximate humid air density and liquid water density of wet air at
    the given dry fraction and temperature.
    
    :arg float airf: Dry air mass fraction in kg/kg.
    :arg float temp: Temperature in K.
    :returns: Humid air density and liquid water density (all in SI
        units).
    """
    earg = _AVL*(1-_TTP/temp) + _BVL*(1-_TTP/temp - numpy.log(temp/_TTP))
    pvap = _PTPE * numpy.exp(earg)
    pres = (pvap*(1-airf) + airf*_EPSW*pvap)/(1-airf)
    dhum = (pres-pvap + _EPSW*pvap)/(_RDRY*temp)
    dliq = _dliq_default(temp,pres)
    return dhum, dliq

def _approx_ap(airf,pres):
    """Approximate TDhDl at AP.
    
    Approximate the temperature, humid air density, and liquid water
    density of wet air at the given dry fraction and pressure.
    
    :arg float airf: Dry air mass fraction in kg/kg.
    :arg float pres: Pressure in Pa.
    :returns: Temperature, humid air density, and liquid water density
        (all in SI units).
    """
    pvap = pres * (1-airf)/(airf*_EPSW + 1-airf)
    v = numpy.log(pvap/_PTPE)/_BVL
    x = maths4.lamb2(v,_AVL/_BVL)
    temp = _TTP/x
    dhum = (pres-pvap + _EPSW*pvap) / (_RDRY*temp)
    dliq = _dliq_default(temp,pres)
    return temp, dhum, dliq

def _approx_ae(airf,entr):
    """Approximate TDhDl at AE.
    
    Approximate the temperature, humid air density, and liquid water
    density of wet air at the given dry fraction and specific entropy.
    
    :arg float airf: Dry air mass fraction in kg/kg.
    :arg float entr: Entropy in J/kg/K.
    :returns: Temperature, humid air density, and liquid water density
        (all in SI units).
    """
    ceff = airf*_CDRY + (1-airf)*_CVAP
    reff = airf*_RDRY + (1-airf)*_RWAT
    gam = ceff/reff
    s0 = airf*_CDRY*numpy.log(_TTP/_TCELS) - airf*_RDRY*numpy.log(_PTPE/_PATM)
    s0 += (1-airf)*_LLVTP/_TTP - airf*_RDRY*numpy.log(_EPSW*airf/(1-airf))
    r = (_AVL-gam)/(gam+_BVL)
    v = (s0 - entr)/(gam+_BVL)/reff
    x = maths4.lamb2(v,r)
    temp = _TTP/x
    dhum, dliq = _approx_at(airf,temp)
    return temp, dhum, dliq

def _diff_tp(a,dh,dl,temp,pres):
    """Calculate wet air disequilibrium at TP.
    
    Calculate both sides of the equations
    
        given pressure = pressure in humid air
        given pressure = pressure in liquid water
        chemical potential of liquid water = potential of water vapour
    
    and their Jacobians with respect to dry fraction, humid air density,
    and liquid water density. Solving these equations gives the
    equilibrium values at the given temperature and pressure.
    
    :arg float a: Dry air mass fraction in kg/kg.
    :arg float dh: Humid air density in kg/m3.
    :arg float dl: Liquid water density in kg/m3.
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :returns: Left-hand side of the equation, right-hand side,
        Jacobian of LHS, and Jacobian of RHS.
    :rtype: tuple(array(float))
    """
    ph = _air_eq_pressure(0,0,0,a,temp,dh)
    pl = _flu_eq_pressure(0,0,temp,dl)
    gv = _air_eq_vappot(0,0,0,a,temp,dh)
    gl = _flu_eq_chempot(0,0,temp,dl)
    lhs = numpy.array([pres, pres, gl])
    rhs = numpy.array([ph, pl, gv])
    
    ph_a = _air_eq_pressure(1,0,0,a,temp,dh)
    ph_d = _air_eq_pressure(0,0,1,a,temp,dh)
    pl_d = _flu_eq_pressure(0,1,temp,dl)
    gv_a = _air_eq_vappot(1,0,0,a,temp,dh)
    gv_d = _air_eq_vappot(0,0,1,a,temp,dh)
    gl_d = _flu_eq_chempot(0,1,temp,dl)
    dlhs = numpy.array([[0.,0.,0.], [0.,0.,0.], [0.,0.,gl_d]])
    drhs = numpy.array([[ph_a,ph_d,0.], [0.,0.,pl_d], [gv_a,gv_d,0.]])
    return lhs, rhs, dlhs, drhs

def _diff_at(dh,dl,airf,temp):
    """Calculate wet air disequilibrium at AT.
    
    Calculate both sides of the equations
    
        pressure in liquid water = pressure in humid air
        chemical potential of liquid water = potential of water vapour
    
    and their Jacobians with respect to humid air density and liquid
    water density. Solving these equations gives equilibrium values at
    the given dry air mass fraction and temperature.
    
    :arg float dh: Humid air density in kg/m3.
    :arg float dl: Liquid water density in kg/m3.
    :arg float airf: Dry air mass fraction in kg/kg.
    :arg float temp: Temperature in K.
    :returns: Left-hand side of the equation, right-hand side,
        Jacobian of LHS, and Jacobian of RHS.
    :rtype: tuple(array(float))
    """
    ph = _air_eq_pressure(0,0,0,airf,temp,dh)
    pl = _flu_eq_pressure(0,0,temp,dl)
    gv = _air_eq_vappot(0,0,0,airf,temp,dh)
    gl = _flu_eq_chempot(0,0,temp,dl)
    lhs = numpy.array([pl, gl])
    rhs = numpy.array([ph, gv])
    
    ph_d = _air_eq_pressure(0,0,1,airf,temp,dh)
    pl_d = _flu_eq_pressure(0,1,temp,dl)
    gv_d = _air_eq_vappot(0,0,1,airf,temp,dh)
    gl_d = _flu_eq_chempot(0,1,temp,dl)
    dlhs = numpy.array([[0.,pl_d], [0.,gl_d]])
    drhs = numpy.array([[ph_d,0.], [gv_d,0.]])
    return lhs, rhs, dlhs, drhs

def _diff_ap(t,dh,dl,airf,pres):
    """Calculate wet air disequilibrium at AP.
    
    Calculate both sides of the equations
    
        given pressure = pressure in humid air
        given pressure = pressure in liquid water
        chemical potential of liquid water = potential of water vapour
    
    and their Jacobians with respect to temperature, humid air density,
    and liquid water density. Solving these equations gives equilibrium
    values at the given dry air mass fraction and pressure.
    
    :arg float t: Temperature in K.
    :arg float dh: Humid air density in kg/m3.
    :arg float dl: Liquid water density in kg/m3.
    :arg float airf: Dry air mass fraction in kg/kg.
    :arg float pres: Pressure in Pa.
    :returns: Left-hand side of the equation, right-hand side,
        Jacobian of LHS, and Jacobian of RHS.
    :rtype: tuple(array(float))
    """
    ph = _air_eq_pressure(0,0,0,airf,t,dh)
    pl = _flu_eq_pressure(0,0,t,dl)
    gv = _air_eq_vappot(0,0,0,airf,t,dh)
    gl = _flu_eq_chempot(0,0,t,dl)
    lhs = numpy.array([pres, pres, gl])
    rhs = numpy.array([ph, pl, gv])
    
    ph_a = _air_eq_pressure(1,0,0,airf,t,dh)
    ph_t = _air_eq_pressure(0,1,0,airf,t,dh)
    ph_d = _air_eq_pressure(0,0,1,airf,t,dh)
    pl_t = _flu_eq_pressure(1,0,t,dl)
    pl_d = _flu_eq_pressure(0,1,t,dl)
    gv_a = _air_eq_vappot(1,0,0,airf,t,dh)
    gv_t = _air_eq_vappot(0,1,0,airf,t,dh)
    gv_d = _air_eq_vappot(0,0,1,airf,t,dh)
    gl_t = _flu_eq_chempot(1,0,t,dl)
    gl_d = _flu_eq_chempot(0,1,t,dl)
    dlhs = numpy.array([[0.,0.,0.], [0.,0.,0.], [gl_t,0.,gl_d]])
    drhs = numpy.array([[ph_t,ph_d,0.], [pl_t,0.,pl_d], [gv_t,gv_d,0.]])
    return lhs, rhs, dlhs, drhs

def _diff_ae(t,dh,dl,airf,entr):
    """Calculate wet air disequilibrium at AE.
    
    Calculate both sides of the equations
    
        pressure in liquid water = pressure in humid air
        chemical potential of liquid water = potential of water vapour
        given entropy = entropy of humid air
    
    and their Jacobians with respect to temperature, humid air density,
    and liquid water density. Solving these equations gives equilibrium
    values at the given dry air mass fraction and specific entropy.
    
    :arg float t: Temperature in K.
    :arg float dh: Humid air density in kg/m3.
    :arg float dl: Liquid water density in kg/m3.
    :arg float airf: Dry air mass fraction in kg/kg.
    :arg float entr: Entropy in J/kg/K.
    :returns: Left-hand side of the equation, right-hand side,
        Jacobian of LHS, and Jacobian of RHS.
    :rtype: tuple(array(float))
    """
    ph = _air_eq_pressure(0,0,0,airf,t,dh)
    pl = _flu_eq_pressure(0,0,t,dl)
    gv = _air_eq_vappot(0,0,0,airf,t,dh)
    gl = _flu_eq_chempot(0,0,t,dl)
    sh = _air_eq_entropy(0,0,0,airf,t,dh)
    lhs = numpy.array([pl, gl, entr])
    rhs = numpy.array([ph, gv, sh])
    
    ph_t = _air_eq_pressure(0,1,0,airf,t,dh)
    ph_d = _air_eq_pressure(0,0,1,airf,t,dh)
    pl_t = _flu_eq_pressure(1,0,t,dl)
    pl_d = _flu_eq_pressure(0,1,t,dl)
    gv_t = _air_eq_vappot(0,1,0,airf,t,dh)
    gv_d = _air_eq_vappot(0,0,1,airf,t,dh)
    gl_t = _flu_eq_chempot(1,0,t,dl)
    gl_d = _flu_eq_chempot(0,1,t,dl)
    sh_t = _air_eq_entropy(0,1,0,airf,t,dh)
    sh_d = _air_eq_entropy(0,0,1,airf,t,dh)
    dlhs = numpy.array([[pl_t,0.,pl_d], [gl_t,0.,gl_d], [0.,0.,0.]])
    drhs = numpy.array([[ph_t,ph_d,0.], [gv_t,gv_d,0.], [sh_t,sh_d,0.]])
    return lhs, rhs, dlhs, drhs

def eq_atpe(airf=None,temp=None,pres=None,entr=None,dhum=None,dliq=None,
    chkvals=False,chktol=_CHKTOL,airf0=None,temp0=None,pres0=None,
    dhum0=None,dliq0=None,chkbnd=False,mathargs=None):
    """Get primary wet air variables at ATPE.
    
    Get the values of all primary variables for wet air at any two of
    the dry fraction, temperature, and pressure, or at dry fraction and
    specific entropy.
    
    If the calculation has already been done, the results can be passed
    to avoid unnecessary repeat calculations. If enough values are
    passed, they will be checked for consistency if chkvals is True.
    
    :arg airf: Dry air mass fraction in kg/kg. If unknown, pass None
        (default) and it will be calculated.
    :type airf: float or None
    :arg temp: Temperature in K. If unknown, pass None (default) and it
        will be calculated.
    :type temp: float or None
    :arg pres: Pressure in Pa. If unknown, pass None (default) and it
        will be calculated.
    :type pres: float or None
    :arg entr: Entropy in J/kg/K.
    :type entr: float or None
    :arg dhum: Humid air density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dhum: float or None
    :arg dliq: Liquid water density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dliq: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg airf0: Initial guess for the dry fraction in kg/kg. If None
        (default) then the appropriate `_approx_*` is used.
    :type airf0: float or None
    :arg temp0: Initial guess for the temperature in K. If None
        (default) then the appropriate `_approx_*` is used.
    :type temp0: float or None
    :arg pres0: Initial guess for the pressure in Pa. If None (default)
        then the appropriate `_approx_*` is used.
    :type pres0: float or None
    :arg dhum0: Initial guess for the humid air density in kg/m3. If
        None (default) then the appropriate `_approx_*` is used.
    :type dhum0: float or None
    :arg dliq0: Initial guess for the liquid water density in kg/m3. If
        None (default) then the appropriate `_approx_*` is used.
    :type dliq0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Dry air mass fraction, temperature, pressure, humid air
        density, and liquid water density (all in SI units).
    :raises ValueError: If not enough values are provided.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    """
    if mathargs is None:
        mathargs = dict()
    if sum(val is not None for val in (airf,temp,pres)) < 2:
        if airf is None and entr is None:
            errmsg = ('Not enough values were provided. Need any two of '
                '(airf,temp,pres) or both (airf,entr)')
            raise ValueError(errmsg)
        x0 = (temp0,dhum0,dliq0)
        fargs = (airf,entr)
        x1 = _newton(_diff_ae,x0,_approx_ae,fargs=fargs,**mathargs)
        temp, dhum, dliq = x1
        pres = air2.pressure(airf,temp,dhum)
    else:
        if all(val is not None for val in (airf,temp,pres)):
            if dhum is None:
                dhum = air3a.eq_atp(airf,temp,pres,dhum0=dhum0,
                    mathargs=mathargs)
            if dliq is None:
                dliq = flu3a.eq_tp_liq(temp,pres,dliq0=dliq0,mathargs=mathargs)
        elif airf is not None and temp is not None:
            x0 = (dhum0,dliq0)
            fargs = (airf,temp)
            x1 = _newton(_diff_at,x0,_approx_at,fargs=fargs,**mathargs)
            dhum, dliq = x1
            pres = air2.pressure(airf,temp,dhum)
        elif airf is not None and pres is not None:
            x0 = (temp0,dhum0,dliq0)
            fargs = (airf,pres)
            x1 = _newton(_diff_ap,x0,_approx_ap,fargs=fargs,**mathargs)
            temp, dhum, dliq = x1
        else:
            x0 = (airf0,dhum0,dliq0)
            fargs = (temp,pres)
            x1 = _newton(_diff_tp,x0,_approx_tp,fargs=fargs,**mathargs)
            airf, dhum, dliq = x1
    
    _chkhumbnds(airf,temp,dhum,chkbnd=chkbnd)
    _chkflubnds(temp,dliq,chkbnd=chkbnd)
    if not chkvals:
        return airf, temp, pres, dhum, dliq
    
    if entr is None:
        entr = air2.entropy(airf,temp,dhum)
    lhs, rhs, __, __ = _diff_ae(temp,dhum,dliq,airf,entr)
    lhs2, rhs2, __, __ = air3a._diff_atp(dhum,airf,temp,pres)
    lhs.append(lhs2)
    rhs.append(rhs2)
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
    return airf, temp, pres, dhum, dliq


## Equilibrium properties
def temperature(airf=None,temp=None,pres=None,entr=None,dhum=None,
    dliq=None,chkvals=False,chktol=_CHKTOL,airf0=None,temp0=None,
    dhum0=None,dliq0=None,chkbnd=False,mathargs=None):
    """Calculate wet air temperature.
    
    Calculate the temperature of wet air.
    
    :arg airf: Dry air mass fraction in kg/kg. If unknown, pass None
        (default) and it will be calculated.
    :type airf: float or None
    :arg temp: Temperature in K. If unknown, pass None (default) and it
        will be calculated.
    :type temp: float or None
    :arg pres: Pressure in Pa. If unknown, pass None (default) and it
        will be calculated.
    :type pres: float or None
    :arg entr: Entropy in J/kg/K.
    :type entr: float or None
    :arg dhum: Humid air density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dhum: float or None
    :arg dliq: Liquid water density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dliq: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg airf0: Initial guess for the dry fraction in kg/kg. If None
        (default) then the appropriate `_approx_*` is used.
    :type airf0: float or None
    :arg temp0: Initial guess for the temperature in K. If None
        (default) then the appropriate `_approx_*` is used.
    :type temp0: float or None
    :arg pres0: Initial guess for the pressure in Pa. If None (default)
        then the appropriate `_approx_*` is used.
    :type pres0: float or None
    :arg dhum0: Initial guess for the humid air density in kg/m3. If
        None (default) then the appropriate `_approx_*` is used.
    :type dhum0: float or None
    :arg dliq0: Initial guess for the liquid water density in kg/m3. If
        None (default) then the appropriate `_approx_*` is used.
    :type dliq0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Temperature in K.
    :raises ValueError: If not enough values are provided.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> temperature(airf=.99,pres=1e5)
    287.078299795
    >>> temperature(airf=.99,entr=100.)
    290.107386673
    """
    airf, temp, pres, dhum, dliq = eq_atpe(airf=airf,temp=temp,pres=pres,
        entr=entr,dhum=dhum,dliq=dliq,chkvals=chkvals,chktol=chktol,airf0=airf0,
        temp0=temp0,dhum0=dhum0,dliq0=dliq0,chkbnd=chkbnd,mathargs=mathargs)
    return temp

def pressure(airf=None,temp=None,pres=None,entr=None,dhum=None,
    dliq=None,chkvals=False,chktol=_CHKTOL,airf0=None,temp0=None,
    dhum0=None,dliq0=None,chkbnd=False,mathargs=None):
    """Calculate wet air pressure.
    
    Calculate the pressure of wet air.
    
    :arg airf: Dry air mass fraction in kg/kg. If unknown, pass None
        (default) and it will be calculated.
    :type airf: float or None
    :arg temp: Temperature in K. If unknown, pass None (default) and it
        will be calculated.
    :type temp: float or None
    :arg pres: Pressure in Pa. If unknown, pass None (default) and it
        will be calculated.
    :type pres: float or None
    :arg entr: Entropy in J/kg/K.
    :type entr: float or None
    :arg dhum: Humid air density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dhum: float or None
    :arg dliq: Liquid water density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dliq: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg airf0: Initial guess for the dry fraction in kg/kg. If None
        (default) then the appropriate `_approx_*` is used.
    :type airf0: float or None
    :arg temp0: Initial guess for the temperature in K. If None
        (default) then the appropriate `_approx_*` is used.
    :type temp0: float or None
    :arg pres0: Initial guess for the pressure in Pa. If None (default)
        then the appropriate `_approx_*` is used.
    :type pres0: float or None
    :arg dhum0: Initial guess for the humid air density in kg/m3. If
        None (default) then the appropriate `_approx_*` is used.
    :type dhum0: float or None
    :arg dliq0: Initial guess for the liquid water density in kg/m3. If
        None (default) then the appropriate `_approx_*` is used.
    :type dliq0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Pressure in Pa.
    :raises ValueError: If not enough values are provided.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> pressure(airf=.99,entr=100.)
    121546.373652
    >>> pressure(airf=.99,temp=300.)
    223057.741750
    """
    airf, temp, pres, dhum, dliq = eq_atpe(airf=airf,temp=temp,pres=pres,
        entr=entr,dhum=dhum,dliq=dliq,chkvals=chkvals,chktol=chktol,airf0=airf0,
        temp0=temp0,dhum0=dhum0,dliq0=dliq0,chkbnd=chkbnd,mathargs=mathargs)
    return pres

def densityair(airf=None,temp=None,pres=None,entr=None,dhum=None,
    dliq=None,chkvals=False,chktol=_CHKTOL,airf0=None,temp0=None,
    dhum0=None,dliq0=None,chkbnd=False,mathargs=None):
    """Calculate wet air density.
    
    Calculate the density of wet air.
    
    :arg airf: Dry air mass fraction in kg/kg. If unknown, pass None
        (default) and it will be calculated.
    :type airf: float or None
    :arg temp: Temperature in K. If unknown, pass None (default) and it
        will be calculated.
    :type temp: float or None
    :arg pres: Pressure in Pa. If unknown, pass None (default) and it
        will be calculated.
    :type pres: float or None
    :arg entr: Entropy in J/kg/K.
    :type entr: float or None
    :arg dhum: Humid air density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dhum: float or None
    :arg dliq: Liquid water density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dliq: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg airf0: Initial guess for the dry fraction in kg/kg. If None
        (default) then the appropriate `_approx_*` is used.
    :type airf0: float or None
    :arg temp0: Initial guess for the temperature in K. If None
        (default) then the appropriate `_approx_*` is used.
    :type temp0: float or None
    :arg pres0: Initial guess for the pressure in Pa. If None (default)
        then the appropriate `_approx_*` is used.
    :type pres0: float or None
    :arg dhum0: Initial guess for the humid air density in kg/m3. If
        None (default) then the appropriate `_approx_*` is used.
    :type dhum0: float or None
    :arg dliq0: Initial guess for the liquid water density in kg/m3. If
        None (default) then the appropriate `_approx_*` is used.
    :type dliq0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Humid air density in kg/m3.
    :raises ValueError: If not enough values are provided.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> densityair(temp=300.,pres=1e5)
    1.14614215827
    >>> densityair(airf=.99,temp=300.)
    2.57657653270
    >>> densityair(airf=.99,pres=1e5)
    1.20675806022
    >>> densityair(airf=.99,entr=100.)
    1.45154665083
    """
    airf, temp, pres, dhum, dliq = eq_atpe(airf=airf,temp=temp,pres=pres,
        entr=entr,dhum=dhum,dliq=dliq,chkvals=chkvals,chktol=chktol,airf0=airf0,
        temp0=temp0,dhum0=dhum0,dliq0=dliq0,chkbnd=chkbnd,mathargs=mathargs)
    return dhum

def densityliq(airf=None,temp=None,pres=None,entr=None,dhum=None,
    dliq=None,chkvals=False,chktol=_CHKTOL,airf0=None,temp0=None,
    dhum0=None,dliq0=None,chkbnd=False,mathargs=None):
    """Calculate wet air liquid water density.
    
    Calculate the density of liquid water in wet air.
    
    :arg airf: Dry air mass fraction in kg/kg. If unknown, pass None
        (default) and it will be calculated.
    :type airf: float or None
    :arg temp: Temperature in K. If unknown, pass None (default) and it
        will be calculated.
    :type temp: float or None
    :arg pres: Pressure in Pa. If unknown, pass None (default) and it
        will be calculated.
    :type pres: float or None
    :arg entr: Entropy in J/kg/K.
    :type entr: float or None
    :arg dhum: Humid air density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dhum: float or None
    :arg dliq: Liquid water density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dliq: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg airf0: Initial guess for the dry fraction in kg/kg. If None
        (default) then the appropriate `_approx_*` is used.
    :type airf0: float or None
    :arg temp0: Initial guess for the temperature in K. If None
        (default) then the appropriate `_approx_*` is used.
    :type temp0: float or None
    :arg pres0: Initial guess for the pressure in Pa. If None (default)
        then the appropriate `_approx_*` is used.
    :type pres0: float or None
    :arg dhum0: Initial guess for the humid air density in kg/m3. If
        None (default) then the appropriate `_approx_*` is used.
    :type dhum0: float or None
    :arg dliq0: Initial guess for the liquid water density in kg/m3. If
        None (default) then the appropriate `_approx_*` is used.
    :type dliq0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Liquid water density in kg/m3.
    :raises ValueError: If not enough values are provided.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> densityliq(temp=300.,pres=1e5)
    996.556340389
    >>> densityliq(airf=.99,temp=300.)
    996.611581662
    >>> densityliq(airf=.99,pres=1e5)
    999.256685197
    >>> densityliq(airf=.99,entr=100.)
    998.794738784
    """
    airf, temp, pres, dhum, dliq = eq_atpe(airf=airf,temp=temp,pres=pres,
        entr=entr,dhum=dhum,dliq=dliq,chkvals=chkvals,chktol=chktol,airf0=airf0,
        temp0=temp0,dhum0=dhum0,dliq0=dliq0,chkbnd=chkbnd,mathargs=mathargs)
    return dliq

def entropy(airf=None,temp=None,pres=None,entr=None,dhum=None,dliq=None,
    chkvals=False,chktol=_CHKTOL,airf0=None,temp0=None,dhum0=None,
    dliq0=None,chkbnd=False,mathargs=None):
    """Calculate wet air entropy.
    
    Calculate the specific entropy of wet air.
    
    :arg airf: Dry air mass fraction in kg/kg. If unknown, pass None
        (default) and it will be calculated.
    :type airf: float or None
    :arg temp: Temperature in K. If unknown, pass None (default) and it
        will be calculated.
    :type temp: float or None
    :arg pres: Pressure in Pa. If unknown, pass None (default) and it
        will be calculated.
    :type pres: float or None
    :arg entr: Entropy in J/kg/K.
    :type entr: float or None
    :arg dhum: Humid air density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dhum: float or None
    :arg dliq: Liquid water density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dliq: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg airf0: Initial guess for the dry fraction in kg/kg. If None
        (default) then the appropriate `_approx_*` is used.
    :type airf0: float or None
    :arg temp0: Initial guess for the temperature in K. If None
        (default) then the appropriate `_approx_*` is used.
    :type temp0: float or None
    :arg pres0: Initial guess for the pressure in Pa. If None (default)
        then the appropriate `_approx_*` is used.
    :type pres0: float or None
    :arg dhum0: Initial guess for the humid air density in kg/m3. If
        None (default) then the appropriate `_approx_*` is used.
    :type dhum0: float or None
    :arg dliq0: Initial guess for the liquid water density in kg/m3. If
        None (default) then the appropriate `_approx_*` is used.
    :type dliq0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Entropy in J/kg/K.
    :raises ValueError: If not enough values are provided.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> entropy(temp=300.,pres=1e5)
    296.711483507
    >>> entropy(airf=.99,temp=300.)
    -41.9991507402
    >>> entropy(airf=.99,pres=1e5)
    145.863545194
    """
    airf, temp, pres, dhum, dliq = eq_atpe(airf=airf,temp=temp,pres=pres,
        entr=entr,dhum=dhum,dliq=dliq,chkvals=chkvals,chktol=chktol,airf0=airf0,
        temp0=temp0,dhum0=dhum0,dliq0=dliq0,chkbnd=chkbnd,mathargs=mathargs)
    entr = air2.entropy(airf,temp,dhum)
    return entr

def densityvap(airf=None,temp=None,pres=None,entr=None,dhum=None,
    dliq=None,chkvals=False,chktol=_CHKTOL,airf0=None,temp0=None,
    dhum0=None,dliq0=None,chkbnd=False,mathargs=None):
    """Calculate wet air water vapour density.
    
    Calculate the partial density of water vapour in wet air.
    
    :arg airf: Dry air mass fraction in kg/kg. If unknown, pass None
        (default) and it will be calculated.
    :type airf: float or None
    :arg temp: Temperature in K. If unknown, pass None (default) and it
        will be calculated.
    :type temp: float or None
    :arg pres: Pressure in Pa. If unknown, pass None (default) and it
        will be calculated.
    :type pres: float or None
    :arg entr: Entropy in J/kg/K.
    :type entr: float or None
    :arg dhum: Humid air density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dhum: float or None
    :arg dliq: Liquid water density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dliq: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg airf0: Initial guess for the dry fraction in kg/kg. If None
        (default) then the appropriate `_approx_*` is used.
    :type airf0: float or None
    :arg temp0: Initial guess for the temperature in K. If None
        (default) then the appropriate `_approx_*` is used.
    :type temp0: float or None
    :arg pres0: Initial guess for the pressure in Pa. If None (default)
        then the appropriate `_approx_*` is used.
    :type pres0: float or None
    :arg dhum0: Initial guess for the humid air density in kg/m3. If
        None (default) then the appropriate `_approx_*` is used.
    :type dhum0: float or None
    :arg dliq0: Initial guess for the liquid water density in kg/m3. If
        None (default) then the appropriate `_approx_*` is used.
    :type dliq0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Water vapour density in kg/m3.
    :raises ValueError: If not enough values are provided.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> densityvap(temp=300.,pres=1e5)
    2.56669393257e-2
    >>> densityvap(airf=.99,temp=300.)
    2.57657653270e-2
    >>> densityvap(airf=.99,pres=1e5)
    1.206758060e-2
    >>> densityvap(airf=.99,entr=100.)
    1.45154665083e-2
    """
    airf, temp, pres, dhum, dliq = eq_atpe(airf=airf,temp=temp,pres=pres,
        entr=entr,dhum=dhum,dliq=dliq,chkvals=chkvals,chktol=chktol,airf0=airf0,
        temp0=temp0,dhum0=dhum0,dliq0=dliq0,chkbnd=chkbnd,mathargs=mathargs)
    dvap = dhum * (1-airf)
    return dvap

def enthalpyevap(airf=None,temp=None,pres=None,entr=None,dhum=None,
    dliq=None,chkvals=False,chktol=_CHKTOL,airf0=None,temp0=None,
    dhum0=None,dliq0=None,chkbnd=False,mathargs=None):
    """Calculate enthalpy of evaporation.
    
    Calculate the specific enthalpy of evaporation.
    
    :arg airf: Dry air mass fraction in kg/kg. If unknown, pass None
        (default) and it will be calculated.
    :type airf: float or None
    :arg temp: Temperature in K. If unknown, pass None (default) and it
        will be calculated.
    :type temp: float or None
    :arg pres: Pressure in Pa. If unknown, pass None (default) and it
        will be calculated.
    :type pres: float or None
    :arg entr: Entropy in J/kg/K.
    :type entr: float or None
    :arg dhum: Humid air density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dhum: float or None
    :arg dliq: Liquid water density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dliq: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg airf0: Initial guess for the dry fraction in kg/kg. If None
        (default) then the appropriate `_approx_*` is used.
    :type airf0: float or None
    :arg temp0: Initial guess for the temperature in K. If None
        (default) then the appropriate `_approx_*` is used.
    :type temp0: float or None
    :arg pres0: Initial guess for the pressure in Pa. If None (default)
        then the appropriate `_approx_*` is used.
    :type pres0: float or None
    :arg dhum0: Initial guess for the humid air density in kg/m3. If
        None (default) then the appropriate `_approx_*` is used.
    :type dhum0: float or None
    :arg dliq0: Initial guess for the liquid water density in kg/m3. If
        None (default) then the appropriate `_approx_*` is used.
    :type dliq0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Enthalpy in J/kg.
    :raises ValueError: If not enough values are provided.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> enthalpyevap(temp=300.,pres=1e5)
    2434585.53919
    >>> enthalpyevap(airf=.99,temp=300.)
    2433111.29696
    >>> enthalpyevap(airf=.99,pres=1e5)
    2465656.38630
    >>> enthalpyevap(airf=.99,entr=100.)
    2458121.74961
    """
    airf, temp, pres, dhum, dliq = eq_atpe(airf=airf,temp=temp,pres=pres,
        entr=entr,dhum=dhum,dliq=dliq,chkvals=chkvals,chktol=chktol,airf0=airf0,
        temp0=temp0,dhum0=dhum0,dliq0=dliq0,chkbnd=chkbnd,mathargs=mathargs)
    fh = _air_f(0,0,0,airf,temp,dhum)
    fh_a = _air_f(1,0,0,airf,temp,dhum)
    fh_t = _air_f(0,1,0,airf,temp,dhum)
    fh_d = _air_f(0,0,1,airf,temp,dhum)
    fh_at = _air_f(1,1,0,airf,temp,dhum)
    fh_ad = _air_f(1,0,1,airf,temp,dhum)
    fh_td = _air_f(0,1,1,airf,temp,dhum)
    fh_dd = _air_f(0,0,2,airf,temp,dhum)
    comph = 2*fh_d + dhum*fh_dd
    hv = (fh - temp*fh_t + dhum*fh_d - airf*fh_a + airf*temp*fh_at
        - airf*temp*dhum * fh_td*fh_ad/comph)
    hl = flu2.enthalpy(temp,dliq)
    hvap = hv - hl
    return hvap


## Thermodynamic functions of 2 variables
def massfractionair(temp,pres,airf=None,dhum=None,dliq=None,
    chkvals=False,chktol=_CHKTOL,airf0=None,dhum0=None,dliq0=None,
    chkbnd=False,mathargs=None):
    """Calculate wet air dry fraction.
    
    Calculate the dry air mass fraction in humid air at saturation. For
    dry fractions below this, the air will be saturated.
    
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg airf: Dry air mass fraction in kg/kg. If unknown, pass None
        (default) and it will be calculated.
    :type airf: float or None
    :arg dhum: Humid air density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dhum: float or None
    :arg dliq: Liquid water density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dliq: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg airf0: Initial guess for the dry fraction in kg/kg. If None
        (default) then the appropriate `_approx_*` is used.
    :type airf0: float or None
    :arg dhum0: Initial guess for the humid air density in kg/m3. If
        None (default) then the appropriate `_approx_*` is used.
    :type dhum0: float or None
    :arg dliq0: Initial guess for the liquid water density in kg/m3. If
        None (default) then the appropriate `_approx_*` is used.
    :type dliq0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Dry air mass fraction in kg/kg.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> massfractionair(300.,1e5)
    0.977605797727
    """
    airf, __, __, dhum, dliq = eq_atpe(temp=temp,pres=pres,airf=airf,dhum=dhum,
        dliq=dliq,chkvals=chkvals,chktol=chktol,airf0=airf0,dhum0=dhum0,
        dliq0=dliq0,chkbnd=chkbnd,mathargs=mathargs)
    return airf

def condensationpressure(airf,temp,dhum=None,dliq=None,chkvals=False,
    chktol=_CHKTOL,dhum0=None,dliq0=None,chkbnd=False,mathargs=None):
    """Calculate the condensation pressure of wet air.
    
    Calculate the condensation pressure, the pressure below which liquid
    water will condense out of humid air.
    
    :arg float airf: Dry air mass fraction in kg/kg.
    :arg float temp: Temperature in K.
    :arg dhum: Humid air density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dhum: float or None
    :arg dliq: Liquid water density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dliq: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg dhum0: Initial guess for the humid air density in kg/m3. If
        None (default) then the appropriate `_approx_*` is used.
    :type dhum0: float or None
    :arg dliq0: Initial guess for the liquid water density in kg/m3. If
        None (default) then the appropriate `_approx_*` is used.
    :type dliq0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Pressure in Pa.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> condensationpressure(0.9,300.)
    23381.2332935
    """
    __, __, pres, dhum, dliq = eq_atpe(airf=airf,temp=temp,dhum=dhum,dliq=dliq,
        chkvals=chkvals,chktol=chktol,dhum0=dhum0,dliq0=dliq0,chkbnd=chkbnd,
        mathargs=mathargs)
    return pres

def dewpoint(airf,pres,temp=None,dhum=None,dliq=None,chkvals=False,
    chktol=_CHKTOL,temp0=None,dhum0=None,dliq0=None,chkbnd=False,
    mathargs=None):
    """Calculate the dew point of wet air.
    
    Calculate the dew point temperature of wet air, the temperature
    below which liquid water will condense out of humid air.
    
    :arg float airf: Dry air mass fraction in kg/kg.
    :arg float pres: Pressure in Pa.
    :arg temp: Temperature in K. If unknown, pass None (default) and it
        will be calculated.
    :type temp: float or None
    :arg dhum: Humid air density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dhum: float or None
    :arg dliq: Liquid water density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dliq: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg temp0: Initial guess for the temperature in K. If None
        (default) then the appropriate `_approx_*` is used.
    :type temp0: float or None
    :arg dhum0: Initial guess for the humid air density in kg/m3. If
        None (default) then the appropriate `_approx_*` is used.
    :type dhum0: float or None
    :arg dliq0: Initial guess for the liquid water density in kg/m3. If
        None (default) then the appropriate `_approx_*` is used.
    :type dliq0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Temperature in K.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> dewpoint(0.99,1e5)
    287.078299795
    """
    __, temp, __, dhum, dliq = eq_atpe(airf=airf,pres=pres,temp=temp,dhum=dhum,
        dliq=dliq,chkvals=chkvals,chktol=chktol,temp0=temp0,dhum0=dhum0,
        dliq0=dliq0,chkbnd=chkbnd,mathargs=mathargs)
    return temp


## Condensation level functions
def _approx_icl(airf,temp,pres,dhum):
    """Approximate TDhDl2 at ATP1.
    
    Approximate the temperature, humid air density, and liquid water
    density at the isentropic condensation level (ICL) of humid air with
    the given dry air mass fraction, in-situ temperature, and in-situ
    pressure.
    
    :arg float airf: Dry air mass fraction in kg/kg.
    :arg float temp: In-situ temperature in K.
    :arg float pres: In-situ pressure in Pa.
    :arg float dhum: In-situ humid air density in kg/m3 (unused).
    :returns: ICL temperature, humid air density, and liquid water
        density (all in SI units).
    """
    reff = airf*_RDRY + (1-airf)*_RWAT
    ceff = airf*_CDRY + (1-airf)*_CVAP
    ginv = ceff/reff
    r = (_AVL+_BVL)/(_BVL+ginv) - 1
    v = (numpy.log((1-airf)/(_EPSW*airf+1-airf)*pres/_PATM)
        + ginv*numpy.log(_TTP/temp)) / (_BVL+ginv)
    x = maths4.lamb2(v,r)
    ticl = _TTP/x
    picl = pres * (ticl/temp)**ginv
    dhicl = picl/(_RDRY*ticl) / (airf + (1-airf)/_EPSW)
    dlicl = _dliq_default(ticl,picl)
    return ticl, dhicl, dlicl

def _diff_icl(t2,dh2,dl2,airf,temp,pres,dhum):
    """Calculate disequilibrium at isentropic condensation level.
    
    Calculate both sides of the equations
    
        ICL liquid water pressure = ICL humid air pressure
        ICL liquid water chemical potential = ICL water vapour potential
        in-situ humid air entropy = ICL humid air entropy
    
    and their Jacobians with respect to the ICL temperature, humid air
    density, and liquid water density. Solving these equations gives ICL
    equilibrium values at the given dry air mass fraction, in-situ
    temperature, and in-situ pressure.
    
    :arg float t2: ICL temperature in K.
    :arg float dh2: ICL humid air density in kg/m3.
    :arg float dl2: ICL liquid water density in kg/m3.
    :arg float airf: Dry air mass fraction in kg/kg.
    :arg float temp: In-situ temperature in K.
    :arg float pres: In-situ pressure in Pa.
    :arg float dhum: In-situ humid air density in kg/m3 (unused).
    :returns: Left-hand side of the equation, right-hand side,
        Jacobian of LHS, and Jacobian of RHS.
    :rtype: tuple(array(float))
    """
    ph2 = _air_eq_pressure(0,0,0,airf,t2,dh2)
    pl2 = _flu_eq_pressure(0,0,t2,dl2)
    gv2 = _air_eq_vappot(0,0,0,airf,t2,dh2)
    gl2 = _flu_eq_chempot(0,0,t2,dl2)
    sh1 = _air_eq_entropy(0,0,0,airf,temp,dhum)
    sh2 = _air_eq_entropy(0,0,0,airf,t2,dh2)
    lhs = numpy.array([pl2, gl2, sh1])
    rhs = numpy.array([ph2, gv2, sh2])
    
    ph2_t = _air_eq_pressure(0,1,0,airf,t2,dh2)
    ph2_d = _air_eq_pressure(0,0,1,airf,t2,dh2)
    pl2_t = _flu_eq_pressure(1,0,t2,dl2)
    pl2_d = _flu_eq_pressure(0,1,t2,dl2)
    gv2_t = _air_eq_vappot(0,1,0,airf,t2,dh2)
    gv2_d = _air_eq_vappot(0,0,1,airf,t2,dh2)
    gl2_t = _flu_eq_chempot(1,0,t2,dl2)
    gl2_d = _flu_eq_chempot(0,1,t2,dl2)
    sh2_t = _air_eq_entropy(0,1,0,airf,t2,dh2)
    sh2_d = _air_eq_entropy(0,0,1,airf,t2,dh2)
    dlhs = numpy.array([[pl2_t,0.,pl2_d], [gl2_t,0.,gl2_d], [0.,0.,0.]])
    drhs = numpy.array([[ph2_t,ph2_d,0.], [gv2_t,gv2_d,0.], [sh2_t,sh2_d,0.]])
    return lhs, rhs, dlhs, drhs

def eq_icl(airf,temp,pres,dhum=None,ticl=None,dhicl=None,dlicl=None,
    chkvals=False,chktol=_CHKTOL,dhum0=None,ticl0=None,dhicl0=None,
    dlicl0=None,chkbnd=False,mathargs=None):
    """Get primary wet air variables at ICL.
    
    Get the values of all primary variables at the isentropic
    condensation level (ICL) of humid air.
    
    If the calculation has already been done, the results can be passed
    to avoid unnecessary repeat calculations. If enough values are
    passed, they will be checked for consistency if chkvals is True.
    
    :arg float airf: Dry air mass fraction in kg/kg.
    :arg float temp: In-situ temperature in K.
    :arg float pres: In-situ pressure in Pa.
    :arg dhum: In-situ humid air density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dhum: float or None
    :arg ticl: ICL temperature in K. If unknown, pass None (default) and
        it will be calculated.
    :type ticl: float or None
    :arg dhicl: ICL humid air density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dhicl: float or None
    :arg dlicl: ICL liquid water density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dlicl: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg dhum0: Initial guess for the in-situ humid air density in
        kg/m3. If None (default) then `air3a._approx_atp` is used.
    :type dhum0: float or None
    :arg ticl0: Initial guess for the ICL temperature in K. If None
        (default) then `_approx_icl` is used.
    :type ticl0: float or None
    :arg dhicl0: Initial guess for the ICL humid air density in kg/m3.
        If None (default) then `_approx_icl` is used.
    :type dhicl0: float or None
    :arg dlicl0: Initial guess for the ICL liquid water density in
        kg/m3. If None (default) then `_approx_icl` is used.
    :type dlicl0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: In-situ humid air density, ICL temperature, ICL humid air
        density, and ICL liquid water density (all in SI units).
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    """
    dhum = air3a.eq_atp(airf,temp,pres,dhum=dhum,dhum0=dhum0,mathargs=mathargs)
    if any(val is None for val in (ticl,dhicl,dlicl)):
        x0 = (ticl0,dhicl0,dlicl0)
        fargs = (airf,temp,pres,dhum)
        if mathargs is None:
            mathargs = dict()
        x1 = _newton(_diff_icl,x0,_approx_icl,fargs=fargs,**mathargs)
        ticl, dhicl, dlicl = x1
    
    _chkhumbnds(airf,temp,dhum,chkbnd=chkbnd)
    _chkhumbnds(airf,ticl,dhicl,chkbnd=chkbnd)
    _chkflubnds(ticl,dlicl,chkbnd=chkbnd)
    if not chkvals:
        return dhum, ticl, dhicl, dlicl
    
    lhs1, rhs1, __, __ = air3a._diff_atp(dhum,airf,temp,pres)
    lhs, rhs, __, __ = _diff_icl(ticl,dhicl,dlicl,airf,temp,pres,dhum)
    lhs.insert(0,lhs1)
    rhs.insert(0,rhs1)
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
    return dhum, ticl, dhicl, dlicl

def icl(airf,temp,pres,dhum=None,ticl=None,dhicl=None,dlicl=None,
    chkvals=False,chktol=_CHKTOL,dhum0=None,ticl0=None,dhicl0=None,
    dlicl0=None,chkbnd=False,mathargs=None):
    """Calculate isentropic condensation level pressure.
    
    Calculate the isentropic condensation level of humid air. Below this
    pressure, liquid water will condense out of humid air.
    
    :arg float airf: Dry air mass fraction in kg/kg.
    :arg float temp: In-situ temperature in K.
    :arg float pres: In-situ pressure in Pa.
    :arg dhum: In-situ humid air density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dhum: float or None
    :arg ticl: ICL temperature in K. If unknown, pass None (default) and
        it will be calculated.
    :type ticl: float or None
    :arg dhicl: ICL humid air density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dhicl: float or None
    :arg dlicl: ICL liquid water density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dlicl: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg dhum0: Initial guess for the in-situ humid air density in
        kg/m3. If None (default) then `air3a._approx_atp` is used.
    :type dhum0: float or None
    :arg ticl0: Initial guess for the ICL temperature in K. If None
        (default) then `_approx_icl` is used.
    :type ticl0: float or None
    :arg dhicl0: Initial guess for the ICL humid air density in kg/m3.
        If None (default) then `_approx_icl` is used.
    :type dhicl0: float or None
    :arg dlicl0: Initial guess for the ICL liquid water density in
        kg/m3. If None (default) then `_approx_icl` is used.
    :type dlicl0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: ICL pressure in Pa.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> icl(0.99,300.,1e5)
    82723.6047631
    """
    dhum, ticl, dhicl, dlicl = eq_icl(airf,temp,pres,dhum=dhum,ticl=ticl,
        dhicl=dhicl,dlicl=dlicl,chkvals=chkvals,chktol=chktol,dhum0=dhum0,
        ticl0=ticl0,dhicl0=dhicl0,dlicl0=dlicl0,chkbnd=chkbnd,mathargs=mathargs)
    picl = air2.pressure(airf,ticl,dhicl)
    return picl

def ict(airf,temp,pres,dhum=None,ticl=None,dhicl=None,dlicl=None,
    chkvals=False,chktol=_CHKTOL,dhum0=None,ticl0=None,dhicl0=None,
    dlicl0=None,chkbnd=False,mathargs=None):
    """Calculate isentropic condensation level temperature.
    
    Calculate the isentropic condensation temperature of humid air.
    Below this temperature, liquid water will condense out of humid air.
    
    :arg float airf: Dry air mass fraction in kg/kg.
    :arg float temp: In-situ temperature in K.
    :arg float pres: In-situ pressure in Pa.
    :arg dhum: In-situ humid air density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dhum: float or None
    :arg ticl: ICL temperature in K. If unknown, pass None (default) and
        it will be calculated.
    :type ticl: float or None
    :arg dhicl: ICL humid air density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dhicl: float or None
    :arg dlicl: ICL liquid water density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dlicl: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg dhum0: Initial guess for the in-situ humid air density in
        kg/m3. If None (default) then `air3a._approx_atp` is used.
    :type dhum0: float or None
    :arg ticl0: Initial guess for the ICL temperature in K. If None
        (default) then `_approx_icl` is used.
    :type ticl0: float or None
    :arg dhicl0: Initial guess for the ICL humid air density in kg/m3.
        If None (default) then `_approx_icl` is used.
    :type dhicl0: float or None
    :arg dlicl0: Initial guess for the ICL liquid water density in
        kg/m3. If None (default) then `_approx_icl` is used.
    :type dlicl0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: ICL temperature in K.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> ict(0.99,300.,1e5)
    284.200207629
    """
    dhum, ticl, dhicl, dlicl = eq_icl(airf,temp,pres,dhum=dhum,ticl=ticl,
        dhicl=dhicl,dlicl=dlicl,chkvals=chkvals,chktol=chktol,dhum0=dhum0,
        ticl0=ticl0,dhicl0=dhicl0,dlicl0=dlicl0,chkbnd=chkbnd,mathargs=mathargs)
    return ticl


## Relative humidity functions
def airffromrh_wmo(rh_wmo,temp,pres,asat=None,dhsat=None,dlsat=None,
    chkvals=False,chktol=_CHKTOL,asat0=None,dhsat0=None,dlsat0=None,
    chkbnd=False,mathargs=None):
    """Calculate dry fraction from WMO RH.
    
    Calculate the dry air mass fraction from the relative humidity. The
    relative humidity used here is defined by the WMO as:
    
        rh_wmo = [(1-airf)/airf] / [(1-asat)/asat]
    
    where asat is the dry air fraction at saturation.
    
    :arg float rh_wmo: Relative humidity, unitless.
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg asat: Dry air mass fraction at saturation in kg/kg. If unknown,
        pass None (default) and it will be calculated.
    :type asat: float or None
    :arg dhsat: Humid air density at saturation in kg/m3. If unknown,
        pass None (default) and it will be calculated.
    :type dhsat: float or None
    :arg dlsat: Liquid water density at saturation in kg/m3. If unknown,
        pass None (default) and it will be calculated.
    :type dlsat: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg asat0: Initial guess for the dry air mass fraction at
        saturation in kg/m3. If None (default) then `_approx_tp` is
        used.
    :type asat0: float or None
    :arg dhsat0: Initial guess for the humid air density at saturation
        in kg/m3. If None (default) then `_approx_tp` is used.
    :type dhsat0: float or None
    :arg dlsat0: Initial guess for the liquid water density at
        saturation in kg/m3. If None (default) then `_approx_tp` is
        used.
    :type dlsat0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Dry air mass fraction in kg/kg.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> airffromrh_wmo(0.8,300.,1e5)
    0.982004037135
    """
    asat, __, __, dhsat, dlsat = eq_atpe(temp=temp,pres=pres,airf=asat,
        dhum=dhsat,dliq=dlsat,chkvals=chkvals,chktol=chktol,airf0=asat0,
        dhum0=dhsat0,dliq0=dlsat0,chkbnd=chkbnd,mathargs=mathargs)
    airf = asat / (rh_wmo*(1-asat) + asat)
    return airf

def rhfromairf_wmo(airf,temp,pres,asat=None,dhsat=None,dlsat=None,
    chkvals=False,chktol=_CHKTOL,asat0=None,dhsat0=None,dlsat0=None,
    chkbnd=False,mathargs=None):
    """Calculate WMO RH from dry fraction.
    
    Calculate the relative humidity from the dry air mass fraction. The
    relative humidity used here is defined by the WMO as:
    
        rh_wmo = [(1-airf)/airf] / [(1-asat)/asat]
    
    where asat is the dry air fraction at saturation.
    
    :arg float airf: Dry air mass fraction in kg/kg.
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg asat: Dry air mass fraction at saturation in kg/kg. If unknown,
        pass None (default) and it will be calculated.
    :type asat: float or None
    :arg dhsat: Humid air density at saturation in kg/m3. If unknown,
        pass None (default) and it will be calculated.
    :type dhsat: float or None
    :arg dlsat: Liquid water density at saturation in kg/m3. If unknown,
        pass None (default) and it will be calculated.
    :type dlsat: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg asat0: Initial guess for the dry air mass fraction at
        saturation in kg/m3. If None (default) then `_approx_tp` is
        used.
    :type asat0: float or None
    :arg dhsat0: Initial guess for the humid air density at saturation
        in kg/m3. If None (default) then `_approx_tp` is used.
    :type dhsat0: float or None
    :arg dlsat0: Initial guess for the liquid water density at
        saturation in kg/m3. If None (default) then `_approx_tp` is
        used.
    :type dlsat0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Relative humidity, unitless.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> rhfromairf_wmo(0.99,300.,1e5)
    0.440953686019
    """
    asat, __, __, dhsat, dlsat = eq_atpe(temp=temp,pres=pres,airf=asat,
        dhum=dhsat,dliq=dlsat,chkvals=chkvals,chktol=chktol,airf0=asat0,
        dhum0=dhsat0,dliq0=dlsat0,chkbnd=chkbnd,mathargs=mathargs)
    rh_wmo = (1-airf) * asat / ((1-asat) * airf)
    return rh_wmo

def airffromrh_cct(rh_cct,temp,pres,asat=None,dhsat=None,dlsat=None,
    chkvals=False,chktol=_CHKTOL,asat0=None,dhsat0=None,dlsat0=None,
    chkbnd=False,mathargs=None):
    """Calculate dry fraction from CCT RH.
    
    Calculate the dry air mass fraction from the relative humidity. The
    relative humidity used here is defined by IUPAC/CCT as:
    
        rh_cct = vapour mol fraction / saturation vapour mol fraction.
    
    :arg float rh_cct: Relative humidity, unitless.
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg asat: Dry air mass fraction at saturation in kg/kg. If unknown,
        pass None (default) and it will be calculated.
    :type asat: float or None
    :arg dhsat: Humid air density at saturation in kg/m3. If unknown,
        pass None (default) and it will be calculated.
    :type dhsat: float or None
    :arg dlsat: Liquid water density at saturation in kg/m3. If unknown,
        pass None (default) and it will be calculated.
    :type dlsat: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg asat0: Initial guess for the dry air mass fraction at
        saturation in kg/m3. If None (default) then `_approx_tp` is
        used.
    :type asat0: float or None
    :arg dhsat0: Initial guess for the humid air density at saturation
        in kg/m3. If None (default) then `_approx_tp` is used.
    :type dhsat0: float or None
    :arg dlsat0: Initial guess for the liquid water density at
        saturation in kg/m3. If None (default) then `_approx_tp` is
        used.
    :type dlsat0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Dry air mass fraction in kg/kg.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> airffromrh_cct(0.8,300.,1e5)
    0.982133277948
    """
    asat, __, __, dhsat, dlsat = eq_atpe(temp=temp,pres=pres,airf=asat,
        dhum=dhsat,dliq=dlsat,chkvals=chkvals,chktol=chktol,airf0=asat0,
        dhum0=dhsat0,dliq0=dlsat0,chkbnd=chkbnd,mathargs=mathargs)
    xsat = convert0.air_molfractionvap(asat)
    airf = convert0.air_massfractiondry(1 - rh_cct*xsat)
    return airf

def rhfromairf_cct(airf,temp,pres,asat=None,dhsat=None,dlsat=None,
    chkvals=False,chktol=_CHKTOL,asat0=None,dhsat0=None,dlsat0=None,
    chkbnd=False,mathargs=None):
    """Calculate CCT RH from dry fraction.
    
    Calculate the relative humidity from the dry air mass fraction. The
    relative humidity used here is defined by IUPAC/CCT as:
    
        rh_cct = vapour mol fraction / saturation vapour mol fraction.
    
    :arg float airf: Dry air mass fraction in kg/kg.
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg asat: Dry air mass fraction at saturation in kg/kg. If unknown,
        pass None (default) and it will be calculated.
    :type asat: float or None
    :arg dhsat: Humid air density at saturation in kg/m3. If unknown,
        pass None (default) and it will be calculated.
    :type dhsat: float or None
    :arg dlsat: Liquid water density at saturation in kg/m3. If unknown,
        pass None (default) and it will be calculated.
    :type dlsat: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg asat0: Initial guess for the dry air mass fraction at
        saturation in kg/m3. If None (default) then `_approx_tp` is
        used.
    :type asat0: float or None
    :arg dhsat0: Initial guess for the humid air density at saturation
        in kg/m3. If None (default) then `_approx_tp` is used.
    :type dhsat0: float or None
    :arg dlsat0: Initial guess for the liquid water density at
        saturation in kg/m3. If None (default) then `_approx_tp` is
        used.
    :type dlsat0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Relative humidity, unitless.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> rhfromairf_cct(0.99,300.,1e5)
    0.449887886959
    """
    asat, __, __, dhsat, dlsat = eq_atpe(temp=temp,pres=pres,airf=asat,
        dhum=dhsat,dliq=dlsat,chkvals=chkvals,chktol=chktol,airf0=asat0,
        dhum0=dhsat0,dliq0=dlsat0,chkbnd=chkbnd,mathargs=mathargs)
    rh_cct = convert0.air_molfractionvap(airf)
    rh_cct /= convert0.air_molfractionvap(asat)
    return rh_cct

