"""Ice-air equilibrium functions.

This module provides functions to get the values of primary variables
for ice-saturated (icy) air, the equilibrium of dry air, water vapour,
and ice.

:Examples:

>>> massfractionair(temp=270.,pres=1e5)
0.997058854720
>>> temperature(airf=.997,pres=1e5)
270.234816126
>>> temperature(airf=.997,entr=100.)
266.514349350
>>> pressure(airf=.997,temp=270.)
98034.4511233
>>> pressure(airf=.997,entr=100.)
72721.4579415
>>> enthalpysubl(temp=270.,pres=1e5)
2833359.27614
>>> enthalpysubl(airf=.997,temp=270.)
2833386.54980
>>> enthalpysubl(airf=.997,pres=1e5)
2833296.51317
>>> enthalpysubl(airf=.997,entr=100.)
2834612.42351
>>> icl(0.997,300.,1e5)
64988.3931838
>>> airffromrh_wmo(0.8,270.,1e5)
0.997645698908
>>> rhfromairf_wmo(0.998,270.,1e5)
0.679365943331

:Functions:

* :func:`eq_atpe`: Calculate icy air equilibrium properties at any two of dry fraction, temperature, and pressure, or at dry fraction and entropy.
* :func:`massfractionair`: Dry air mass fraction in humid air at ice-air equilibrium.
* :func:`temperature`: Temperature at ice-air equilibrium.
* :func:`pressure`: Pressure at ice-air equilibrium.
* :func:`densityair`: Humid air density at ice-air equilibrium.
* :func:`densityvap`: Water vapour density at ice-air equilibrium.
* :func:`densityice`: Ice density at ice-air equilibrium.
* :func:`enthalpysubl`: Enthalpy of sublimation.
* :func:`condensationpressure`: Pressure at which ice will condense out of humid air.
* :func:`sublimationpressure`: Partial pressure of water vapour at the condensation pressure.
* :func:`eq_icl`: Calculate icy air equilibrium properties at the isentropic condensation level from in-situ values.
* :func:`icl`: Isentropic condensation (sublimation) pressure.
* :func:`ict`: Isentropic condensation (sublimation) temperature.
* :func:`airffromrh_wmo`: Calculate dry air fraction from relative humidity (WMO definition).
* :func:`rhfromairf_wmo`: Calculate relative humidity (WMO definition) from dry air fraction.
* :func:`airffromrh_cct`: Calculate dry air fraction from relative humidity (IUPAC/CCT definition).
* :func:`rhfromairf_cct`: Calculate relative humidity (IUPAC/CCT definition) from dry air fraction.

"""

__all__ = ['eq_atpe','massfractionair','temperature','pressure','densityair','densityvap','densityice','enthalpysubl','condensationpressure','sublimationpressure','eq_icl','icl','ict','airffromrh_wmo','rhfromairf_wmo','airffromrh_cct','rhfromairf_cct']

import warnings
import numpy
import constants0
import convert0
import ice1
import air2
import ice2
import maths3
import air3a
import maths4

_CHKTOL = constants0.CHKTOL
_MDRY = constants0.MDRY
_MWAT = constants0.MWAT
_RWAT = constants0.RWAT
_RDRY = constants0.RDRY
_PATM = constants0.PATM
_TCELS = constants0.TCELS
_TTP = constants0.TTP
_PTPE = constants0.PTPE
_LLVTP = constants0.LLVTP
_LILTP = constants0.LILTP
_CDRY = constants0.CDRY
_CVAP = constants0.CVAP
_CICE = constants0.CICE
_chkhumbnds = constants0.chkhumbnds
_chkicebnds = constants0.chkicebnds
_ice_g = ice1.ice_g
_air_f = air2.air_f
_eq_pressure = air2.eq_pressure
_eq_vappot = air2.eq_vappot
_eq_entropy = air2.eq_entropy
_newton = maths3.newton
_EPSW = _MWAT/_MDRY
_AVI = (_LLVTP + _LILTP)/(_RWAT*_TTP)
_BVI = (_CICE-_CVAP)/_RWAT
_RAB = _AVI/_BVI


## Equilibrium functions
def _approx_tp(temp,pres):
    """Approximate ADh at TP.
    
    Approximate the dry air mass fraction and humid air density for icy
    air at the given temperature and pressure. Based on an approximation
    with constant heat capacities.
    
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :returns: Dry fraction in kg/kg and humid air density in kg/m3.
    """
    earg = _AVI*(1-_TTP/temp) + _BVI*(1-_TTP/temp-numpy.log(temp/_TTP))
    pvap = _PTPE*numpy.exp(earg)
    airf = (pres-pvap)/(pres-pvap + _EPSW*pvap)
    dhum = (pres-pvap + _EPSW*pvap)/(_RDRY*temp)
    return airf, dhum

def _approx_ap(airf,pres):
    """Approximate TDh at AP.
    
    Approximate the temperature and humid air density for icy air at the
    given dry air mass fraction and pressure. Based on an approximation
    with constant heat capacities.
    
    :arg float airf: Dry air mass fraction in kg/kg.
    :arg float pres: Pressure in Pa.
    :returns: Temperature in K and humid air density in kg/m3.
    """
    pvap = pres * (1-airf)/(airf*_EPSW + 1-airf)
    v = numpy.log(pvap/_PTPE)/_BVI
    x = maths4.lamb2(v,_RAB)
    temp = _TTP/x
    dhum = (pres-pvap + _EPSW*pvap) / (_RDRY*temp)
    return temp, dhum

def _approx_at(airf,temp):
    """Approximate PDh at AT.
    
    Approximate the pressure and humid air density for icy air at the
    given dry air mass fraction and temperature. Based on an
    approximation with constant heat capacities.
    
    :arg float airf: Dry air mass fraction in kg/kg.
    :arg float temp: Temperature in K.
    :returns: Pressure in Pa and humid air density in kg/m3.
    """
    earg = _AVI*(1-_TTP/temp) + _BVI*(1-_TTP/temp-numpy.log(temp/_TTP))
    pvap = _PTPE * numpy.exp(earg)
    pres = pvap * (airf*_EPSW + 1-airf) / (1-airf)
    dhum = (pres-pvap + _EPSW*pvap) / (_RDRY*temp)
    return pres, dhum

def _approx_ae(airf,entr):
    """Approximate TPDh at AE.
    
    Approximate the temperature, pressure, and humid air density for icy
    air at the given dry air mass fraction and specific entropy. Uses an
    approximation based on constant heat capacities.
    
    :arg float airf: Dry air mass fraction in kg/kg.
    :arg float entr: Entropy in J/kg/K.
    :returns: Temperature, pressure, and humid air density (all in SI
        units).
    """
    ceff = airf*_CDRY + (1-airf)*_CVAP
    reff = airf*_RDRY + (1-airf)*_RWAT
    gam = ceff/reff
    s0 = airf*_CDRY*numpy.log(_TTP/_TCELS) - airf*_RDRY*numpy.log(_PTPE/PATM)
    s0 += (1-airf)*_LLVTP/_TTP - airf*_RDRY*numpy.log(_EPSW*airf/(1-airf))
    r = (_AVI-gam)/(gam+_BVI)
    v = (s0 - entr)/(gam+_BVI)/reff
    x = maths4.lamb2(v,r)
    temp = _TTP/x
    pres, dhum = _approx_at(airf,temp)
    return temp, pres, dhum

def _diff_tp(a,dh,temp,pres):
    """Calculate ice-air disequilibrium at TP.
    
    Calculate both sides of the equations
    
        given pressure = pressure of humid air
        chemical potential of ice = potential of water vapour
    
    and their Jacobians with respect to dry air mass fraction and humid
    air density. Solving these equations gives the mass fraction and
    humid air density at the given temperature and pressure.
    
    :arg float a: Dry air mass fraction in kg/kg.
    :arg float dh: Humid air density in kg/m3.
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :returns: Left-hand side of the equation, right-hand side,
        Jacobian of LHS, and Jacobian of RHS.
    :rtype: tuple(array(float))
    """
    ph = _eq_pressure(0,0,0,a,temp,dh)
    gv = _eq_vappot(0,0,0,a,temp,dh)
    gi = _ice_g(0,0,temp,pres)
    lhs = numpy.array([pres, gi])
    rhs = numpy.array([ph, gv])
    
    ph_a = _eq_pressure(1,0,0,a,temp,dh)
    ph_d = _eq_pressure(0,0,1,a,temp,dh)
    gv_a = _eq_vappot(1,0,0,a,temp,dh)
    gv_d = _eq_vappot(0,0,1,a,temp,dh)
    dlhs = numpy.array([[0.,0.], [0.,0.]])
    drhs = numpy.array([[ph_a,ph_d], [gv_a,gv_d]])
    return lhs, rhs, dlhs, drhs

def _diff_ap(t,dh,airf,pres):
    """Calculate ice-air disequilibrium at AP.
    
    Calculate both sides of the equations
    
        given pressure = pressure of humid air
        chemical potential of ice = potential of water vapour
    
    and their Jacobians with respect to temperature and humid air
    density. Solving these equations gives the temperature and humid air
    density at the given dry air mass fraction and pressure.

    :arg float t: Temperature in K.
    :arg float dh: Humid air density in kg/m3.
    :arg float airf: Dry air mass fraction in kg/kg.
    :arg float pres: Pressure in Pa.
    :returns: Left-hand side of the equation, right-hand side,
        Jacobian of LHS, and Jacobian of RHS.
    :rtype: tuple(array(float))
    """
    gv = _eq_vappot(0,0,0,airf,t,dh)
    gi = _ice_g(0,0,t,pres)
    ph = _eq_pressure(0,0,0,airf,t,dh)
    lhs = numpy.array([pres, gi])
    rhs = numpy.array([ph, gv])
    
    ph_t = _eq_pressure(0,1,0,airf,t,dh)
    ph_d = _eq_pressure(0,0,1,airf,t,dh)
    gi_t = _ice_g(1,0,t,pres)
    gv_t = _eq_vappot(0,1,0,airf,t,dh)
    gv_d = _eq_vappot(0,0,1,airf,t,dh)
    dlhs = numpy.array([[0.,0.], [gi_t,0.]])
    drhs = numpy.array([[ph_t,ph_d], [gv_t,gv_d]])
    return lhs, rhs, dlhs, drhs

def _diff_at(p,dh,airf,temp):
    """Calculate ice-air disequilibrium at AT.
    
    Calculate both sides of the equations
    
        given pressure = pressure of humid air
        chemical potential of ice = potential of water vapour
    
    and their Jacobians with respect to pressure and humid air density.
    Solving these equations gives the pressure and humid air density at
    the given dry air mass fraction and temperature.

    :arg float p: Pressure in Pa.
    :arg float dh: Humid air density in kg/m3.
    :arg float airf: Dry air mass fraction in kg/kg.
    :arg float temp: Temperature in K.
    :returns: Left-hand side of the equation, right-hand side,
        Jacobian of LHS, and Jacobian of RHS.
    :rtype: tuple(array(float))
    """
    gv = _eq_vappot(0,0,0,airf,temp,dh)
    gi = _ice_g(0,0,temp,p)
    ph = _eq_pressure(0,0,0,airf,temp,dh)
    lhs = numpy.array([p, gi])
    rhs = numpy.array([ph, gv])
    
    ph_d = _eq_pressure(0,0,1,airf,temp,dh)
    gi_p = _ice_g(0,1,temp,p)
    gv_d = _eq_vappot(0,0,1,airf,temp,dh)
    dlhs = numpy.array([[1.,0.], [gi_p,0.]])
    drhs = numpy.array([[0.,ph_d], [0.,gv_d]])
    return lhs, rhs, dlhs, drhs

def _diff_ae(t,p,dh,airf,entr):
    """Calculate ice-air disequilibrium at AE.
    
    Calculate both sides of the equations
    
        given pressure = pressure of humid air
        chemical potential of ice = potential of water vapour
        given entropy = entropy of humid air
    
    and their Jacobians with respect to temperature, pressure, and humid
    air density. Solving these equations gives the temperature,
    pressure, and humid air density at the given dry air mass fraction
    and specific entropy.
    
    :arg float t: Temperature in K.
    :arg float p: Pressure in Pa.
    :arg float dh: Humid air density in kg/m3.
    :arg float airf: Dry air mass fraction in kg/kg.
    :arg float entr: Entropy in J/kg/K.
    :returns: Left-hand side of the equation, right-hand side,
        Jacobian of LHS, and Jacobian of RHS.
    :rtype: tuple(array(float))
    """
    ph = _eq_pressure(0,0,0,airf,t,dh)
    gi = _ice_g(0,0,t,p)
    gv = _eq_vappot(0,0,0,airf,t,dh)
    sh = _eq_entropy(0,0,0,airf,t,dh)
    lhs = numpy.array([p, gi, entr])
    rhs = numpy.array([ph, gv, sh])
    
    ph_t = _eq_pressure(0,1,0,airf,t,dh)
    ph_d = _eq_pressure(0,0,1,airf,t,dh)
    gi_t = _ice_g(1,0,t,p)
    gi_p = _ice_g(0,1,t,p)
    gv_t = _eq_vappot(0,1,0,airf,t,dh)
    gv_d = _eq_vappot(0,0,1,airf,t,dh)
    sh_t = _eq_entropy(0,1,0,airf,t,dh)
    sh_d = _eq_entropy(0,0,1,airf,t,dh)
    dlhs = numpy.array([[0.,1.,0.], [gi_t,gi_p,0.], [0.,0.,0.]])
    drhs = numpy.array([[ph_t,0.,ph_d], [gv_t,0.,gv_d], [sh_t,0.,sh_d]])
    return lhs, rhs, dlhs, drhs

def eq_atpe(airf=None,temp=None,pres=None,entr=None,dhum=None,
    chkvals=False,chktol=_CHKTOL,airf0=None,temp0=None,pres0=None,
    dhum0=None,chkbnd=False,mathargs=None):
    """Get primary ice-air variables at ATPE.
    
    Get the values of all primary variables for icy air at any two of
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
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Dry air mass fraction, temperature, pressure, and humid
        air density (all in SI units).
    :raises ValueError: If not enough values are provided.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    """
    if sum(val is not None for val in (airf,temp,pres)) < 2:
        if airf is None or entr is None:
            errmsg = ('Not enough values were provided. Need any two of '
                '(airf,temp,pres) or both (airf,entr)')
            raise ValueError(errmsg)
        x0 = (temp0,pres0,dhum0)
        fargs = (airf,entr)
        if mathargs is None:
            mathargs = dict()
        x1 = _newton(_diff_ae,x0,_approx_ae,fargs=fargs,**mathargs)
        temp, pres, dhum = x1
    else:
        if all(val is not None for val in (airf,temp,pres)):
            if dhum is None:
                dhum = air3a.eq_atp(airf,temp,pres,dhum0=dhum0,
                    mathargs=mathargs)
        elif airf is not None and temp is not None:
            x0 = (pres0,dhum0)
            fargs = (airf,temp)
            if mathargs is None:
                mathargs = dict()
            x1 = _newton(_diff_at,x0,_approx_at,fargs=fargs,**mathargs)
            pres, dhum = x1
        elif airf is not None and pres is not None:
            x0 = (temp0,dhum0)
            fargs = (airf,pres)
            if mathargs is None:
                mathargs = dict()
            x1 = _newton(_diff_ap,x0,_approx_ap,fargs=fargs,**mathargs)
            temp, dhum = x1
        else:
            x0 = (pres0,dhum0)
            fargs = (airf,temp)
            if mathargs is None:
                mathargs = dict()
            x1 = _newton(_diff_at,x0,_approx_at,fargs=fargs,**mathargs)
            pres, dhum = x1
    
    _chkhumbnds(airf,temp,dhum,chkbnd=chkbnd)
    _chkicebnds(temp,pres,chkbnd=chkbnd)
    if not chkvals:
        return airf, temp, pres, dhum
    
    if entr is None:
        entr = air2.entropy(airf,temp,dhum)
    lhs, rhs, __, __ = _diff_ae(temp,pres,dhum,airf,entr)
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
    return airf, temp, pres, dhum


## Thermodynamic properties
def massfractionair(airf=None,temp=None,pres=None,entr=None,dhum=None,
    chkvals=False,chktol=_CHKTOL,airf0=None,temp0=None,pres0=None,
    dhum0=None,chkbnd=False,mathargs=None):
    """Calculate icy air dry air fraction.
    
    Calculate the dry air mass fraction of icy air.
    
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
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Dry air mass fraction in kg/kg.
    :raises ValueError: If not enough values are provided.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> massfractionair(temp=270.,pres=1e5)
    0.997058854720
    """
    airf, temp, pres, dhum = eq_atpe(airf=airf,temp=temp,pres=pres,entr=entr,
        dhum=dhum,chkvals=chkvals,chktol=chktol,airf0=airf0,temp0=temp0,
        pres0=pres0,dhum0=dhum0,chkbnd=chkbnd,mathargs=mathargs)
    return airf

def temperature(airf=None,temp=None,pres=None,entr=None,dhum=None,
    chkvals=False,chktol=_CHKTOL,airf0=None,temp0=None,pres0=None,
    dhum0=None,chkbnd=False,mathargs=None):
    """Calculate icy air temperature.
    
    Calculate the temperature of icy air.
    
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
    
    >>> temperature(airf=.997,pres=1e5)
    270.234816126
    >>> temperature(airf=.997,entr=100.)
    266.514349350
    """
    airf, temp, pres, dhum = eq_atpe(airf=airf,temp=temp,pres=pres,entr=entr,
        dhum=dhum,chkvals=chkvals,chktol=chktol,airf0=airf0,temp0=temp0,
        pres0=pres0,dhum0=dhum0,chkbnd=chkbnd,mathargs=mathargs)
    return temp

def pressure(airf=None,temp=None,pres=None,entr=None,dhum=None,
    chkvals=False,chktol=_CHKTOL,airf0=None,temp0=None,pres0=None,
    dhum0=None,chkbnd=False,mathargs=None):
    """Calculate icy air pressure.
    
    Calculate the pressure of icy air.
    
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
    
    >>> pressure(airf=.997,temp=270.)
    98034.4511233
    >>> pressure(airf=.997,entr=100.)
    72721.4579415
    """
    airf, temp, pres, dhum = eq_atpe(airf=airf,temp=temp,pres=pres,entr=entr,
        dhum=dhum,chkvals=chkvals,chktol=chktol,airf0=airf0,temp0=temp0,
        pres0=pres0,dhum0=dhum0,chkbnd=chkbnd,mathargs=mathargs)
    return pres

def densityair(airf=None,temp=None,pres=None,entr=None,dhum=None,
    chkvals=False,chktol=_CHKTOL,airf0=None,temp0=None,pres0=None,
    dhum0=None,chkbnd=False,mathargs=None):
    """Calculate icy air humid air density.
    
    Calculate the density of humid air in icy air.
    
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
    
    >>> densityair(temp=270.,pres=1e5)
    1.28880078014
    >>> densityair(airf=.997,temp=270.)
    1.26340801028
    >>> densityair(airf=.997,pres=1e5)
    1.28763121402
    >>> densityair(airf=.997,entr=100.)
    0.949325026119
    """
    airf, temp, pres, dhum = eq_atpe(airf=airf,temp=temp,pres=pres,entr=entr,
        dhum=dhum,chkvals=chkvals,chktol=chktol,airf0=airf0,temp0=temp0,
        pres0=pres0,dhum0=dhum0,chkbnd=chkbnd,mathargs=mathargs)
    return dhum

def densityvap(airf=None,temp=None,pres=None,entr=None,dhum=None,
    chkvals=False,chktol=_CHKTOL,airf0=None,temp0=None,pres0=None,
    dhum0=None,chkbnd=False,mathargs=None):
    """Calculate icy air vapour density.
    
    Calculate the partial density of water vapour in icy air.
    
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
    
    >>> densityvap(temp=270.,pres=1e5)
    3.79055033080e-3
    >>> densityvap(airf=.997,temp=270.)
    3.79022403085e-3
    >>> densityvap(airf=.997,pres=1e5)
    3.86289364206e-3
    >>> densityvap(airf=.997,entr=100.)
    2.84797507836e-3
    """
    airf, temp, pres, dhum = eq_atpe(airf=airf,temp=temp,pres=pres,entr=entr,
        dhum=dhum,chkvals=chkvals,chktol=chktol,airf0=airf0,temp0=temp0,
        pres0=pres0,dhum0=dhum0,chkbnd=chkbnd,mathargs=mathargs)
    dvap = dhum * (1 - airf)
    return dvap

def densityice(airf=None,temp=None,pres=None,entr=None,dhum=None,
    chkvals=False,chktol=_CHKTOL,airf0=None,temp0=None,pres0=None,
    dhum0=None,chkbnd=False,mathargs=None):
    """Calculate icy air ice density.
    
    Calculate the density of ice in icy air.
    
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
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Ice density in kg/m3.
    :raises ValueError: If not enough values are provided.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> densityice(temp=270.,pres=1e5)
    917.181167192
    >>> densityice(airf=.997,temp=270.)
    917.180955861
    >>> densityice(airf=.997,pres=1e5)
    917.147060527
    >>> densityice(airf=.997,entr=100.)
    917.681749114
    """
    airf, temp, pres, dhum = eq_atpe(airf=airf,temp=temp,pres=pres,entr=entr,
        dhum=dhum,chkvals=chkvals,chktol=chktol,airf0=airf0,temp0=temp0,
        pres0=pres0,dhum0=dhum0,chkbnd=chkbnd,mathargs=mathargs)
    dice = ice2.density(temp,pres)
    return dice

def enthalpysubl(airf=None,temp=None,pres=None,entr=None,dhum=None,
    chkvals=False,chktol=_CHKTOL,airf0=None,temp0=None,pres0=None,
    dhum0=None,chkbnd=False,mathargs=None):
    """Calculate icy air enthalpy of sublimation.
    
    Calculate the enthalpy of sublimation in icy air. This is the
    enthalpy (per unit mass of water) released when water vapour freezes
    into ice.
    
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
    
    >>> enthalpysubl(temp=270.,pres=1e5)
    2833359.27614
    >>> enthalpysubl(airf=.997,temp=270.)
    2833386.54980
    >>> enthalpysubl(airf=.997,pres=1e5)
    2833296.51317
    >>> enthalpysubl(airf=.997,entr=100.)
    2834612.42351
    """
    airf, temp, pres, dhum = eq_atpe(airf=airf,temp=temp,pres=pres,entr=entr,
        dhum=dhum,chkvals=chkvals,chktol=chktol,airf0=airf0,temp0=temp0,
        pres0=pres0,dhum0=dhum0,chkbnd=chkbnd,mathargs=mathargs)
    fh = _air_f(0,0,0,airf,temp,dhum)
    fh_a = _air_f(1,0,0,airf,temp,dhum)
    fh_t = _air_f(0,1,0,airf,temp,dhum)
    fh_d = _air_f(0,0,1,airf,temp,dhum)
    fh_at = _air_f(1,1,0,airf,temp,dhum)
    fh_ad = _air_f(1,0,1,airf,temp,dhum)
    fh_td = _air_f(0,1,1,airf,temp,dhum)
    fh_dd = _air_f(0,0,2,airf,temp,dhum)
    comph = 2*fh_d + dhum*fh_dd
    hv = fh - temp*fh_t + dhum*fh_d - airf*fh_a + airf*temp*fh_at
    hv -= airf*temp*dhum * fh_td*fh_ad/comph
    hi = ice2.enthalpy(temp,pres)
    hsubl = hv - hi
    return hsubl


## Thermodynamic functions of 2 variables
def frostpoint(airf,pres,temp=None,dhum=None,chkvals=False,
    chktol=_CHKTOL,temp0=None,dhum0=None,chkbnd=False,mathargs=None):
    """Calculate frost point of humid air.
    
    Calculate the frost point temperature, the temperature below which
    ice will condense out of humid air.
    
    :arg float airf: Dry air mass fraction in kg/kg.
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
        (default) then the appropriate `_approx_*` is used.
    :type temp0: float or None
    :arg dhum0: Initial guess for the humid air density in kg/m3. If
        None (default) then the appropriate `_approx_*` is used.
    :type dhum0: float or None
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
    
    >>> frostpoint(0.997,1e5)
    270.234816126
    """
    __, temp, __, dhum = eq_atpe(airf=airf,pres=pres,temp=temp,dhum=dhum,
        chkvals=chkvals,chktol=chktol,temp0=temp0,dhum0=dhum0,chkbnd=chkbnd,
        mathargs=mathargs)
    return temp

def condensationpressure(airf,temp,pres=None,dhum=None,chkvals=False,
    chktol=_CHKTOL,pres0=None,dhum0=None,chkbnd=False,mathargs=None):
    """Calculate the condensation pressure of humid air.
    
    Calculate the condensation pressure, the pressure below which ice
    will condense out of humid air.
    
    :arg float airf: Dry air mass fraction in kg/kg.
    :arg float temp: Temperature in K.
    :arg pres: Pressure in Pa. If unknown, pass None (default) and it
        will be calculated.
    :type pres: float or None
    :arg dhum: Humid air density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dhum: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg pres0: Initial guess for the pressure in Pa. If None (default)
        then the appropriate `_approx_*` is used.
    :type pres0: float or None
    :arg dhum0: Initial guess for the humid air density in kg/m3. If
        None (default) then the appropriate `_approx_*` is used.
    :type dhum0: float or None
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
    
    >>> ice_air_condensationpressure(0.997,270.)
    98034.4511233
    """
    __, __, pres, dhum = eq_atpe(airf=airf,temp=temp,pres=pres,dhum=dhum,
        chkvals=chkvals,chktol=chktol,pres0=pres0,dhum0=dhum0,chkbnd=chkbnd,
        mathargs=mathargs)
    return pres

def sublimationpressure(temp,pres,airf=None,dhum=None,chkvals=False,
    chktol=_CHKTOL,airf0=None,dhum0=None,chkbnd=False,mathargs=None):
    """Calculate the sublimation pressure of humid air.
    
    Calculate the sublimation pressure, the partial pressure of water
    vapour at the point at which ice will start to deposit out of humid
    air.
    
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg airf: Dry air mass fraction in kg/kg. If unknown, pass None
        (default) and it will be calculated.
    :type airf: float or None
    :arg entr: Entropy in J/kg/K.
    :type entr: float or None
    :arg dhum: Humid air density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dhum: float or None
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
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Vapour pressure in Pa.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> sublimationpressure(270.,1e5)
    472.041823975
    """
    airf, __, __, dhum = eq_atpe(temp=temp,pres=pres,airf=airf,dhum=dhum,
        chkvals=chkvals,chktol=chktol,airf0=airf0,dhum0=dhum0,chkbnd=chkbnd,
        mathargs=mathargs)
    xvap = convert0.air_molfractionvap(airf)
    pvap = pres * xvap
    return pvap


## Condensation level functions
def _approx_icl(airf,temp,pres):
    """Approximate TPDh2 at ATP1.
    
    Approximate the temperature, pressure, and humid air density at the
    isentropic condensation level (ICL) of humid air with the given
    in-situ dry air mass fraction, temperature, and pressure.
    
    :arg float airf: Dry air mass fraction in kg/kg.
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :returns: ICL temperature, pressure, and humid air density (all in
        SI units).
    """
    reff = airf*_RDRY + (1-airf)*_RWAT
    ceff = airf*_CDRY + (1-airf)*_CVAP
    ginv = ceff/reff
    r = (_AVI+_BVI)/(_BVI+ginv) - 1
    v = (numpy.log((1-airf)/(_EPSW*airf+1-airf)*pres/PATM)
        + ginv*numpy.log(_TTP/temp)) / (_BVI+ginv)
    x = maths4.lamb2(v,r)
    ticl = _TTP/x
    picl = pres * (ticl/temp)**gieff
    dhicl = pres/(_RDRY*temp) / (airf + (1-airf)/_EPSW)
    return ticl, picl, dhicl

def _diff_icl(t2,p2,dh2,airf,temp,pres,dhum):
    """Calculate ICL disequilibrium.
    
    Calculate both sides of the equations
    
        given icl pressure = final pressure of humid air
        chemical potential of ice at icl = vapour potential at icl
        in-situ humid air entropy = humid air entropy at icl
    
    and their Jacobians with respect to the isentropic condensation
    level (ICL) temperature, pressure, and humid air density. Solving
    these equations gives the temperature, pressure, and humid air
    density at the ICL for humid air with the given in-situ dry air mass
    fraction, temperature, and pressure.
    
    :arg float t2: ICL temperature in K.
    :arg float p2: ICL pressure in Pa.
    :arg float dh2: ICL humid air density in kg/m3.
    :arg float airf: Dry air mass fraction in kg/kg.
    :arg float temp: In-situ temperature in K.
    :arg float pres: In-situ pressure in Pa.
    :arg float dhum: In-situ humid air density in kg/m3.
    :returns: Left-hand side of the equation, right-hand side,
        Jacobian of LHS, and Jacobian of RHS.
    :rtype: tuple(array(float))
    """
    ph2 = _eq_pressure(0,0,0,airf,t2,dh2)
    gi2 = _ice_g(0,0,t2,p2)
    gv2 = _eq_vappot(0,0,0,airf,t2,dh2)
    sh1 = _eq_entropy(0,0,0,airf,temp,dhum)
    sh2 = _eq_entropy(0,0,0,airf,t2,dh2)
    lhs = numpy.array([p2, gi2, sh1])
    rhs = numpy.array([ph2, gv2, sh2])
    
    ph2_t = _eq_pressure(0,1,0,airf,t2,dh2)
    ph2_d = _eq_pressure(0,0,1,airf,t2,dh2)
    gi2_t = _ice_g(1,0,t2,p2)
    gi2_p = _ice_g(0,1,t2,p2)
    gv2_t = _eq_vappot(0,1,0,airf,t2,dh2)
    gv2_d = _eq_vappot(0,0,1,airf,t2,dh2)
    sh2_t = _eq_entropy(0,1,0,airf,t2,dh2)
    sh2_d = _eq_entropy(0,0,1,airf,t2,dh2)
    dlhs = numpy.array([[0.,1.,0.], [gi2_t,gi2_p,0.], [0.,0.,0.]])
    drhs = numpy.array([[ph2_t,0.,ph2_d], [gv2_t,0.,gv2_d], [sh2_t,0.,sh2_d]])
    return lhs, rhs, dlhs, drhs

def eq_icl(airf,temp,pres,dhum=None,ticl=None,picl=None,dhicl=None,
    chkvals=False,chktol=_CHKTOL,dhum0=None,ticl0=None,picl0=None,
    dhicl0=None,chkbnd=False,mathargs=None):
    """Get primary icy air variables at the ICL.
    
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
    :arg picl: ICL pressure in Pa. If unknown, pass None (default) and
        it will be calculated.
    :type picl: float or None
    :arg dhicl: ICL humid air density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dhicl: float or None
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
    :arg picl0: Initial guess for the ICL pressure in Pa. If None
        (default) then `_approx_icl` is used.
    :type picl0: float or None
    :arg dhicl0: Initial guess for the ICL humid air density in kg/m3.
        If None (default) then `_approx_icl` is used.
    :type dhicl0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: In-situ humid air density, ICL temperature, ICL pressure,
        and ICL humid air density (all in SI units).
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    """
    dhum = air3a.eq_atp(airf,temp,pres,dhum=dhum,chkvals=chkvals,chktol=chktol,
        dhum0=dhum0,chkbnd=chkbnd,mathargs=mathargs)
    if any(val is None for val in (ticl,picl,dhicl)):
        x0 = (ticl0,picl0,dhicl0)
        fargs = (airf,temp,pres,dhum)
        if mathargs is None:
            mathargs = dict()
        x1 = _newton(_diff_icl,x0,_approx_icl,fargs=fargs,**mathargs)
        ticl, picl, dhicl = x1
    
    _chkhumbnds(airf,temp,dhum,chkbnd=chkbnd)
    _chkhumbnds(airf,ticl,dhicl,chkbnd=chkbnd)
    _chkicebnds(ticl,picl,chkbnd=chkbnd)
    if not chkvals:
        return dhum, ticl, picl, dhicl
    
    lhs1, rhs1, __, __ = air3a._diff_atp(dhum,airf,temp,pres)
    lhs2, rhs2, __, __ = _diff_icl(ticl,picl,dhicl,airf,temp,pres,dhum)
    lhs = [lhs1] + lhs2.tolist()
    rhs = [rhs1] + rhs2.tolist()
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
    return dhum, ticl, picl, dhicl

def icl(airf,temp,pres,dhum=None,ticl=None,picl=None,dhicl=None,
    chkvals=False,chktol=_CHKTOL,dhum0=None,ticl0=None,picl0=None,
    dhicl0=None,chkbnd=False,mathargs=None):
    """Calculate isentropic condensation pressure.
    
    Calculate the pressure at the isentropic condensation (deposition)
    level. Below this pressure, ice will deposit out of humid air.
    
    :arg float airf: Dry air mass fraction in kg/kg.
    :arg float temp: In-situ temperature in K.
    :arg float pres: In-situ pressure in Pa.
    :arg dhum: In-situ humid air density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dhum: float or None
    :arg ticl: ICL temperature in K. If unknown, pass None (default) and
        it will be calculated.
    :type ticl: float or None
    :arg picl: ICL pressure in Pa. If unknown, pass None (default) and
        it will be calculated.
    :type picl: float or None
    :arg dhicl: ICL humid air density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dhicl: float or None
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
    :arg picl0: Initial guess for the ICL pressure in Pa. If None
        (default) then `_approx_icl` is used.
    :type picl0: float or None
    :arg dhicl0: Initial guess for the ICL humid air density in kg/m3.
        If None (default) then `_approx_icl` is used.
    :type dhicl0: float or None
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
    
    >>> icl(0.997,300.,1e5)
    64988.3931838
    """
    dhum, ticl, picl, dhicl = eq_icl(airf,temp,pres,dhum=dhum,ticl=ticl,
        picl=picl,dhicl=dhicl,chkvals=chkvals,chktol=chktol,dhum0=dhum0,
        ticl0=ticl0,picl0=picl0,dhicl0=dhicl0,chkbnd=chkbnd,mathargs=mathargs)
    return picl

def ict(airf,temp,pres,dhum=None,ticl=None,picl=None,dhicl=None,
    chkvals=False,chktol=_CHKTOL,dhum0=None,ticl0=None,picl0=None,
    dhicl0=None,chkbnd=False,mathargs=None):
    """Calculate isentropic condensation temperature.
    
    Calculate the temperature at the isentropic condensation
    (deposition) level. Below this temperature, ice will deposit out of
    humid air.
    
    :arg float airf: Dry air mass fraction in kg/kg.
    :arg float temp: In-situ temperature in K.
    :arg float pres: In-situ pressure in Pa.
    :arg dhum: In-situ humid air density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dhum: float or None
    :arg ticl: ICL temperature in K. If unknown, pass None (default) and
        it will be calculated.
    :type ticl: float or None
    :arg picl: ICL pressure in Pa. If unknown, pass None (default) and
        it will be calculated.
    :type picl: float or None
    :arg dhicl: ICL humid air density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dhicl: float or None
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
    :arg picl0: Initial guess for the ICL pressure in Pa. If None
        (default) then `_approx_icl` is used.
    :type picl0: float or None
    :arg dhicl0: Initial guess for the ICL humid air density in kg/m3.
        If None (default) then `_approx_icl` is used.
    :type dhicl0: float or None
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
    
    >>> ict(0.997,300.,1e5)
    265.224998411
    """
    dhum, ticl, picl, dhicl = eq_icl(airf,temp,pres,dhum=dhum,ticl=ticl,
        picl=picl,dhicl=dhicl,chkvals=chkvals,chktol=chktol,dhum0=dhum0,
        ticl0=ticl0,picl0=picl0,dhicl0=dhicl0,chkbnd=chkbnd,mathargs=mathargs)
    return ticl


## Relative humidity functions
def ice_air_airffromrh_wmo(rh_wmo,temp,pres,asat=None,dhsat=None,chkvals=False,
    chktol=_CHKTOL,asat0=None,dhsat0=None,chkbnd=False,mathargs=None):
    """Calculate air fraction from WMO relative humidity.
    
    Calculate the dry air mass fraction from the relative humidity. The
    relative humidity used here is defined by the WMO as:
        rh_wmo = (1-A)/A / ((1-A_sat/A_sat)).
    
    Args:
        rh_wmo (float): Relative humidity, unitless.
        temp (float): Temperature in K.
        pres (float): Pressure in Pa.
        asat (float or None, optional): Dry air mass fraction at saturation in
            kg/kg. If unknown, pass None (default) and it will be calculated.
        dhsat (float or None, optional): Humid air density at saturation in
            kg/m3. If unknown, pass None (default) and it will be calculated.
        chkvals (boolean, optional): If True (default) and results are given,
            warnings will be raised if the results are out of equilibrium.
        chktol (float, optional): Relative tolerance to use for raising warnings
            when results are given (default 1e-8).
        asat0 (float or None, optional): Initial guess for the dry air mass
            fraction at saturation in kg/kg. If None (default) then approx_tp is
            used.
        dhsat0 (float or None, optional): Initial guess for the humid air
            density at saturation in kg/m3. If None (default) then approx_tp is
            used.
        chkbnd (boolean, optional): If True (default), warnings are raised when
            the given values are valid but outside the recommended bounds.
        mathargs (dict or None, optional): Keyword arguments to pass to the
            root-finder _newton (e.g. maxiter, rtol). If None (default), no
            arguments are passed and all parameters will be default.
    
    Returns:
        airf (float): Dry air mass fraction in kg/kg.
    
    Examples:
        >>> ice_air_airffromrh_wmo(0.8,270.,1e5)
        0.997645698908
    """
    
    res = ice_air_massfraction_air(airf=asat,temp=temp,pres=pres,dhum=dhsat,
        chkvals=chkvals,chktol=chktol,airf0=asat0,dhum0=dhsat0,chkbnd=chkbnd,
        mathargs=mathargs)
    if asat is None:
        asat = res
    airf = asat / (rh_wmo*(1-asat) + asat)
    return airf

def ice_air_rhfromairf_wmo(airf,temp,pres,asat=None,dhsat=None,chkvals=False,
    chktol=_CHKTOL,asat0=None,dhsat0=None,chkbnd=False,mathargs=None):
    """Calculate WMO relative humidity from air fraction.
    
    Calculate the relative humidity from the dry air mass fraction. The
    relative humidity used here is defined by the WMO as:
        rh_wmo = (1-A)/A / ((1-A_sat/A_sat)).
    
    Args:
        airf (float): Dry air mass fraction in kg/kg.
        temp (float): Temperature in K.
        pres (float): Pressure in Pa.
        asat (float or None, optional): Dry air mass fraction at saturation in
            kg/kg. If unknown, pass None (default) and it will be calculated.
        dhsat (float or None, optional): Humid air density at saturation in
            kg/m3. If unknown, pass None (default) and it will be calculated.
        chkvals (boolean, optional): If True (default) and results are given,
            warnings will be raised if the results are out of equilibrium.
        chktol (float, optional): Relative tolerance to use for raising warnings
            when results are given (default 1e-8).
        asat0 (float or None, optional): Initial guess for the dry air mass
            fraction at saturation in kg/kg. If None (default) then approx_tp is
            used.
        dhsat0 (float or None, optional): Initial guess for the humid air
            density at saturation in kg/m3. If None (default) then approx_tp is
            used.
        chkbnd (boolean, optional): If True (default), warnings are raised when
            the given values are valid but outside the recommended bounds.
        mathargs (dict or None, optional): Keyword arguments to pass to the
            root-finder _newton (e.g. maxiter, rtol). If None (default), no
            arguments are passed and all parameters will be default.
    
    Returns:
        rh_wmo (float): Relative humidity, unitless.
    
    Examples:
        >>> ice_air_rhfromairf_wmo(0.998,270.,1e5)
        0.679365943331
    """
    
    res = ice_air_massfraction_air(airf=asat,temp=temp,pres=pres,dhum=dhsat,
        chkvals=chkvals,chktol=chktol,airf0=asat0,dhum0=dhsat0,chkbnd=chkbnd,
        mathargs=mathargs)
    if asat is None:
        asat = res
    rh_wmo = (1-airf) * asat / ((1-asat) * airf)
    return rh_wmo

def ice_air_airffromrh_cct(rh_cct,temp,pres,asat=None,dhsat=None,chkvals=False,
    chktol=_CHKTOL,asat0=None,dhsat0=None,chkbnd=False,mathargs=None):
    """Calculate CCT relative humidity from air fraction.
    
    Calculate the relative humidity from the dry air mass fraction. The
    relative humidity used here is defined by the CCT/IUPAC as:
        rh_cct = vapour mol fraction / saturation vapour mol fraction.
    
    Args:
        rh_cct (float): Relative humidity, unitless.
        temp (float): Temperature in K.
        pres (float): Pressure in Pa.
        asat (float or None, optional): Dry air mass fraction at saturation in
            kg/kg. If unknown, pass None (default) and it will be calculated.
        dhsat (float or None, optional): Humid air density at saturation in
            kg/m3. If unknown, pass None (default) and it will be calculated.
        chkvals (boolean, optional): If True (default) and results are given,
            warnings will be raised if the results are out of equilibrium.
        chktol (float, optional): Relative tolerance to use for raising warnings
            when results are given (default 1e-8).
        asat0 (float or None, optional): Initial guess for the dry air mass
            fraction at saturation in kg/kg. If None (default) then approx_tp is
            used.
        dhsat0 (float or None, optional): Initial guess for the humid air
            density at saturation in kg/m3. If None (default) then approx_tp is
            used.
        chkbnd (boolean, optional): If True (default), warnings are raised when
            the given values are valid but outside the recommended bounds.
        mathargs (dict or None, optional): Keyword arguments to pass to the
            root-finder _newton (e.g. maxiter, rtol). If None (default), no
            arguments are passed and all parameters will be default.
    
    Returns:
        airf (float): Dry air mass fraction in kg/kg.
    
    Examples:
        >>> ice_air_airffromrh_cct(0.8,270.,1e5)
        0.997647924743
    """
    
    res = ice_air_massfraction_air(airf=asat,temp=temp,pres=pres,dhum=dhsat,
        chkvals=chkvals,chktol=chktol,airf0=asat0,dhum0=dhsat0,chkbnd=chkbnd,
        mathargs=mathargs)
    if asat is None:
        asat = res
    xsat = convert0.air_molfractionvap(asat)
    airf = convert0.air_massfractiondry(1 - rh_cct*xsat)
    return airf

def ice_air_rhfromairf_cct(airf,temp,pres,asat=None,dhsat=None,chkvals=False,
    chktol=_CHKTOL,asat0=None,dhsat0=None,chkbnd=False,mathargs=None):
    """Calculate air fraction from CCT relative humidity.
    
    Calculate the dry air mass fraction from the relative humidity. The
    relative humidity used here is defined by the CCT/IUPAC as:
        rh_cct = vapour mol fraction / saturation vapour mol fraction.
    
    Args:
        airf (float): Dry air mass fraction in kg/kg.
        temp (float): Temperature in K.
        pres (float): Pressure in Pa.
        asat (float or None, optional): Dry air mass fraction at saturation in
            kg/kg. If unknown, pass None (default) and it will be calculated.
        dhsat (float or None, optional): Humid air density at saturation in
            kg/m3. If unknown, pass None (default) and it will be calculated.
        chkvals (boolean, optional): If True (default) and results are given,
            warnings will be raised if the results are out of equilibrium.
        chktol (float, optional): Relative tolerance to use for raising warnings
            when results are given (default 1e-8).
        asat0 (float or None, optional): Initial guess for the dry air mass
            fraction at saturation in kg/kg. If None (default) then approx_tp is
            used.
        dhsat0 (float or None, optional): Initial guess for the humid air
            density at saturation in kg/m3. If None (default) then approx_tp is
            used.
        chkbnd (boolean, optional): If True (default), warnings are raised when
            the given values are valid but outside the recommended bounds.
        mathargs (dict or None, optional): Keyword arguments to pass to the
            root-finder _newton (e.g. maxiter, rtol). If None (default), no
            arguments are passed and all parameters will be default.
    
    Returns:
        rh_cct (float): Relative humidity, unitless.
    
    Examples:
        >>> ice_air_rhfromairf_cct(0.998,270.,1e5)
        0.680395740553
    """
    
    res = ice_air_massfraction_air(airf=asat,temp=temp,pres=pres,dhum=dhsat,
        chkvals=chkvals,chktol=chktol,airf0=asat0,dhum0=dhsat0,chkbnd=chkbnd,
        mathargs=mathargs)
    if asat is None:
        asat = res
    rh_cct = convert0.air_molfractionvap(airf)
    rh_cct /= convert0.air_molfractionvap(asat)
    return rh_cct


