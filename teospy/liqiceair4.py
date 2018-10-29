"""Liquid water-ice-humid air equilibrium functions.

This module provides functions to get the values of primary variables
for humid air in equilibrium with both liquid water and ice (wet-icy
air). This requires temperatures very close to the triple point
(273.16 K). The primary variables can be either the mass fractions of
dry air, liquid water, and ice; or the dry air mass fraction, total
entropy, and wet fraction of the condensates.

:Examples:

>>> pressure(airf=.99)
38338.9622424
>>> temperature(airf=.99)
273.157198087
>>> airfraction(temp=273.155)
0.994366063923
>>> pressure(temp=273.155)
67931.60108
>>> airfraction(pres=1e4)
0.961024307544
>>> temperature(pres=1e4)
273.159302793
>>> pressure(wair=.1,wliq=.2,wice=.3)
706.817425301
>>> temperature(wair=.1,wliq=.2,wice=.3)
273.159992933
>>> density(.1,wliq=.2,wice=.3)
0.0121364037568
>>> enthalpy(.1,wliq=.2,wice=.3)
900361.135280
>>> entropy(.1,wliq=.2,wice=.3)
3496.16306903
>>> airfraction(wair=.99,entr=0.,wetf=.5)
0.996583352944
>>> pressure(wair=.99,entr=0.,wetf=.5)
112016.075795
>>> temperature(wair=.99,entr=0.,wetf=.5)
273.151724970
>>> density(.99,entr=0.,wetf=.5)
1.43611528680
>>> enthalpy(.99,entr=0.,wetf=.5)
7356.12943724
>>> liquidfraction(.99,entr=0.,wetf=.5)
3.30296152581e-3
>>> solidfraction(.99,entr=0.,wetf=.5)
3.30296152581e-3
>>> vapourfraction(.99,entr=0.,wetf=.5)
3.39407694837e-3
>>> iml(.99,100.)
81605.5557729
>>> ifl(.99,100.)
83234.7314358

:Functions:

* :func:`eq_atp`: Calculate wet-icy air equilibrium properties from any
  of the humid air dry fraction, temperature, or pressure.
* :func:`eq_wefli`: Calculate wet-icy air equilibrium properties from
  either the mass fractions of dry air, liquid water, and ice; or from
  the dry air fraction, entropy, and wet fraction of the condensates.
* :func:`eq_all`: Calculate wet-icy air equilibrium properties. This
  function is just a common wrapper for `eq_atp` and `eq_wefli`.
* :func:`airfraction`: Humid air dry fraction at equilibrium.
* :func:`pressure`: Pressure at equilibrium.
* :func:`temperature`: Temperature at equilibrium.
* :func:`density`: Total wet-icy air density.
* :func:`dryairfraction`: Total dry fraction in wet-icy air.
* :func:`enthalpy`: Specific enthalpy of wet-icy air.
* :func:`entropy`: Specific entropy of wet-icy air.
* :func:`liquidfraction`: Mass fraction of liquid water in wet-icy air.
* :func:`solidfraction`: Mass fraction of ice in wet-icy air.
* :func:`vapourfraction`: Mass fraction of water vapour in wet-icy air.
* :func:`iml`: Isentropic melting level of wet-icy air.
* :func:`ifl`: Isentropic freezing level of wet-icy air.

"""

__all__ = ['eq_atp','eq_wefli','eq_all',
    'airfraction','pressure','temperature',
    'density','dryairfraction','enthalpy','entropy','liquidfraction',
    'solidfraction','vapourfraction',
    'iml','ifl']

import numpy
import warnings
import constants0
import flu1
import ice1
import air2
import flu2
import ice2
import maths3
import flu3a
import iceliq4

_CHKTOL = constants0.CHKTOL
_RWAT = constants0.RWAT
_RDRY = constants0.RDRY
_PATM = constants0.PATM
_TCELS = constants0.TCELS
_TTP = constants0.TTP
_PTPE = constants0.PTPE
_LLVTP = constants0.LLVTP
_LILTP = -constants0.LILTP
_DLTP = constants0.DLTP
_DITP = constants0.DITP
_CDRY = constants0.CDRY
_CVAP = constants0.CVAP
_CLIQ = constants0.CLIQ
_CICE = constants0.CICE
_EPSW = _RDRY/_RWAT
_AVL = _LLVTP / (_RWAT*_TTP)
_ALI = _LILTP/(_DITP**(-1) - _DLTP**(-1))/_PTPE

_chkflubnds = constants0.chkflubnds
_chkhumbnds = constants0.chkhumbnds
_chkicebnds = constants0.chkicebnds
_flu_f = flu1.flu_f
_ice_g = ice1.ice_g
_air_f = air2.air_f
_air_eq_pressure = air2.eq_pressure
_air_eq_vappot = air2.eq_vappot
_flu_eq_chempot = flu2.eq_chempot
_flu_eq_pressure = flu2.eq_pressure
_newton = maths3.newton
_dliq_default = flu3a._dliq_default


### Equilibrium functions
def _approx_a(airf):
    """Approximate TPDhDl at A.
    
    Approximate the temperature, pressure, humid air density, and liquid
    water density of wet-icy air at the given humid air dry fraction.
    
    :arg float airf: Humid air dry fraction in kg/kg.
    :returns: Temperature, pressure, humid air density, and liquid water
        density (all in SI units).
    """
    temp = _TTP * (1 - _EPSW*airf/(1-airf)/_ALI)
    pres = _PTPE*(_EPSW*airf + 1-airf)/(1-airf)
    dhum = pres/(_RDRY*temp) / (airf + (1-airf)/_EPSW)
    dliq = _dliq_default(temp,pres)
    return temp, pres, dhum, dliq

def _approx_t(temp):
    """Approximate APDhDl at T.
    
    Approximate the humid air dry fraction, pressure, humid air density,
    and liquid water density of wet-icy air at the given temperature.
    
    :arg float temp: Temperature in K.
    :returns: Humid air dry fraction, pressure, humid air density, and
        liquid water density (all in SI units).
    """
    #pres = _PTPE*(1 + _ALI*(1-temp/_TTP))
    #dliq = _dliq_default(temp,pres)
    pres, dliq = iceliq4._approx_t(temp)
    pvap = _PTPE * (1 - _AVL*(1 - temp/_TTP))
    airf = (pres-pvap)/(pres-pvap + _EPSW*pvap)
    dhum = pres/(_RDRY*temp) / (airf + (1-airf)/_EPSW)
    return airf, pres, dhum, dliq

def _approx_p(pres):
    """Approximate ATDhDl at P.
    
    Approximate the humid air dry fraction, temperature, humid air
    density, and liquid water density of wet-icy air at the given
    pressure.
    
    :arg float pres: Pressure in Pa.
    :returns: Humid air dry fraction, temperature, humid air density,
        and liquid water density (all in SI units).
    """
    #temp = min(_TTP*(1 - (pres/_PTPE-1)/_ALI), _TTP-_CHKTOL)
    #dliq = _dliq_default(temp,pres)
    temp, dliq = iceliq4._approx_p(pres)
    pvap = _PTPE * (1 - _AVL*(1 - temp/_TTP))
    airf = (pres-pvap)/(pres-pvap + _EPSW*pvap)
    dhum = pres/(_RDRY*temp) / (airf + (1-airf)/_EPSW)
    return airf, temp, dhum, dliq

def _approx_wef(wair,entr,wetf):
    """Approximate ATPDhDl at WEF.
    
    Approximate the humid air dry fraction, temperature, pressure, humid
    air density, and liquid water density of wet-icy air at the given
    total dry fraction, specific entropy, and wet fraction of
    condensate.
    
    :arg float wair: Total dry fraction in kg/kg.
    :arg float entr: Entropy in J/kg/K.
    :arg float wetf: Wet fraction in kg/kg.
    :returns: Humid air dry fraction, temperature, pressure, humid air
        density, and liquid water density (all in SI units).
    """
    earg = ((entr + (1-wair)*(1-wetf)*_LILTP/_TTP
        - wair*_RDRY*numpy.log(_PATM/_PTPE)) / (wair*_RDRY))
    coeff = (_LLVTP + (1-wetf)*_LILTP)/(_RWAT*_TTP)
    w = coeff * numpy.exp(earg)
    z = w - w**2 + 1.5*w**3
    pres = _PTPE + coeff*_PTPE/z
    airf, temp, dhum, dliq = _approx_p(pres)
    return airf, temp, pres, dhum, dliq

def _diff_a(t,p,dh,dl,airf):
    """Calculate wet-icy air disequilibrium at A.
    
    Calculate both sides of the equations
    
        given pressure = pressure in humid air
        given pressure = pressure of liquid water
        chemical potential of liquid water = potential of ice
        chemical potential of liquid water = potential of water vapour
    
    and their Jacobians with respect to temperature, pressure, humid air
    density, and liquid water density. Solving these equations gives
    equilibrium values at the given humid air dry fraction.
    
    :arg float t: Temperature in K.
    :arg float p: Pressure in Pa.
    :arg float dh: Humid air density in kg/m3.
    :arg float dl: Liquid water density in kg/m3.
    :arg float airf: Humid air dry fraction in kg/kg.
    :returns: Left-hand side of the equation, right-hand side,
        Jacobian of LHS, and Jacobian of RHS.
    :rtype: tuple(array(float))
    """
    ph = _air_eq_pressure(0,0,0,airf,t,dh)
    pl = _flu_eq_pressure(0,0,t,dl)
    gl = _flu_eq_chempot(0,0,t,dl)
    gi = _ice_g(0,0,t,p)
    gv = _air_eq_vappot(0,0,0,airf,t,dh)
    lhs = numpy.array([p, p, gl, gl])
    rhs = numpy.array([ph, pl, gi, gv])
    
    ph_t = _air_eq_pressure(0,1,0,airf,t,dh)
    ph_d = _air_eq_pressure(0,0,1,airf,t,dh)
    pl_t = _flu_eq_pressure(1,0,t,dl)
    pl_d = _flu_eq_pressure(0,1,t,dl)
    gl_t = _flu_eq_chempot(1,0,t,dl)
    gl_d = _flu_eq_chempot(0,1,t,dl)
    gi_t = _ice_g(1,0,t,p)
    gi_p = _ice_g(0,1,t,p)
    gv_t = _air_eq_vappot(0,1,0,airf,t,dh)
    gv_d = _air_eq_vappot(0,0,1,airf,t,dh)
    dlhs = numpy.array([[0.,1.,0.,0.], [0.,1.,0.,0.], [gl_t,0.,0.,gl_d],
        [gl_t,0.,0.,gl_d]])
    drhs = numpy.array([[ph_t,0.,ph_d,0.], [pl_t,0.,0.,pl_d],
        [gi_t,gi_p,0.,0.], [gv_t,0.,gv_d,0.]])
    return lhs, rhs, dlhs, drhs

def _diff_t(a,p,dh,dl,temp):
    """Calculate wet-icy air disequilibrium at T.
    
    Calculate both sides of the equations
    
        given pressure = pressure in humid air
        given pressure = pressure of liquid water
        chemical potential of liquid water = potential of ice
        chemical potential of liquid water = potential of water vapour
    
    and their Jacobians with respect to humid air dry fraction,
    pressure, humid air density, and liquid water density. Solving these
    equations gives equilibrium values at the given temperature.

    :arg float a: Humid air dry fraction in kg/kg.
    :arg float p: Pressure in Pa.
    :arg float dh: Humid air density in kg/m3.
    :arg float dl: Liquid water density in kg/m3.
    :arg float temp: Temperature in K.
    :returns: Left-hand side of the equation, right-hand side,
        Jacobian of LHS, and Jacobian of RHS.
    :rtype: tuple(array(float))
    """
    ph = _air_eq_pressure(0,0,0,a,temp,dh)
    pl = _flu_eq_pressure(0,0,temp,dl)
    gl = _flu_eq_chempot(0,0,temp,dl)
    gi = _ice_g(0,0,temp,p)
    gv = _air_eq_vappot(0,0,0,a,temp,dh)
    lhs = numpy.array([p, p, gl, gl])
    rhs = numpy.array([ph, pl, gi, gv])
    
    ph_a = _air_eq_pressure(1,0,0,a,temp,dh)
    ph_d = _air_eq_pressure(0,0,1,a,temp,dh)
    pl_d = _flu_eq_pressure(0,1,temp,dl)
    gl_d = _flu_eq_chempot(0,1,temp,dl)
    gi_p = _ice_g(0,1,temp,p)
    gv_a = _air_eq_vappot(1,0,0,a,temp,dh)
    gv_d = _air_eq_vappot(0,0,1,a,temp,dh)
    dlhs = numpy.array([[0.,1.,0.,0.], [0.,1.,0.,0.], [0.,0.,0.,gl_d],
        [0.,0.,0.,gl_d]])
    drhs = numpy.array([[ph_a,0.,ph_d,0.], [0.,0.,0.,pl_d], [0.,gi_p,0.,0.],
        [gv_a,0.,gv_d,0.]])
    return lhs, rhs, dlhs, drhs

def _diff_p(a,t,dh,dl,pres):
    """Calculate wet-icy air disequilibrium at pressure.
    
    Calculate both sides of the equations
    
        given pressure = pressure in humid air
        given pressure = pressure of liquid water
        chemical potential of liquid water = potential of ice
        chemical potential of liquid water = potential of water vapour
    
    and their Jacobians with respect to dry air mass fraction in humid
    air, temperature, humid air density, and liquid water density.
    Solving these equations gives equilibrium values at the given
    pressure.
    
    :arg float a: Humid air dry fraction in kg/kg.
    :arg float t: Temperature in K.
    :arg float dh: Humid air density in kg/m3.
    :arg float dl: Liquid water density in kg/m3.
    :arg float pres: Pressure in Pa.
    :returns: Left-hand side of the equation, right-hand side,
        Jacobian of LHS, and Jacobian of RHS.
    :rtype: tuple(array(float))
    """
    ph = _air_eq_pressure(0,0,0,a,t,dh)
    pl = _flu_eq_pressure(0,0,t,dl)
    gl = _flu_eq_chempot(0,0,t,dl)
    gi = _ice_g(0,0,t,pres)
    gv = _air_eq_vappot(0,0,0,a,t,dh)
    lhs = numpy.array([pres, pres, gl, gl])
    rhs = numpy.array([ph, pl, gi, gv])
    
    ph_a = _air_eq_pressure(1,0,0,a,t,dh)
    ph_t = _air_eq_pressure(0,1,0,a,t,dh)
    ph_d = _air_eq_pressure(0,0,1,a,t,dh)
    pl_t = _flu_eq_pressure(1,0,t,dl)
    pl_d = _flu_eq_pressure(0,1,t,dl)
    gl_t = _flu_eq_chempot(1,0,t,dl)
    gl_d = _flu_eq_chempot(0,1,t,dl)
    gi_t = _ice_g(1,0,t,pres)
    gv_a = _air_eq_vappot(1,0,0,a,t,dh)
    gv_t = _air_eq_vappot(0,1,0,a,t,dh)
    gv_d = _air_eq_vappot(0,0,1,a,t,dh)
    dlhs = numpy.array([[0.,0.,0.,0.], [0.,0.,0.,0.], [0.,gl_t,0.,gl_d],
        [0.,gl_t,0.,gl_d]])
    drhs = numpy.array([[ph_a,ph_t,ph_d,0.], [0.,pl_t,0.,pl_d], [0.,gi_t,0.,0.],
        [gv_a,gv_t,gv_d,0.]])
    return lhs, rhs, dlhs, drhs

def _diff_wef(a,t,p,dh,dl,wair,entr,wetf):
    """Calculate wet-icy air disequilibrium at WEF.
    
    Calculate both sides of the equations
    
        given pressure = pressure in humid air
        given pressure = pressure of liquid water
        chemical potential of liquid water = potential of ice
        chemical potential of liquid water = potential of water vapour
        given entropy = entropy of wet-icy air
    
    and their Jacobians with respect to humid air dry fraction,
    temperature, pressure, humid air density, and liquid water density.
    Solving these equations gives equilibrium values at the given total
    dry fraction, specific entropy, and wet fraction of condensate.
    
    :arg float a: Humid air dry fraction in kg/kg.
    :arg float t: Temperature in K.
    :arg float p: Pressure in Pa.
    :arg float dh: Humid air density in kg/m3.
    :arg float dl: Liquid water density in kg/m3.
    :arg float wair: Total dry fraction in kg/kg.
    :arg float entr: Entropy in J/kg/K.
    :arg float wetf: Wet fraction of condensate in kg/kg.
    :returns: Left-hand side of the equation, right-hand side,
        Jacobian of LHS, and Jacobian of RHS.
    :rtype: tuple(array(float))
    """
    ph = _air_eq_pressure(0,0,0,a,t,dh)
    pl = _flu_eq_pressure(0,0,t,dl)
    gl = _flu_eq_chempot(0,0,t,dl)
    gi = _ice_g(0,0,t,p)
    gv = _air_eq_vappot(0,0,0,a,t,dh)
    sh = -_air_f(0,1,0,a,t,dh)
    sl = -_flu_f(1,0,t,dl)
    si = -_ice_g(1,0,t,p)
    s = wair/a*sh + wetf*(1-wair/a)*sl + (1-wetf)*(1-wair/a)*si
    lhs = numpy.array([p, p, gl, gl, entr])
    rhs = numpy.array([ph, pl, gi, gv, s])
    
    ph_a = _air_eq_pressure(1,0,0,a,t,dh)
    ph_t = _air_eq_pressure(0,1,0,a,t,dh)
    ph_d = _air_eq_pressure(0,0,1,a,t,dh)
    pl_t = _flu_eq_pressure(1,0,t,dl)
    pl_d = _flu_eq_pressure(0,1,t,dl)
    gl_t = _flu_eq_chempot(1,0,t,dl)
    gl_d = _flu_eq_chempot(0,1,t,dl)
    gi_t = _ice_g(1,0,t,p)
    gi_p = _ice_g(0,1,t,p)
    gv_a = _air_eq_vappot(1,0,0,a,t,dh)
    gv_t = _air_eq_vappot(0,1,0,a,t,dh)
    gv_d = _air_eq_vappot(0,0,1,a,t,dh)
    sh_a = -_air_f(1,1,0,a,t,dh)
    sh_t = -_air_f(0,2,0,a,t,dh)
    sh_d = -_air_f(0,1,1,a,t,dh)
    sl_t = -_flu_f(2,0,t,dl)
    sl_d = -_flu_f(1,1,t,dl)
    si_t = -_ice_g(2,0,t,p)
    si_p = -_ice_g(1,1,t,p)
    s_a = -wair/a**2*(sh - a*sh_a - wetf*sl - (1-wetf)*si)
    s_t = wair/a*sh_t + wetf*(1-wair/a)*sl_t + (1-wetf)*(1-wair/a)*si_t
    s_p = (1-wetf)*(1-wair/a)*si_p
    s_dh = wair/a*sh_d
    s_dl = wetf*(1-wair/a)*sl_d
    dlhs = numpy.array([[0.,0.,1.,0.,0.], [0.,0.,1.,0.,0.],
        [0.,gl_t,0.,0.,gl_d], [0.,gl_t,0.,0.,gl_d], [0.,0.,0.,0.,0.]])
    drhs = numpy.array([[ph_a,ph_t,0.,ph_d,0.], [0.,pl_t,0.,0.,pl_d],
        [0.,gi_t,gi_p,0.,0.], [gv_a,gv_t,0.,gv_d,0.], [s_a,s_t,s_p,s_dh,s_dl]])
    return lhs, rhs, dlhs, drhs

def eq_atp(airf=None,temp=None,pres=None,dhum=None,dliq=None,
    chkvals=False,chktol=_CHKTOL,airf0=None,temp0=None,pres0=None,
    dhum0=None,dliq0=None,chkbnd=False,mathargs=None):
    """Get primary wet-icy air variables at ATP.
    
    Get the values of all primary variables for wet-icy air at any of
    the humid air dry fraction, temperature, or pressure. At least one
    of these values must be provided.
    
    If the calculation has already been done, the results can be passed
    to avoid unnecessary repeat calculations. If enough values are
    passed, they will be checked for consistency if chkvals is True.
    
    :arg airf: Humid air dry fraction in kg/kg. If unknown, pass None
        (default) and it will be calculated.
    :type airf: float or None
    :arg temp: Temperature in K. If unknown, pass None (default) and it
        will be calculated.
    :type temp: float or None
    :arg pres: Pressure in Pa. If unknown, pass None (default) and it
        will be calculated.
    :type pres: float or None
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
    :arg airf0: Initial guess for the humid air dry fraction in kg/kg.
        If None (default) then the appropriate `_approx_*` is used.
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
    :returns: Humid air dry fraction, temperature, pressure, humid air
        density, and liquid water density (all in SI units).
    :raises ValueError: If no values are provided.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    """
    if all(val is None for val in (airf,temp,pres)):
        errmsg = 'Must provide one of airf, temp, or pres'
        raise ValueError(errmsg)
    if mathargs is None:
        mathargs = dict()
    
    if airf is not None:
        if any(val is None for val in (temp,pres,dhum,dliq)):
            x0 = (temp0,pres0,dhum0,dliq0)
            fargs = (airf,)
            x1 = _newton(_diff_a,x0,_approx_a,fargs=fargs,**mathargs)
            temp, pres, dhum, dliq = x1
    elif temp is not None:
        if any(val is None for val in (airf,pres,dhum,dliq)):
            x0 = (airf0,pres0,dhum0,dliq0)
            fargs = (temp,)
            x1 = _newton(_diff_t,x0,_approx_t,fargs=fargs,**mathargs)
            airf, pres, dhum, dliq = x1
    else:
        if any(val is None for val in (airf,temp,dhum,dliq)):
            x0 = (airf0,temp0,dhum0,dliq0)
            fargs = (pres,)
            x1 = _newton(_diff_p,x0,_approx_p,fargs=fargs,**mathargs)
            airf, temp, dhum, dliq = x1

    _chkflubnds(temp,dliq,chkbnd=chkbnd)
    _chkhumbnds(airf,temp,dhum,chkbnd=chkbnd)
    _chkicebnds(temp,pres,chkbnd=chkbnd)
    if not chkvals:
        return airf, temp, pres, dhum, dliq
    
    lhs, rhs, __, __ = _diff_a(temp,pres,dhum,dliq,airf)
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

def eq_wefli(wair,entr=None,wetf=None,wliq=None,wice=None,airf=None,
    temp=None,pres=None,dhum=None,dliq=None,chkvals=False,
    chktol=_CHKTOL,airf0=None,temp0=None,pres0=None,dhum0=None,
    dliq0=None,chkbnd=False,mathargs=None):
    """Get equilibrium values at WEF or WLI.
    
    Get the values of all primary variables for wet-icy air with the
    given properties. The properties can be either the total dry
    fraction, specific entropy, and wet fraction of condensate; or the
    total mass fractions of dry air, liquid water, and ice.
    
    If the calculation has already been done, the results can be passed
    to avoid unnecessary repeat calculations. If enough values are
    passed, they will be checked for consistency if chkvals is True.
    
    :arg float wair: Total dry fraction in kg/kg.
    :arg entr: Entropy in J/kg/K.
    :type entr: float or None
    :arg wetf: Wet fraction of condensate in kg/kg.
    :type wetf: float or None
    :arg wliq: Mass fraction of liquid water in kg/kg.
    :type wliq: float or None
    :arg wice: Mass fraction of ice in kg/kg.
    :type wice: float or None
    :arg airf: Humid air dry fraction in kg/kg. If unknown, pass None
        (default) and it will be calculated.
    :type airf: float or None
    :arg temp: Temperature in K. If unknown, pass None (default) and it
        will be calculated.
    :type temp: float or None
    :arg pres: Pressure in Pa. If unknown, pass None (default) and it
        will be calculated.
    :type pres: float or None
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
    :arg airf0: Initial guess for the humid air dry fraction in kg/kg.
        If None (default) then the appropriate `_approx_*` is used.
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
    :returns: Humid air dry fraction, temperature, pressure, humid air
        density, and liquid water density (all in SI units).
    :raises ValueError: If not enough values are provided.
    :raises ValueError: If all mass fractions are provided but their sum
        is larger than 1.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    """
    cond1 = (entr is None or wetf is None)
    cond2 = (wliq is None or wice is None)
    if cond1 and cond2:
        errmsg = ('Not enough values were provided. Please provide wair and '
            'either (entr,wetf) or (wliq,wice).')
        raise ValueError(errmsg)
    if mathargs is None:
        mathargs = dict()
    
    if not cond1:
        if any(val is None for val in (airf,temp,pres,dhum,dliq)):
            x0 = (airf0,temp0,pres0,dhum0,dliq0)
            fargs = (wair,entr,wetf)
            x1 = _newton(_diff_wef,x0,_approx_wef,fargs=fargs,**mathargs)
            airf, temp, pres, dhum, dliq = x1
    else:
        wtot = wair + wliq + wice
        if wtot >= 1:
            errmsg = ('The mass fractions {0} sum to more than '
                '1').format((wair,wliq,wice))
            raise ValueError(errmsg)
        if airf is None:
            airf = wair / (1-wliq-wice)
        if any(val is None for val in (temp,pres,dhum,dliq)):
            x0 = (temp0,pres0,dhum0,dliq0)
            fargs = (airf,)
            x1 = _newton(_diff_a,x0,_approx_a,fargs=fargs,**mathargs)
            temp, pres, dhum, dliq = x1
    
    _chkflubnds(temp,dliq,chkbnd=chkbnd)
    _chkhumbnds(airf,temp,dhum,chkbnd=chkbnd)
    _chkicebnds(temp,pres,chkbnd=chkbnd)
    if not chkvals:
        return airf, temp, pres, dhum, dliq
    
    if entr is None:
        sh = -_air_f(0,1,0,airf,temp,dhum)
        sl = -_flu_f(1,0,temp,dliq)
        si = -_ice_g(1,0,temp,pres)
        entr = wair/airf*sh + wliq*sl + wice*si
    if wetf is None:
        wetf = wliq/(wliq + wice)
    lhs, rhs, __, __ = _diff_wef(airf,temp,pres,dhum,dliq,wair,entr,wetf)
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

def eq_all(wair=None,entr=None,wetf=None,wliq=None,wice=None,airf=None,
    temp=None,pres=None,dhum=None,dliq=None,chkvals=False,
    chktol=_CHKTOL,airf0=None,temp0=None,pres0=None,dhum0=None,
    dliq0=None,chkbnd=False,mathargs=None):
    """Get equilibrium values at ATP, WEF, or WLI.
    
    Get the values of all primary variables for wet-icy air with the
    given properties. The properties can be: any of the humid air dry
    fraction, temperature or pressure; or the total dry fraction,
    specific entropy, and wet fraction of condensate; or the total mass
    fractions of dry air, liquid water, and ice. This function only
    serves as a common wrapper for :func:`eq_atp` and :func:`eq_wefli`,
    which handle these cases separately.
    
    If the calculation has already been done, the results can be passed
    to avoid unnecessary repeat calculations. If enough values are
    passed, they will be checked for consistency if chkvals is True.
    
    :arg wair: Total dry fraction in kg/kg.
    :type wair: float or None
    :arg entr: Entropy in J/kg/K.
    :type entr: float or None
    :arg wetf: Wet fraction of condensate in kg/kg.
    :type wetf: float or None
    :arg wliq: Mass fraction of liquid water in kg/kg.
    :type wliq: float or None
    :arg wice: Mass fraction of ice in kg/kg.
    :type wice: float or None
    :arg airf: Humid air dry fraction in kg/kg. If unknown, pass None
        (default) and it will be calculated.
    :type airf: float or None
    :arg temp: Temperature in K. If unknown, pass None (default) and it
        will be calculated.
    :type temp: float or None
    :arg pres: Pressure in Pa. If unknown, pass None (default) and it
        will be calculated.
    :type pres: float or None
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
    :arg airf0: Initial guess for the humid air dry fraction in kg/kg.
        If None (default) then the appropriate `_approx_*` is used.
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
    :returns: Humid air dry fraction, temperature, pressure, humid air
        density, and liquid water density (all in SI units).
    :raises ValueError: If not enough values are provided.
    :raises ValueError: If all mass fractions are provided but their sum
        is larger than 1.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    """
    if all(val is None for val in (airf,temp,pres)):
        res = eq_wefli(wair,entr=entr,wetf=wetf,wliq=wliq,wice=wice,airf0=airf0,
            temp0=temp0,pres0=pres0,dhum0=dhum0,dliq0=dliq0,chkbnd=chkbnd,
            mathargs=mathargs)
    else:
        res = eq_atp(airf=airf,temp=temp,pres=pres,dhum=dhum,dliq=dliq,
            chkvals=chkvals,chktol=chktol,airf0=airf0,temp0=temp0,pres0=pres0,
            dhum0=dhum0,dliq0=dliq0,chkbnd=chkbnd,mathargs=mathargs)
    return res


## Thermodynamic properties not needing mass fractions
def airfraction(wair=None,entr=None,wetf=None,wliq=None,wice=None,
    airf=None,temp=None,pres=None,dhum=None,dliq=None,chkvals=False,
    chktol=_CHKTOL,airf0=None,temp0=None,pres0=None,dhum0=None,
    dliq0=None,chkbnd=False,mathargs=None):
    """Calculate wet-icy air humid air dry fraction.
    
    Calculate the mass fraction of dry air in humid air for wet-icy air.
    
    :arg wair: Total dry fraction in kg/kg.
    :type wair: float or None
    :arg entr: Entropy in J/kg/K.
    :type entr: float or None
    :arg wetf: Wet fraction of condensate in kg/kg.
    :type wetf: float or None
    :arg wliq: Mass fraction of liquid water in kg/kg.
    :type wliq: float or None
    :arg wice: Mass fraction of ice in kg/kg.
    :type wice: float or None
    :arg airf: Humid air dry fraction in kg/kg. If unknown, pass None
        (default) and it will be calculated.
    :type airf: float or None
    :arg temp: Temperature in K. If unknown, pass None (default) and it
        will be calculated.
    :type temp: float or None
    :arg pres: Pressure in Pa. If unknown, pass None (default) and it
        will be calculated.
    :type pres: float or None
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
    :arg airf0: Initial guess for the humid air dry fraction in kg/kg.
        If None (default) then the appropriate `_approx_*` is used.
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
    :returns: Humid air dry fraction in kg/kg.
    :raises ValueError: If not enough values are provided.
    :raises ValueError: If all mass fractions are provided but their sum
        is larger than 1.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> airfraction(pres=1e4)
    0.961024307544
    >>> airfraction(temp=273.155)
    0.994366063923
    >>> airfraction(wair=.99,entr=0.,wetf=.5)
    0.996583352944
    """
    airf, temp, pres, dhum, dliq = eq_all(wair=wair,entr=entr,wetf=wetf,
        wliq=wliq,wice=wice,airf=airf,temp=temp,pres=pres,dhum=dhum,dliq=dliq,
        chkvals=chkvals,chktol=chktol,airf0=airf0,temp0=temp0,pres0=pres0,
        dhum0=dhum0,dliq0=dliq0,chkbnd=chkbnd,mathargs=mathargs)
    return airf

def pressure(wair=None,entr=None,wetf=None,wliq=None,wice=None,
    airf=None,temp=None,pres=None,dhum=None,dliq=None,chkvals=False,
    chktol=_CHKTOL,airf0=None,temp0=None,pres0=None,dhum0=None,
    dliq0=None,chkbnd=False,mathargs=None):
    """Calculate wet-icy air pressure.
    
    Calculate the pressure of wet-icy air.
    
    :arg wair: Total dry fraction in kg/kg.
    :type wair: float or None
    :arg entr: Entropy in J/kg/K.
    :type entr: float or None
    :arg wetf: Wet fraction of condensate in kg/kg.
    :type wetf: float or None
    :arg wliq: Mass fraction of liquid water in kg/kg.
    :type wliq: float or None
    :arg wice: Mass fraction of ice in kg/kg.
    :type wice: float or None
    :arg airf: Humid air dry fraction in kg/kg. If unknown, pass None
        (default) and it will be calculated.
    :type airf: float or None
    :arg temp: Temperature in K. If unknown, pass None (default) and it
        will be calculated.
    :type temp: float or None
    :arg pres: Pressure in Pa. If unknown, pass None (default) and it
        will be calculated.
    :type pres: float or None
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
    :arg airf0: Initial guess for the humid air dry fraction in kg/kg.
        If None (default) then the appropriate `_approx_*` is used.
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
    :raises ValueError: If all mass fractions are provided but their sum
        is larger than 1.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> pressure(airf=.99)
    38338.9622424
    >>> pressure(temp=273.155)
    67931.60108
    >>> pressure(wair=.99,entr=0.,wetf=.5)
    112016.075795
    >>> pressure(wair=.1,wliq=.2,wice=.3)
    706.817425301
    """
    airf, temp, pres, dhum, dliq = eq_all(wair=wair,entr=entr,wetf=wetf,
        wliq=wliq,wice=wice,airf=airf,temp=temp,pres=pres,dhum=dhum,dliq=dliq,
        chkvals=chkvals,chktol=chktol,airf0=airf0,temp0=temp0,pres0=pres0,
        dhum0=dhum0,dliq0=dliq0,chkbnd=chkbnd,mathargs=mathargs)
    return pres

def temperature(wair=None,entr=None,wetf=None,wliq=None,wice=None,
    airf=None,temp=None,pres=None,dhum=None,dliq=None,chkvals=False,
    chktol=_CHKTOL,airf0=None,temp0=None,pres0=None,dhum0=None,
    dliq0=None,chkbnd=False,mathargs=None):
    """Calculate wet-icy air temperature.
    
    Calculate the temperature of wet-icy air.
    
    :arg wair: Total dry fraction in kg/kg.
    :type wair: float or None
    :arg entr: Entropy in J/kg/K.
    :type entr: float or None
    :arg wetf: Wet fraction of condensate in kg/kg.
    :type wetf: float or None
    :arg wliq: Mass fraction of liquid water in kg/kg.
    :type wliq: float or None
    :arg wice: Mass fraction of ice in kg/kg.
    :type wice: float or None
    :arg airf: Humid air dry fraction in kg/kg. If unknown, pass None
        (default) and it will be calculated.
    :type airf: float or None
    :arg temp: Temperature in K. If unknown, pass None (default) and it
        will be calculated.
    :type temp: float or None
    :arg pres: Pressure in Pa. If unknown, pass None (default) and it
        will be calculated.
    :type pres: float or None
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
    :arg airf0: Initial guess for the humid air dry fraction in kg/kg.
        If None (default) then the appropriate `_approx_*` is used.
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
    :raises ValueError: If all mass fractions are provided but their sum
        is larger than 1.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> temperature(airf=.99)
    273.157198087
    >>> temperature(pres=1e4)
    273.159302793
    >>> temperature(wair=.99,entr=0.,wetf=.5)
    273.151724970
    >>> temperature(wair=.1,wliq=.2,wice=.3)
    273.159992933
    """
    airf, temp, pres, dhum, dliq = eq_all(wair=wair,entr=entr,wetf=wetf,
        wliq=wliq,wice=wice,airf=airf,temp=temp,pres=pres,dhum=dhum,dliq=dliq,
        chkvals=chkvals,chktol=chktol,airf0=airf0,temp0=temp0,pres0=pres0,
        dhum0=dhum0,dliq0=dliq0,chkbnd=chkbnd,mathargs=mathargs)
    return temp


## Equilibrium properties requiring total mass fractions
def density(wair,entr=None,wetf=None,wliq=None,wice=None,airf=None,
    temp=None,pres=None,dhum=None,dliq=None,chkvals=False,
    chktol=_CHKTOL,airf0=None,temp0=None,pres0=None,dhum0=None,
    dliq0=None,chkbnd=False,mathargs=None):
    """Calculate wet-icy air density.
    
    Calculate the density of wet-icy air.
    
    :arg float wair: Total dry fraction in kg/kg.
    :arg entr: Entropy in J/kg/K.
    :type entr: float or None
    :arg wetf: Wet fraction of condensate in kg/kg.
    :type wetf: float or None
    :arg wliq: Mass fraction of liquid water in kg/kg.
    :type wliq: float or None
    :arg wice: Mass fraction of ice in kg/kg.
    :type wice: float or None
    :arg airf: Humid air dry fraction in kg/kg. If unknown, pass None
        (default) and it will be calculated.
    :type airf: float or None
    :arg temp: Temperature in K. If unknown, pass None (default) and it
        will be calculated.
    :type temp: float or None
    :arg pres: Pressure in Pa. If unknown, pass None (default) and it
        will be calculated.
    :type pres: float or None
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
    :arg airf0: Initial guess for the humid air dry fraction in kg/kg.
        If None (default) then the appropriate `_approx_*` is used.
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
    :returns: Density in kg/m3.
    :raises ValueError: If not enough values are provided.
    :raises ValueError: If all mass fractions are provided but their sum
        is larger than 1.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> density(0.99,entr=0.,wetf=.5)
    1.43611528680
    >>> density(.1,wliq=.2,wice=.3)
    0.0121364037568
    """
    airf, temp, pres, dhum, dliq = eq_wefli(wair,entr=entr,wetf=wetf,wliq=wliq,
        wice=wice,airf=airf,temp=temp,pres=pres,dhum=dhum,dliq=dliq,
        chkvals=chkvals,chktol=chktol,airf0=airf0,temp0=temp0,pres0=pres0,
        dhum0=dhum0,dliq0=dliq0,chkbnd=chkbnd,mathargs=mathargs)
    if wliq is None or wice is None:
        sh = -_air_f(0,1,0,airf,temp,dhum)
        si = -_ice_g(1,0,temp,pres)
        if wetf == 0:
            wliq = 0.
            wice = (entr - wair/airf*sh)/si
        else:
            sl = -_flu_f(1,0,temp,dliq)
            wliq = (entr - wair/airf*sh)/(sl + (1-wetf)/wetf*si)
            wice = wliq * (1-wetf)/wetf
    
    whum = wair/airf
    vi = _ice_g(0,1,temp,pres,chkbnd=chkbnd)
    vtot = whum/dhum + wliq/dliq + wice*vi
    rho = vtot**(-1)
    return rho

def dryairfraction(wair,entr=None,wetf=None,wliq=None,wice=None,
    airf=None,temp=None,pres=None,dhum=None,dliq=None,chkvals=False,
    chktol=_CHKTOL,airf0=None,temp0=None,pres0=None,dhum0=None,
    dliq0=None,chkbnd=False,mathargs=None):
    """Calculate wet-icy air total dry fraction.
    
    Calculate the total mass fraction of dry air in wet-icy air.
    
    :arg float wair: Total dry fraction in kg/kg.
    :arg entr: Entropy in J/kg/K.
    :type entr: float or None
    :arg wetf: Wet fraction of condensate in kg/kg.
    :type wetf: float or None
    :arg wliq: Mass fraction of liquid water in kg/kg.
    :type wliq: float or None
    :arg wice: Mass fraction of ice in kg/kg.
    :type wice: float or None
    :arg airf: Humid air dry fraction in kg/kg. If unknown, pass None
        (default) and it will be calculated.
    :type airf: float or None
    :arg temp: Temperature in K. If unknown, pass None (default) and it
        will be calculated.
    :type temp: float or None
    :arg pres: Pressure in Pa. If unknown, pass None (default) and it
        will be calculated.
    :type pres: float or None
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
    :arg airf0: Initial guess for the humid air dry fraction in kg/kg.
        If None (default) then the appropriate `_approx_*` is used.
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
    :returns: Total dry fraction in kg/kg.
    :raises ValueError: If not enough values are provided.
    :raises ValueError: If all mass fractions are provided but their sum
        is larger than 1.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    """
    return wair

def enthalpy(wair,entr=None,wetf=None,wliq=None,wice=None,airf=None,
    temp=None,pres=None,dhum=None,dliq=None,chkvals=False,
    chktol=_CHKTOL,airf0=None,temp0=None,pres0=None,dhum0=None,
    dliq0=None,chkbnd=False,mathargs=None):
    """Calculate wet-icy air enthalpy.
    
    Calculate the specific enthalpy of wet-icy air.
    
    :arg float wair: Total dry fraction in kg/kg.
    :arg entr: Entropy in J/kg/K.
    :type entr: float or None
    :arg wetf: Wet fraction of condensate in kg/kg.
    :type wetf: float or None
    :arg wliq: Mass fraction of liquid water in kg/kg.
    :type wliq: float or None
    :arg wice: Mass fraction of ice in kg/kg.
    :type wice: float or None
    :arg airf: Humid air dry fraction in kg/kg. If unknown, pass None
        (default) and it will be calculated.
    :type airf: float or None
    :arg temp: Temperature in K. If unknown, pass None (default) and it
        will be calculated.
    :type temp: float or None
    :arg pres: Pressure in Pa. If unknown, pass None (default) and it
        will be calculated.
    :type pres: float or None
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
    :arg airf0: Initial guess for the humid air dry fraction in kg/kg.
        If None (default) then the appropriate `_approx_*` is used.
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
    :raises ValueError: If all mass fractions are provided but their sum
        is larger than 1.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> enthalpy(.99,entr=0.,wetf=.5)
    7356.12943724
    >>> enthalpy(.1,wliq=.2,wice=.3)
    900361.135280
    """
    airf, temp, pres, dhum, dliq = eq_wefli(wair,entr=entr,wetf=wetf,wliq=wliq,
        wice=wice,airf=airf,temp=temp,pres=pres,dhum=dhum,dliq=dliq,
        chkvals=chkvals,chktol=chktol,airf0=airf0,temp0=temp0,pres0=pres0,
        dhum0=dhum0,dliq0=dliq0,chkbnd=chkbnd,mathargs=mathargs)
    if wliq is None or wice is None:
        sh = -_air_f(0,1,0,airf,temp,dhum)
        si = -_ice_g(1,0,temp,pres)
        if wetf == 0:
            wliq = 0.
            wice = (entr - wair/airf*sh)/si
        else:
            sl = -_flu_f(1,0,temp,dliq)
            wliq = (entr - wair/airf*sh)/(sl + (1-wetf)/wetf*si)
            wice = wliq * (1-wetf)/wetf
    
    whum = wair/airf
    hh = air2.enthalpy(airf,temp,dhum)
    hl = flu2.enthalpy(temp,dliq)
    hi = ice2.enthalpy(temp,pres)
    h = whum*hh + wliq*hl + wice*hi
    return h

def entropy(wair,entr=None,wetf=None,wliq=None,wice=None,airf=None,
    temp=None,pres=None,dhum=None,dliq=None,chkvals=False,
    chktol=_CHKTOL,airf0=None,temp0=None,pres0=None,dhum0=None,
    dliq0=None,chkbnd=False,mathargs=None):
    """Calculate wet-icy air entropy.
    
    Calculate the specific entropy of wet-icy air.
    
    :arg float wair: Total dry fraction in kg/kg.
    :arg entr: Entropy in J/kg/K.
    :type entr: float or None
    :arg wetf: Wet fraction of condensate in kg/kg.
    :type wetf: float or None
    :arg wliq: Mass fraction of liquid water in kg/kg.
    :type wliq: float or None
    :arg wice: Mass fraction of ice in kg/kg.
    :type wice: float or None
    :arg airf: Humid air dry fraction in kg/kg. If unknown, pass None
        (default) and it will be calculated.
    :type airf: float or None
    :arg temp: Temperature in K. If unknown, pass None (default) and it
        will be calculated.
    :type temp: float or None
    :arg pres: Pressure in Pa. If unknown, pass None (default) and it
        will be calculated.
    :type pres: float or None
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
    :arg airf0: Initial guess for the humid air dry fraction in kg/kg.
        If None (default) then the appropriate `_approx_*` is used.
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
    :raises ValueError: If all mass fractions are provided but their sum
        is larger than 1.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> entropy(.1,wliq=.2,wice=.3)
    3496.16306903
    """
    if entr is not None:
        return entr
    airf, temp, pres, dhum, dliq = eq_wefli(wair,wliq=wliq,wice=wice,airf=airf,
        temp=temp,pres=pres,dhum=dhum,dliq=dliq,chkvals=chkvals,chktol=chktol,
        airf0=airf0,temp0=temp0,pres0=pres0,dhum0=dhum0,dliq0=dliq0,
        chkbnd=chkbnd,mathargs=mathargs)
    
    whum = wair/airf
    sh = -_air_f(0,1,0,airf,temp,dhum)
    sl = -_flu_f(1,0,temp,dliq)
    si = -_ice_g(1,0,temp,pres)
    s = whum*sh + wliq*sl + wice*si
    return s

def liquidfraction(wair,entr=None,wetf=None,wliq=None,wice=None,
    airf=None,temp=None,pres=None,dhum=None,dliq=None,chkvals=False,
    chktol=_CHKTOL,airf0=None,temp0=None,pres0=None,dhum0=None,
    dliq0=None,chkbnd=False,mathargs=None):
    """Calculate wet-icy air liquid water fraction.
    
    Calculate the mass fraction of liquid water in wet-icy air.
    
    :arg float wair: Total dry fraction in kg/kg.
    :arg entr: Entropy in J/kg/K.
    :type entr: float or None
    :arg wetf: Wet fraction of condensate in kg/kg.
    :type wetf: float or None
    :arg wliq: Mass fraction of liquid water in kg/kg.
    :type wliq: float or None
    :arg wice: Mass fraction of ice in kg/kg.
    :type wice: float or None
    :arg airf: Humid air dry fraction in kg/kg. If unknown, pass None
        (default) and it will be calculated.
    :type airf: float or None
    :arg temp: Temperature in K. If unknown, pass None (default) and it
        will be calculated.
    :type temp: float or None
    :arg pres: Pressure in Pa. If unknown, pass None (default) and it
        will be calculated.
    :type pres: float or None
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
    :arg airf0: Initial guess for the humid air dry fraction in kg/kg.
        If None (default) then the appropriate `_approx_*` is used.
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
    :returns: Liquid water fraction in kg/kg.
    :raises ValueError: If not enough values are provided.
    :raises ValueError: If all mass fractions are provided but their sum
        is larger than 1.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> liquidfraction(.99,entr=0.,wetf=.5)
    3.30296152581e-3
    """
    if wliq is not None:
        return wliq
    airf, temp, pres, dhum, dliq = eq_wefli(wair,entr=entr,wetf=wetf,wliq=wliq,
        wice=wice,airf=airf,temp=temp,pres=pres,dhum=dhum,dliq=dliq,
        chkvals=chkvals,chktol=chktol,airf0=airf0,temp0=temp0,pres0=pres0,
        dhum0=dhum0,dliq0=dliq0,chkbnd=chkbnd,mathargs=mathargs)
    
    if wetf == 0:
        wliq = 0.
    else:
        sh = -_air_f(0,1,0,airf,temp,dhum)
        si = -_ice_g(1,0,temp,pres)
        sl = -_flu_f(1,0,temp,dliq)
        wliq = (entr - wair/airf*sh)/(sl + (1-wetf)/wetf*si)
    return wliq

def solidfraction(wair,entr=None,wetf=None,wliq=None,wice=None,
    airf=None,temp=None,pres=None,dhum=None,dliq=None,chkvals=False,
    chktol=_CHKTOL,airf0=None,temp0=None,pres0=None,dhum0=None,
    dliq0=None,chkbnd=False,mathargs=None):
    """Calculate wet-icy air ice mass fraction.
    
    Calculate the mass fraction of ice in wet-icy air.
    
    :arg float wair: Total dry fraction in kg/kg.
    :arg entr: Entropy in J/kg/K.
    :type entr: float or None
    :arg wetf: Wet fraction of condensate in kg/kg.
    :type wetf: float or None
    :arg wliq: Mass fraction of liquid water in kg/kg.
    :type wliq: float or None
    :arg wice: Mass fraction of ice in kg/kg.
    :type wice: float or None
    :arg airf: Humid air dry fraction in kg/kg. If unknown, pass None
        (default) and it will be calculated.
    :type airf: float or None
    :arg temp: Temperature in K. If unknown, pass None (default) and it
        will be calculated.
    :type temp: float or None
    :arg pres: Pressure in Pa. If unknown, pass None (default) and it
        will be calculated.
    :type pres: float or None
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
    :arg airf0: Initial guess for the humid air dry fraction in kg/kg.
        If None (default) then the appropriate `_approx_*` is used.
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
    :returns: Ice fraction in kg/kg.
    :raises ValueError: If not enough values are provided.
    :raises ValueError: If all mass fractions are provided but their sum
        is larger than 1.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> solidfraction(.99,entr=0.,wetf=.5)
    3.30296152581e-3
    """
    if wice is not None:
        return wice
    airf, temp, pres, dhum, dliq = eq_wefli(wair,entr=entr,wetf=wetf,wliq=wliq,
        wice=wice,airf=airf,temp=temp,pres=pres,dhum=dhum,dliq=dliq,
        chkvals=chkvals,chktol=chktol,airf0=airf0,temp0=temp0,pres0=pres0,
        dhum0=dhum0,dliq0=dliq0,chkbnd=chkbnd,mathargs=mathargs)
    
    sh = -_air_f(0,1,0,airf,temp,dhum)
    si = -_ice_g(1,0,temp,pres)
    if wetf == 0:
        wliq = 0.
        wice = (entr - wair/airf*sh)/si
    else:
        sl = -_flu_f(1,0,temp,dliq)
        wliq = (entr - wair/airf*sh)/(sl + (1-wetf)/wetf*si)
        wice = wliq * (1-wetf)/wetf
    return wice

def vapourfraction(wair,entr=None,wetf=None,wliq=None,wice=None,
    airf=None,temp=None,pres=None,dhum=None,dliq=None,chkvals=False,
    chktol=_CHKTOL,airf0=None,temp0=None,pres0=None,dhum0=None,
    dliq0=None,chkbnd=False,mathargs=None):
    """Calculate wet-icy air water vapour fraction.
    
    Calculate the mass fraction of water vapour in wet-icy air.
    
    :arg float wair: Total dry fraction in kg/kg.
    :arg entr: Entropy in J/kg/K.
    :type entr: float or None
    :arg wetf: Wet fraction of condensate in kg/kg.
    :type wetf: float or None
    :arg wliq: Mass fraction of liquid water in kg/kg.
    :type wliq: float or None
    :arg wice: Mass fraction of ice in kg/kg.
    :type wice: float or None
    :arg airf: Humid air dry fraction in kg/kg. If unknown, pass None
        (default) and it will be calculated.
    :type airf: float or None
    :arg temp: Temperature in K. If unknown, pass None (default) and it
        will be calculated.
    :type temp: float or None
    :arg pres: Pressure in Pa. If unknown, pass None (default) and it
        will be calculated.
    :type pres: float or None
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
    :arg airf0: Initial guess for the humid air dry fraction in kg/kg.
        If None (default) then the appropriate `_approx_*` is used.
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
    :returns: Water vapour fraction in kg/kg.
    :raises ValueError: If not enough values are provided.
    :raises ValueError: If all mass fractions are provided but their sum
        is larger than 1.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> vapourfraction(.99,entr=0.,wetf=.5)
    3.39407694837e-3
    """
    if wliq is not None and wice is not None:
        wvap = 1 - wair - wliq - wice
        return wvap
    
    airf, temp, pres, dhum, dliq = eq_wefli(wair,entr=entr,wetf=wetf,wliq=wliq,
        wice=wice,airf=airf,temp=temp,pres=pres,dhum=dhum,dliq=dliq,
        chkvals=chkvals,chktol=chktol,airf0=airf0,temp0=temp0,pres0=pres0,
        dhum0=dhum0,dliq0=dliq0,chkbnd=chkbnd,mathargs=mathargs)
    sh = -_air_f(0,1,0,airf,temp,dhum)
    si = -_ice_g(1,0,temp,pres)
    if wetf == 0:
        wliq = 0.
        wice = (entr - wair/airf*sh)/si
    else:
        sl = -_flu_f(1,0,temp,dliq)
        wliq = (entr - wair/airf*sh)/(sl + (1-wetf)/wetf*si)
        wice = wliq * (1-wetf)/wetf
    wvap = 1 - wair - wliq - wice
    return wvap


## Isentropic melting and freezing levels
def iml(wair,entr,airf0=None,temp0=None,pres0=None,dhum0=None,
    dliq0=None,chkbnd=False,mathargs=None):
    """Calculate isentropic melting level.
    
    Calculate the isentropic melting level of wet-icy air, the pressure
    below which all condensed water will be ice at equilibrium.
    
    :arg float wair: Total dry fraction in kg/kg.
    :arg float entr: Entropy in J/kg/K.
    :arg airf0: Initial guess for the humid air dry fraction in kg/kg.
        If None (default) then the appropriate `_approx_*` is used.
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
    :returns: Pressure at the isentropic melting level in Pa.
    
    :Examples:
    
    >>> iml(.99,100.)
    81605.5557729
    """
    wetf = 0.
    airf, temp, pres, dhum, dliq = eq_wefli(wair,entr=entr,wetf=wetf,
        airf0=airf0,temp0=temp0,pres0=pres0,dhum0=dhum0,dliq0=dliq0,
        chkbnd=chkbnd,mathargs=mathargs)
    return pres

def ifl(wair,entr,airf0=None,temp0=None,pres0=None,dhum0=None,
    dliq0=None,chkbnd=False,mathargs=None):
    """Calculate isentropic freezing level.
    
    Calculate the isentropic freezing level of wet-icy air, the pressure
    above which all condensed water will be liquid water at equilibrium.
    
    :arg float wair: Total dry fraction in kg/kg.
    :arg float entr: Entropy in J/kg/K.
    :arg airf0: Initial guess for the humid air dry fraction in kg/kg.
        If None (default) then the appropriate `_approx_*` is used.
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
    :returns: Pressure at the isentropic freezing level in Pa.
    
    :Examples:
    
    >>> ifl(.99,100.)
    83234.7314358
    """
    wetf = 1.
    airf, temp, pres, dhum, dliq = eq_wefli(wair,entr=entr,wetf=wetf,
        airf0=airf0,temp0=temp0,pres0=pres0,dhum0=dhum0,dliq0=dliq0,
        chkbnd=chkbnd,mathargs=mathargs)
    return pres

