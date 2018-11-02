"""Icy air enthalpy and potential temperature.

This module provides the specific enthalpy of ice-saturated (icy) air.
It also provides the potential temperature and other adiabatic
properties. The primary variables are the total dry air mass fraction,
entropy, and pressure. As in :mod:`iceair4b`, the total dry fraction
here is the ratio of dry air mass to total parcel mass, including ice.

:Examples:

>>> iceair_h(0,0,0,0.5,1e5,entr=-600.)
-164588.1060
>>> iceair_h(0,1,0,0.5,1e5,entr=-600.)
271.4499944
>>> iceair_h(0,1,1,0.5,1e5,entr=-600.)
2.268401085e-04
>>> temperature(0.9,1e5,entr=-100.)
270.3836801
>>> cp(0.9,1e5,entr=-100.)
1768.514397
>>> lapserate(0.9,1e5,entr=-100.)
4.420925682e-04
>>> potdensity(0.9,230.,5e4,1e5)
1.450481104
>>> potenthalpy(0.9,230.,5e4,1e5)
-35781.25645
>>> pottemp(0.9,230.,5e4,1e5)
266.1052089

:Functions:

* :func:`eq_wpte`: Calculate icy air equilibrium properties at total dry
  fraction, pressure, and either entropy or temperature.
* :func:`iceair_h`: Icy air specific enthalpy with derivatives.
* :func:`cp`: Icy air isobaric heat capacity.
* :func:`density`: Icy air density.
* :func:`kappa_s`: Icy air isentropic compressibility.
* :func:`lapserate`: Icy air adiabatic lapse rate.
* :func:`temperature`: Icy air temperature.
* :func:`eq_pot`: Calculate icy air equilibrium properties under
  adiabatic displacement.
* :func:`potdensity`: Icy air potential density.
* :func:`potenthalpy`: Icy air potential enthalpy.
* :func:`pottemp`: Icy air potential temperature.

"""

__all__ = ['eq_wpte','liqair_h',
    'cp','density','kappa_s','lapserate','temperature',
    'eq_pot','potdensity','potenthalpy','pottemp']

import numpy
import warnings
import constants0
import ice1
import air2
import ice2
import maths3
import air3c
import maths4
import iceair4a
import iceair4b

_CHKTOL = constants0.CHKTOL
_RWAT = constants0.RWAT
_RDRY = constants0.RDRY
_PATM = constants0.PATM
_TCELS = constants0.TCELS
_TTP = constants0.TTP
_PTPE = constants0.PTPE
_LLVTP = constants0.LLVTP
_LILTP = constants0.LILTP
_CICE = constants0.CICE
_CVAP = constants0.CVAP
_CDRY = constants0.CDRY
_EPSW = _RDRY / _RWAT
_AVI = (_LLVTP+_LILTP) / (_RWAT*_TTP)
_BVI = (_CICE-_CVAP) / _RWAT
_chkhumbnds = constants0.chkhumbnds
_chkicebnds = constants0.chkicebnds
_ice_g = ice1.ice_g
_air_f = air2.air_f
_eq_pressure = air2.eq_pressure
_eq_vappot = air2.eq_vappot
_newton = maths3.newton


## Equilibrium functions
def _approx_wep(wair,entr,pres):
    """Approximate ATDh at WEP.
    
    Approximate the humid air dry fraction, temperature, and humid air
    density at the given total dry fraction, entropy, and pressure.
    
    :arg float wair: Total dry fraction in kg/kg.
    :arg float entr: Entropy in J/kg/K.
    :arg float pres: Pressure in Pa.
    :returns: Humid air dry fraction, temperature, and humid air density
        (all in SI units).
    """
    pvmax = pres * (1-wair) / (1-wair + _EPSW*wair)
    if pvmax >= _PTPE:
        # Saturation would start at _TTP; use saturated heat capacity at _TTP
        a_t = (pres - _PTPE)/(pres - _PTPE + _EPSW*_PTPE)
        s_t = (wair*_CDRY*numpy.log(_TTP/_TCELS) - (1-wair)*_LILTP/_TTP
            - wair*_RDRY*numpy.log((pres-_PTPE)/_PATM)
            + wair*_RWAT*_EPSW*_PTPE/(pres-_PTPE)*_AVI)
        c_t = (wair*_CDRY + wair*(1-a_t)/a_t*_CVAP + (1-wair/a_t)*_CICE
            + wair*_RWAT*(1-a_t)*(_EPSW*a_t + 1-a_t)/_EPSW/a_t**2 * _AVI**2)
        temp = _TTP * numpy.exp(-(s_t-entr)/c_t)
    else:
        # Get approximate saturation temperature
        v = numpy.log(pres*(1-wair)/(_PTPE*(_EPSW*wair + 1-wair)))/_BVI
        r = _AVI/_BVI
        x = maths4.lamb2(v,r)
        tsat = _TTP/x
        ssat = (wair * (_CDRY*numpy.log(tsat/_TCELS)
                - _RDRY*numpy.log((pres-pvmax)/_PATM))
            + (1-wair) * (_CVAP*numpy.log(tsat/_TTP) + _LLVTP/_TTP
                - _RWAT*numpy.log(pvmax/_PTPE)))
        
        if entr >= ssat:
            ceff = wair*_CDRY + (1-wair)*_CVAP
            temp = _TTP * numpy.exp((entr-ssat)/ceff)
        else:
            csat = (wair*_CDRY + (1-wair)*_CVAP
                + (1-wair)*_RWAT*pres/(pres-pvmax)
                    * ((_AVI+_BVI)*_TTP/tsat - _BVI)**2)
            temp = tsat * numpy.exp(-(ssat-entr)/csat)
    pvap = _PTPE * numpy.exp((_AVI+_BVI)*(1 - _TTP/temp)
        - _BVI*numpy.log(temp/_TTP))
    airf = (pres - pvap) / (pres - pvap + _EPSW*pvap)
    dhum = pres/(_RDRY*temp) / (airf + (1-airf)/_EPSW)
    return airf, temp, dhum

def _diff_wep(a,t,dh,wair,entr,pres):
    """Calculate icy air disequilibrium at WEP.
    
    Calculate both sides of the equations
    
        given pressure = pressure of humid air
        chemical potential of ice = potential of water vapour
        given entropy = entropy of icy air
    
    and their Jacobians with respect to humid air dry fraction,
    temperature, and humid air density. Solving these equations produces
    equilibrium values at the given total dry fraction, entropy, and
    pressure.
    
    :arg float a: Humid air dry fraction in kg/kg.
    :arg float t: Temperature in K.
    :arg float dh: Humid air density in kg/m3.
    :arg float wair: Total dry fraction in kg/kg.
    :arg float entr: Entropy in J/kg/K.
    :arg float pres: Pressure in Pa.
    :returns: Left-hand side of the equation, right-hand side,
        Jacobian of LHS, and Jacobian of RHS.
    :rtype: tuple(array(float))
    """
    ph = _eq_pressure(0,0,0,a,t,dh)
    gi = _ice_g(0,0,t,pres)
    gv = _eq_vappot(0,0,0,a,t,dh)
    sh = -_air_f(0,1,0,a,t,dh)
    si = -_ice_g(1,0,t,pres)
    stot = wair/a*sh + (1-wair/a)*si
    lhs = numpy.array([pres, gi, entr])
    rhs = numpy.array([ph, gv, stot])
    
    ph_a = _eq_pressure(1,0,0,a,t,dh)
    ph_t = _eq_pressure(0,1,0,a,t,dh)
    ph_d = _eq_pressure(0,0,1,a,t,dh)
    gi_t = _ice_g(1,0,t,pres)
    gv_a = _eq_vappot(1,0,0,a,t,dh)
    gv_t = _eq_vappot(0,1,0,a,t,dh)
    gv_d = _eq_vappot(0,0,1,a,t,dh)
    sh_a = -_air_f(1,1,0,a,t,dh)
    sh_t = -_air_f(0,2,0,a,t,dh)
    sh_d = -_air_f(0,1,1,a,t,dh)
    si_t = -_ice_g(2,0,t,pres)
    s_a = -wair/a**2 * (sh - a*sh_a - si)
    s_t = wair/a*sh_t + (1-wair/a)*si_t
    s_d = wair/a*sh_d
    dlhs = numpy.array([[0.,0.,0.], [0.,gi_t,0.], [0.,0.,0.]])
    drhs = numpy.array([[ph_a,ph_t,ph_d], [gv_a,gv_t,gv_d], [s_a,s_t,s_d]])
    return lhs, rhs, dlhs, drhs

def eq_wpte(wair,pres,entr=None,airf=None,temp=None,dhum=None,
    chkvals=False,chktol=_CHKTOL,airf0=None,temp0=None,dhum0=None,
    chkbnd=False,mathargs=None):
    """Get primary icy air variables at WP and T/E.
    
    Get the values of all primary variables for icy air at the given
    total dry fraction, pressure, and either entropy or temperature.
    
    If the calculation has already been done, the results can be passed
    to avoid unnecessary repeat calculations. If enough values are
    passed, they will be checked for consistency if chkvals is True.
    
    :arg float wair: Total dry fraction in kg/kg.
    :arg float pres: Pressure in Pa.
    :arg entr: Entropy in J/kg/K. If None (default) then `temp` must be
        provided.
    :type entr: float or None
    :arg airf: Humid air dry fraction in kg/kg. If unknown, pass None
        (default) and it will be calculated.
    :type airf: float or None
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
    :arg airf0: Initial guess for the humid air dry fraction in kg/kg.
        If None (default) then `_approx_wep` is used.
    :type airf0: float or None
    :arg temp0: Initial guess for the temperature in K. If None
        (default) then `_approx_wep` is used.
    :type temp0: float or None
    :arg dhum0: Initial guess for the humid air density in kg/m3. If
        None (default) then `_approx_wep` is used.
    :type dhum0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Humid air dry fraction, temperature, and humid air density
        (all in SI units).
    :raises ValueError: If both entr and temp are None.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    """
    if entr is None and temp is None:
        errmsg = 'Must provide at least one of entr or temp'
        raise ValueError(errmsg)
    if temp is None:
        x0 = (airf0,temp0,dhum0)
        fargs = (wair,entr,pres)
        if mathargs is None:
            mathargs = dict()
        x1 = _newton(_diff_wep,x0,_approx_wep,fargs=fargs,**mathargs)
        airf, temp, dhum = x1
    elif any(val is None for val in (airf,dhum)):
        airf, __, __, dhum = iceair4a.eq_atpe(temp=temp,pres=pres,
            chkvals=chkvals,chktol=chktol,airf0=airf0,dhum0=dhum0,
            mathargs=mathargs)
    
    _chkhumbnds(airf,temp,dhum,chkbnd=chkbnd)
    _chkicebnds(temp,pres,chkbnd=chkbnd)
    if not chkvals:
        return airf, temp, dhum
    
    if entr is None:
        entr = iceair4b.entropy(wair,temp,pres,airf=airf,dhum=dhum)
    lhs, rhs, __, __ = _diff_wep(airf,temp,dhum,wair,entr,pres)
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
    return airf, temp, dhum


## Enthalpy
def iceair_h(drvw,drve,drvp,wair,pres,entr=None,temp=None,airf=None,
    dhum=None,chkvals=False,chktol=_CHKTOL,airf0=None,temp0=None,
    dhum0=None,chkbnd=False,mathargs=None):
    """Calculate icy air enthalpy with derivatives.
    
    Calculate the specific enthalpy of icy air or its derivatives with
    respect to total dry fraction, entropy, and pressure.
    
    :arg int drvw: Number of total dry fraction derivatives.
    :arg int drve: Number of entropy derivatives.
    :arg int drvp: Number of pressure derivatives.
    :arg float wair: Total dry fraction in kg/kg.
    :arg float pres: Pressure in Pa.
    :arg entr: Entropy in J/kg/K. If None (default) then `temp` must be
        provided.
    :type entr: float or None
    :arg airf: Humid air dry fraction in kg/kg. If unknown, pass None
        (default) and it will be calculated.
    :type airf: float or None
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
    :arg airf0: Initial guess for the humid air dry fraction in kg/kg.
        If None (default) then `_approx_wep` is used.
    :type airf0: float or None
    :arg temp0: Initial guess for the temperature in K. If None
        (default) then `_approx_wep` is used.
    :type temp0: float or None
    :arg dhum0: Initial guess for the humid air density in kg/m3. If
        None (default) then `_approx_wep` is used.
    :type dhum0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Enthalpy in units of
        (J/kg) / (kg/kg)^drvw / (J/kg/K)^drve / Pa^drvp.
    :raises ValueError: If both entr and temp are None.
    :raises RuntimeWarning: If air with the given parameters would be
        unsaturated.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    :raises ValueError: If any of (drvw,drve,drvp) are negative, or if
        (drvw+drve+drvp)>2.
    
    :Examples:
    
    >>> iceair_h(0,0,0,0.5,1e5,entr=-600.)
    -164588.1060
    >>> iceair_h(1,0,0,0.5,1e5,entr=-600.)
    543.0166476
    >>> iceair_h(0,1,0,0.5,1e5,entr=-600.)
    271.4499944
    >>> iceair_h(0,0,1,0.5,1e5,entr=-600.)
    0.3919818785
    >>> iceair_h(2,0,0,0.5,1e5,entr=-600.)
    224806.0619
    >>> iceair_h(1,1,0,0.5,1e5,entr=-600.)
    -177.3368083
    >>> iceair_h(1,0,1,0.5,1e5,entr=-600.)
    0.4942223289
    >>> iceair_h(0,2,0,0.5,1e5,entr=-600.)
    0.1398909945
    >>> iceair_h(0,1,1,0.5,1e5,entr=-600.)
    2.268401085e-04
    >>> iceair_h(0,0,2,0.5,1e5,entr=-600.)
    -3.569822871e-06
    """
    if any(drv < 0 for drv in (drvw,drve,drvp)) or (drvw+drve+drvp) > 2:
        errmsg = 'Derivatives {0} not recognized'.format((drvw,drve,drvp))
        raise ValueError(errmsg)
    airf, temp, dhum = eq_wpte(wair,pres,entr=entr,temp=temp,airf=airf,
        dhum=dhum,chkvals=chkvals,chktol=chktol,airf0=airf0,temp0=temp0,
        dhum0=dhum0,chkbnd=chkbnd,mathargs=mathargs)
    if airf <= wair:
        warnmsg = 'Air with the given parameters is unsaturated'
        warnings.warn(warnmsg,RuntimeWarning)
        if entr is None:
            dhum = air3a.eq_atp(airf,temp,pres,dhum=dhum,chkvals=chkvals,
                chktol=chktol,dhum0=dhum0,chkbnd=chkbnd,mathargs=mathargs)
            entr = air3b.entropy(wair,temp,pres,dhum=dhum)
        elif temp is None:
            temp, dhum = air3c.eq_aep(airf,entr,pres,temp=temp,dhum=dhum,
                chkvals=chkvals,chktol=chktol,temp0=temp0,dhum0=dhum0,
                chkbnd=chkbnd,mathargs=mathargs)
        h = air3c.air_h(drvw,drvt,drvp,wair,entr,pres,temp=temp,dhum=dhum)
        return h
    
    w = wair / airf
    if (drvw,drve,drvp) == (0,0,0):
        hh = air2.enthalpy(airf,temp,dhum)
        hi = ice2.enthalpy(temp,pres)
        h = w*hh + (1-w)*hi
        return h
    elif (drvw,drve,drvp) == (1,0,0):
        fh_a = _air_f(1,0,0,airf,temp,dhum)
        h_w = fh_a
        return h_w
    elif (drvw,drve,drvp) == (0,1,0):
        h_e = temp
        return h_e
    elif (drvw,drve,drvp) == (0,0,1):
        gi_p = _ice_g(0,1,temp,pres)
        h_p = w/dhum + (1-w)*gi_p
        return h_p
    
    # Higher-order derivatives require inversion
    if entr is None:
        entr = iceair4b.entropy(wair,temp,pres,airf=airf,dhum=dhum)
    __, __, dlhs, drhs = _diff_wep(airf,temp,dhum,wair,entr,pres)
    jac = drhs - dlhs
    if (drvw,drve,drvp) == (2,0,0):
        fh_t = _air_f(0,1,0,airf,temp,dhum)
        fh_aa = _air_f(2,0,0,airf,temp,dhum)
        fh_at = _air_f(1,1,0,airf,temp,dhum)
        fh_ad = _air_f(1,0,1,airf,temp,dhum)
        gi_t = _ice_g(1,0,temp,pres)
        rhs_w = numpy.array([0.,0.,(-fh_t + gi_t)/airf])
        x_w = numpy.linalg.solve(jac,-rhs_w)
        h_wx = numpy.array([fh_aa,fh_at,fh_ad])
        h_ww = h_wx.dot(x_w)
        return h_ww
    elif (drvw,drve,drvp) == (1,1,0):
        fh_t = _air_f(0,1,0,airf,temp,dhum)
        gi_t = _ice_g(1,0,temp,pres)
        rhs_w = numpy.array([0.,0.,(-fh_t + gi_t)/airf])
        x_w = numpy.linalg.solve(jac,-rhs_w)
        temp_w = x_w[1]
        h_we = temp_w
        return h_we
    elif (drvw,drve,drvp) == (1,0,1):
        fh_aa = _air_f(2,0,0,airf,temp,dhum)
        fh_at = _air_f(1,1,0,airf,temp,dhum)
        fh_ad = _air_f(1,0,1,airf,temp,dhum)
        gi_p = _ice_g(0,1,temp,pres)
        gi_tp = _ice_g(1,1,temp,pres)
        lhs_p = numpy.array([1.,gi_p,0.])
        rhs_p = numpy.array([0.,0.,-(1-w)*gi_tp])
        x_p = numpy.linalg.solve(jac,lhs_p-rhs_p)
        h_wx = numpy.array([fh_aa,fh_at,fh_ad])
        h_wp = h_wx.dot(x_p)
        return h_wp
    elif (drvw,drve,drvp) == (0,2,0):
        lhs_e = numpy.array([0.,0.,1.])
        x_e = numpy.linalg.solve(jac,lhs_e)
        t_e = x_e[1]
        h_ee = t_e
        return h_ee
    elif (drvw,drve,drvp) == (0,1,1):
        gi_p = _ice_g(0,1,temp,pres)
        gi_tp = _ice_g(1,1,temp,pres)
        lhs_p = numpy.array([1.,gi_p,0.])
        rhs_p = numpy.array([0.,0.,-(1-w)*gi_tp])
        x_p = numpy.linalg.solve(jac,lhs_p-rhs_p)
        t_p = x_p[1]
        h_ep = t_p
        return h_ep
    elif (drvw,drve,drvp) == (0,0,2):
        gi_p = _ice_g(0,1,temp,pres)
        gi_tp = _ice_g(1,1,temp,pres)
        gi_pp = _ice_g(0,2,temp,pres)
        lhs_p = numpy.array([1.,gi_p,0.])
        rhs_p = numpy.array([0.,0.,-(1-w)*gi_tp])
        x_p = numpy.linalg.solve(jac,lhs_p-rhs_p)
        h_pa = -w/airf*(dhum**(-1) - gi_p)
        h_pt = (1-w)*gi_tp
        h_pd = -w/dhum**2
        h_px = numpy.array([h_pa,h_pt,h_pd])
        h_pp = (1-w)*gi_pp + h_px.dot(x_p)
        return h_pp
    
    # Should not have made it this far!
    errmsg = 'Derivatives {0} not recognized'.format((drvw,drve,drvp))
    raise ValueError(errmsg)


## Thermodynamic properties
def cp(wair,pres,entr=None,temp=None,airf=None,dhum=None,chkvals=False,
    chktol=_CHKTOL,airf0=None,temp0=None,dhum0=None,chkbnd=False,
    mathargs=None):
    """Calculate icy air isobaric heat capacity.
    
    Calculate the isobaric (constant pressure) heat capacity of icy air.
    
    :arg float wair: Total dry fraction in kg/kg.
    :arg float pres: Pressure in Pa.
    :arg entr: Entropy in J/kg/K. If None (default) then `temp` must be
        provided.
    :type entr: float or None
    :arg airf: Humid air dry fraction in kg/kg. If unknown, pass None
        (default) and it will be calculated.
    :type airf: float or None
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
    :arg airf0: Initial guess for the humid air dry fraction in kg/kg.
        If None (default) then `_approx_wep` is used.
    :type airf0: float or None
    :arg temp0: Initial guess for the temperature in K. If None
        (default) then `_approx_wep` is used.
    :type temp0: float or None
    :arg dhum0: Initial guess for the humid air density in kg/m3. If
        None (default) then `_approx_wep` is used.
    :type dhum0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Heat capacity in J/kg/K.
    :raises ValueError: If both entr and temp are None.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> cp(0.9,1e5,entr=-100.)
    1768.514397
    """
    airf, temp, dhum = eq_wpte(wair,pres,entr=entr,temp=temp,airf=airf,
        dhum=dhum,chkvals=chkvals,chktol=chktol,airf0=airf0,temp0=temp0,
        dhum0=dhum0,chkbnd=chkbnd,mathargs=mathargs)
    h_s = temp
    h_ss = iceair_h(0,2,0,wair,pres,temp=temp,airf=airf,dhum=dhum)
    cp = h_s/h_ss
    return cp

def density(wair,pres,entr=None,temp=None,airf=None,dhum=None,
    chkvals=False,chktol=_CHKTOL,airf0=None,temp0=None,dhum0=None,
    chkbnd=False,mathargs=None):
    """Calculate icy air density.
    
    Calculate the total density of icy air; ratio of the mass of dry
    air, water vapour, and ice to the total parcel volume.
    
    :arg float wair: Total dry fraction in kg/kg.
    :arg float pres: Pressure in Pa.
    :arg entr: Entropy in J/kg/K. If None (default) then `temp` must be
        provided.
    :type entr: float or None
    :arg airf: Humid air dry fraction in kg/kg. If unknown, pass None
        (default) and it will be calculated.
    :type airf: float or None
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
    :arg airf0: Initial guess for the humid air dry fraction in kg/kg.
        If None (default) then `_approx_wep` is used.
    :type airf0: float or None
    :arg temp0: Initial guess for the temperature in K. If None
        (default) then `_approx_wep` is used.
    :type temp0: float or None
    :arg dhum0: Initial guess for the humid air density in kg/m3. If
        None (default) then `_approx_wep` is used.
    :type dhum0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Density in kg/m3.
    :raises ValueError: If both entr and temp are None.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> density(0.9,1e5,entr=-100.)
    1.4253189599
    """
    airf, temp, dhum = eq_wpte(wair,pres,entr=entr,temp=temp,airf=airf,
        dhum=dhum,chkvals=chkvals,chktol=chktol,airf0=airf0,temp0=temp0,
        dhum0=dhum0,chkbnd=chkbnd,mathargs=mathargs)
    h_p = iceair_h(0,0,1,wair,pres,temp=temp,airf=airf,dhum=dhum)
    dens = h_p**(-1)
    return dens

def kappa_s(wair,pres,entr=None,temp=None,airf=None,dhum=None,
    chkvals=False,chktol=_CHKTOL,airf0=None,temp0=None,dhum0=None,
    chkbnd=False,mathargs=None):
    """Calculate icy air isentropic compressibility.
    
    Calculate the isentropic compressibility of icy air.
    
    :arg float wair: Total dry fraction in kg/kg.
    :arg float pres: Pressure in Pa.
    :arg entr: Entropy in J/kg/K. If None (default) then `temp` must be
        provided.
    :type entr: float or None
    :arg airf: Humid air dry fraction in kg/kg. If unknown, pass None
        (default) and it will be calculated.
    :type airf: float or None
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
    :arg airf0: Initial guess for the humid air dry fraction in kg/kg.
        If None (default) then `_approx_wep` is used.
    :type airf0: float or None
    :arg temp0: Initial guess for the temperature in K. If None
        (default) then `_approx_wep` is used.
    :type temp0: float or None
    :arg dhum0: Initial guess for the humid air density in kg/m3. If
        None (default) then `_approx_wep` is used.
    :type dhum0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Compressibility in 1/Pa.
    :raises ValueError: If both entr and temp are None.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> kappa_s(0.9,1e5,entr=-100.)
    8.231417515e-06
    """
    airf, temp, dhum = eq_wpte(wair,pres,entr=entr,temp=temp,airf=airf,
        dhum=dhum,chkvals=chkvals,chktol=chktol,airf0=airf0,temp0=temp0,
        dhum0=dhum0,chkbnd=chkbnd,mathargs=mathargs)
    h_p = iceair_h(0,0,1,wair,pres,temp=temp,airf=airf,dhum=dhum)
    h_pp = iceair_h(0,0,2,wair,pres,temp=temp,airf=airf,dhum=dhum)
    kappa = -h_pp / h_p
    return kappa

def lapserate(wair,pres,entr=None,temp=None,airf=None,dhum=None,
    chkvals=False,chktol=_CHKTOL,airf0=None,temp0=None,dhum0=None,
    chkbnd=False,mathargs=None):
    """Calculate icy air adiabatic lapse rate.
    
    Calculate the adiabatic lapse rate of icy air. This is the proper
    'moist' lapse rate when the air is ice-saturated.
    
    :arg float wair: Total dry fraction in kg/kg.
    :arg float pres: Pressure in Pa.
    :arg entr: Entropy in J/kg/K. If None (default) then `temp` must be
        provided.
    :type entr: float or None
    :arg airf: Humid air dry fraction in kg/kg. If unknown, pass None
        (default) and it will be calculated.
    :type airf: float or None
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
    :arg airf0: Initial guess for the humid air dry fraction in kg/kg.
        If None (default) then `_approx_wep` is used.
    :type airf0: float or None
    :arg temp0: Initial guess for the temperature in K. If None
        (default) then `_approx_wep` is used.
    :type temp0: float or None
    :arg dhum0: Initial guess for the humid air density in kg/m3. If
        None (default) then `_approx_wep` is used.
    :type dhum0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Lapse rate in K/Pa.
    :raises ValueError: If both entr and temp are None.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> lapserate(0.9,1e5,entr=-100.)
    4.420925682e-04
    """
    h_sp = iceair_h(0,1,1,wair,pres,entr=entr,airf=airf,temp=temp,dhum=dhum,
        chkvals=chkvals,chktol=chktol,airf0=airf0,temp0=temp0,dhum0=dhum0,
        chkbnd=chkbnd,mathargs=mathargs)
    gamma = h_sp
    return gamma

def temperature(wair,pres,entr=None,temp=None,airf=None,dhum=None,
    chkvals=False,chktol=_CHKTOL,airf0=None,temp0=None,dhum0=None,
    chkbnd=False,mathargs=None):
    """Calculate icy air temperature.
    
    Calculate the temperature of icy air.
    
    :arg float wair: Total dry fraction in kg/kg.
    :arg float pres: Pressure in Pa.
    :arg entr: Entropy in J/kg/K. If None (default) then `temp` must be
        provided.
    :type entr: float or None
    :arg airf: Humid air dry fraction in kg/kg. If unknown, pass None
        (default) and it will be calculated.
    :type airf: float or None
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
    :arg airf0: Initial guess for the humid air dry fraction in kg/kg.
        If None (default) then `_approx_wep` is used.
    :type airf0: float or None
    :arg temp0: Initial guess for the temperature in K. If None
        (default) then `_approx_wep` is used.
    :type temp0: float or None
    :arg dhum0: Initial guess for the humid air density in kg/m3. If
        None (default) then `_approx_wep` is used.
    :type dhum0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Temperature in K.
    :raises ValueError: If both entr and temp are None.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> temperature(0.9,1e5,entr=-100.)
    270.3836801
    """
    airf, temp, dhum = eq_wpte(wair,pres,entr=entr,temp=temp,airf=airf,
        dhum=dhum,chkvals=chkvals,chktol=chktol,airf0=airf0,temp0=temp0,
        dhum0=dhum0,chkbnd=chkbnd,mathargs=mathargs)
    return temp


## Potential thermodynamic functions
def _approx_pot(wair,temp,pres,ppot,airf,dhum):
    """Approximate ATDh2 at WT1P1P2.
    
    Approximate the humid air dry fraction, temperature, and humid air
    density after adiabatic displacement.
    
    :arg float wair: Total dry fraction in kg/kg.
    :arg float temp: In-situ temperature in K.
    :arg float pres: In-situ pressure in Pa.
    :arg float ppot: Potential pressure in Pa.
    :arg float airf: In-situ humid air dry fraction in kg/kg (unused).
    :arg float dhum: In-situ humid air density in kg/m3 (unused).
    :returns: Potential humid air dry fraction, potential temperature,
        and potential humid air density (all in SI units).
    """
    pvsat0 = _PTPE*numpy.exp((_AVI+_BVI)*(1 - _TTP/temp)
        + _BVI*numpy.log(_TTP/temp))
    pvmax0 = pres * (1-wair)/(_EPSW*wair + 1-wair)
    if pvmax0 > pvsat0:
        # Parcel starts saturated
        pv0 = pvsat0
        a0 = (pres-pv0) / (pres-pv0 + _EPSW*pv0)
        ceff0 = (wair*_CDRY + wair*(1-a0)/a0*_CVAP + (1-wair/a0)*_CICE
            + wair*_RWAT*(1-a0)/a0*pres/(pres-pv0)
                * ((_AVI+_BVI)*_TTP/temp - _BVI)**2)
        reff0 = wair*(_RDRY + _RWAT*(1-a0)/a0
            + _RWAT*(1-a0)/a0*pres/(pres-pv0) * ((_AVI+_BVI)*_TTP/temp - _BVI))
        pvmaxt = pvmax0 * (_TTP/temp)**(ceff0/reff0)
        ginv0 = ceff0/reff0
        r = (_AVI+_BVI)/(ginv0+_BVI) - 1
        v = numpy.log((_TTP/temp)**ginv0 * pvmax0/_PTPE)/(ginv0+_BVI)
        if pvmaxt > _PTPE or v <= r:
            # Parcel is always ice-saturated
            tpot = temp * (ppot/pres)**(reff0/ceff0)
            pv2 = _PTPE*numpy.exp((_AVI+_BVI)*(1 - _TTP/tpot)
                + _BVI*numpy.log(_TTP/tpot))
            apot = (ppot-pv2) / (ppot-pv2 + _EPSW*pv2)
        else:
            # Find where parcel de-saturates
            x = maths4.lamb2(v,r)
            ticl = _TTP/x
            picl = pres * (ticl/temp)**ginv
            if ppot < picl:
                # Parcel ends saturated
                tpot = temp * (ppot/pres)**(reff0/ceff0)
                pv2 = _PTPE*numpy.exp((_AVI+_BVI)*(1 - _TTP/tpot)
                    + _BVI*numpy.log(_TTP/tpot))
                apot = (ppot-pv2) / (ppot-pv2 + _EPSW*pv2)
            else:
                # Parcel ends unsaturated
                p1 = picl
                t1 = ticl
                ceff1 = wair*_CDRY + (1-wair)*_CVAP
                reff1 = wair*_RDRY + (1-wair)*_RWAT
                tpot = t1 * (ppot/p1)**(reff1/ceff1)
                apot = wair
    else:
        # Parcel starts unsaturated
        ticl, picl, __ = iceair4a._approx_icl(wair,temp,pres,dhum)
        if ppot < picl:
            # Parcel ends saturated
            p1 = picl
            t1 = ticl
            pv1 = _PTPE*numpy.exp((_AVI+_BVI)*(1 - _TTP/t1)
                + _BVI*numpy.log(_TTP/t1))
            a1 = (p1-pv1) / (p1-pv1 + _EPSW*pv1)
            ceff1 = (wair*_CDRY + (1-wair)*_CVAP
                + (1-wair)*_RWAT*p1/(p1-pv1) * ((_AVI+_BVI)*_TTP/t1 - _BVI)**2)
            reff1 = (wair*_RDRY + (1-wair)*_RWAT
                + (1-wair)*_RWAT*p1/(p1-pv1) * ((_AVI+_BVI)*_TTP/t1 - _BVI))
            tpot = t1 * (ppot/p1)**(reff1/ceff1)
            pv2 = _PTPE*numpy.exp((_AVI+_BVI)*(1 - _TTP/tpot)
                + _BVI*numpy.log(_TTP/tpot))
            apot = (ppot-pv2) / (ppot-pv2 + _EPSW*pv2)
        else:
            # Parcel ends unsaturated
            ceff1 = wair*_CDRY + (1-wair)*_CVAP
            reff1 = wair*_RDRY + (1-wair)*_RWAT
            tpot = temp * (ppot/pres)**(reff1/ceff1)
            apot = wair
    dhpot = ppot/(_RDRY*tpot) / (apot + (1-apot)/_EPSW)
    return apot, tpot, dhpot

def _diff_pot(a2,t2,d2,wair,temp,pres,ppot,airf,dhum):
    """Calculate icy air disequilibrium at WT1P1P2.
    
    Calculate both sides of the equations
    
        given potential pressure = potential humid air pressure
        potential ice chemical potential = potential vapour potential
        initial entropy of icy air = potential entropy of icy air
    
    and their Jacobians with respect to potential humid air dry
    fraction, potential temperature, and potential humid air density.
    Solving these equations produces equilibrium values at the given
    total dry fraction, in-situ temperature, in-situ pressure, and
    potential pressure.
    
    :arg float a2: Potential humid air dry fraction in kg/kg.
    :arg float t2: Potential temperature in K.
    :arg float d2: Potential humid air density in kg/m3.
    :arg float wair: Total dry fraction in kg/kg.
    :arg float temp: In-situ temperature in K.
    :arg float pres: In-situ pressure in Pa.
    :arg float ppot: Potential pressure in Pa.
    :arg float airf: In-situ humid air dry fraction in kg/kg.
    :arg float dhum: In-situ humid air density in kg/m3.
    :returns: Left-hand side of the equations, right-hand side,
        Jacobian of LHS, and Jacobian of RHS.
    :rtype: tuple(array(float))
    """
    ph2 = _eq_pressure(0,0,0,a2,t2,d2)
    gi2 = _ice_g(0,0,t2,ppot)
    gv2 = _eq_vappot(0,0,0,a2,t2,d2)
    sh1 = -_air_f(0,1,0,airf,temp,dhum)
    si1 = -_ice_g(1,0,temp,pres)
    s1 = wair/airf*sh1 + (1-wair/airf)*si1
    sh2 = -_air_f(0,1,0,a2,t2,d2)
    si2 = -_ice_g(1,0,t2,ppot)
    s2 = wair/a2*sh2 + (1-wair/a2)*si2
    lhs = numpy.array([ppot, gi2, s1])
    rhs = numpy.array([ph2, gv2, s2])
    
    ph2_a = _eq_pressure(1,0,0,a2,t2,d2)
    ph2_t = _eq_pressure(0,1,0,a2,t2,d2)
    ph2_d = _eq_pressure(0,0,1,a2,t2,d2)
    gi2_t = _ice_g(1,0,t2,ppot)
    gv2_a = _eq_vappot(1,0,0,a2,t2,d2)
    gv2_t = _eq_vappot(0,1,0,a2,t2,d2)
    gv2_d = _eq_vappot(0,0,1,a2,t2,d2)
    sh2_a = -_air_f(1,1,0,a2,t2,d2)
    sh2_t = -_air_f(0,2,0,a2,t2,d2)
    sh2_d = -_air_f(0,1,1,a2,t2,d2)
    si2_t = -_ice_g(2,0,t2,ppot)
    s2_a = -wair/a2**2*(sh2 - a2*sh2_a - si2)
    s2_t = wair/a2*sh2_t + (1-wair/a2)*si2_t
    s2_d = wair/a2*sh2_d
    dlhs = numpy.array([[0.,0.,0.], [0.,gi2_t,0.], [0.,0.,0.]])
    drhs = numpy.array([[ph2_a,ph2_t,ph2_d], [gv2_a,gv2_t,gv2_d],
        [s2_a,s2_t,s2_d]])
    return lhs, rhs, dlhs, drhs

def eq_pot(wair,temp,pres,ppot,airf=None,dhum=None,apot=None,tpot=None,
    dhpot=None,chkvals=False,chktol=_CHKTOL,airf0=None,dhum0=None,
    apot0=None,tpot0=None,dhpot0=None,chkbnd=False,mathargs=None):
    """Get primary variables at WT1P1P2.
    
    Get the values of the equilibrium potential humid air dry fraction,
    potential temperature, and potential humid air density for the given
    total dry fraction, in-situ temperature, in-situ pressure, and
    potential pressure.
    
    If the calculation has already been done, the results can be passed
    to avoid unnecessary repeat calculations. If enough values are
    passed, they will be checked for consistency if chkvals is True.
    
    :arg float wair: Total dry fraction in kg/kg.
    :arg float temp: In-situ temperature in K.
    :arg float pres: In-situ pressure in Pa.
    :arg float ppot: Potential pressure in Pa.
    :arg airf: In-situ humid air dry fraction in kg/kg. If unknown, pass
        None (default) and it will be calculated.
    :type airf: float or None
    :arg dhum: In-situ humid air density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dhum: float or None
    :arg apot: Potential humid air dry fraction in kg/kg. If unknown,
        pass None (default) and it will be calculated.
    :type apot: float or None
    :arg tpot: Potential temperature in K. If unknown, pass None
        (default) and it will be calculated.
    :type tpot: float or None
    :arg dhpot: Potential humid air density in kg/m3. If unknown, pass
        None (default) and it will be calculated.
    :type dhpot: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg airf0: Initial guess for the in-situ humid air dry fraction in
        kg/kg. If None (default) then `iceair4a._approx_tp` is used.
    :type airf0: float or None
    :arg dhum0: Initial guess for the in-situ humid air density in
        kg/m3. If None (default) then `iceair4a._approx_tp` is used.
    :type dhum0: float or None
    :arg apot0: Initial guess for the potential humid air dry fraction
        in kg/kg. If None (default) then `_approx_pot` is used.
    :type apot0: float or None
    :arg temp0: Initial guess for the potential temperature in K. If
        None (default) then `_approx_pot` is used.
    :type temp0: float or None
    :arg dhpot0: Initial guess for the potential humid air density in
        kg/m3. If None (default) then `_approx_pot` is used.
    :type dhpot0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: In-situ humid air dry fraction, in-situ humid air density,
        potential humid air dry fraction, potential temperature, and
        potential humid air density (all in SI units).
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    """
    if any(val is None for val in (airf,dhum)):
        airf, __, __, dhum = iceair4a.eq_atpe(temp=temp,pres=pres,airf0=airf0,
            dhum0=dhum0,mathargs=mathargs)
    if any(val is None for val in (apot,tpot,dhpot)):
        x0 = (apot0,tpot0,dhpot0)
        fargs = (wair,temp,pres,ppot,airf,dhum)
        if mathargs is None:
            mathargs = dict()
        x1 = _newton(_diff_pot,x0,_approx_pot,fargs=fargs,**mathargs)
        apot, tpot, dhpot = x1
    
    _chkhumbnds(airf,temp,dhum,chkbnd=chkbnd)
    _chkhumbnds(apot,tpot,dhpot,chkbnd=chkbnd)
    _chkicebnds(tpot,ppot,chkbnd=chkbnd)
    if not chkvals:
        return airf, dhum, apot, tpot, dhpot
    
    lhs1, rhs1, __, __ = iceair4a._diff_tp(airf,dhum,temp,pres)
    lhs2, rhs2, __, __ = _diff_pot(apot,tpot,dhpot,wair,temp,pres,ppot,airf,
        dhum)
    lhs = numpy.concatenate((lhs1,lhs2))
    rhs = numpy.concatenate((rhs1,rhs2))
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
    return airf, dhum, apot, tpot, dhpot

def potdensity(wair,temp,pres,ppot,airf=None,dhum=None,apot=None,
    tpot=None,dhpot=None,chkvals=False,chktol=_CHKTOL,airf0=None,
    dhum0=None,apot0=None,tpot0=None,dhpot0=None,chkbnd=False,
    mathargs=None):
    """Calculate icy air potential density.
    
    Calculate the potential density of icy air.
    
    :arg float wair: Total dry fraction in kg/kg.
    :arg float temp: In-situ temperature in K.
    :arg float pres: In-situ pressure in Pa.
    :arg float ppot: Potential pressure in Pa.
    :arg airf: In-situ humid air dry fraction in kg/kg. If unknown, pass
        None (default) and it will be calculated.
    :type airf: float or None
    :arg dhum: In-situ humid air density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dhum: float or None
    :arg apot: Potential humid air dry fraction in kg/kg. If unknown,
        pass None (default) and it will be calculated.
    :type apot: float or None
    :arg tpot: Potential temperature in K. If unknown, pass None
        (default) and it will be calculated.
    :type tpot: float or None
    :arg dhpot: Potential humid air density in kg/m3. If unknown, pass
        None (default) and it will be calculated.
    :type dhpot: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg airf0: Initial guess for the in-situ humid air dry fraction in
        kg/kg. If None (default) then `iceair4a._approx_tp` is used.
    :type airf0: float or None
    :arg dhum0: Initial guess for the in-situ humid air density in
        kg/m3. If None (default) then `iceair4a._approx_tp` is used.
    :type dhum0: float or None
    :arg apot0: Initial guess for the potential humid air dry fraction
        in kg/kg. If None (default) then `_approx_pot` is used.
    :type apot0: float or None
    :arg temp0: Initial guess for the potential temperature in K. If
        None (default) then `_approx_pot` is used.
    :type temp0: float or None
    :arg dhpot0: Initial guess for the potential humid air density in
        kg/m3. If None (default) then `_approx_pot` is used.
    :type dhpot0: float or None
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
    
    >>> potdensity(0.9,230.,5e4,1e5)
    1.45048110422
    """
    airf, dhum, apot, tpot, dhpot = eq_pot(wair,temp,pres,ppot,airf=airf,
        dhum=dhum,apot=apot,tpot=tpot,dhpot=dhpot,chkvals=chkvals,chktol=chktol,
        airf0=airf0,dhum0=dhum0,apot0=apot0,tpot0=tpot0,dhpot0=dhpot0,
        chkbnd=chkbnd,mathargs=mathargs)
    hp_p = iceair_h(0,0,1,wair,ppot,temp=tpot,airf=apot,dhum=dhpot)
    dpot = hp_p**(-1)
    return dpot

def potenthalpy(wair,temp,pres,ppot,airf=None,dhum=None,apot=None,
    tpot=None,dhpot=None,chkvals=False,chktol=_CHKTOL,airf0=None,
    dhum0=None,apot0=None,tpot0=None,dhpot0=None,chkbnd=False,
    mathargs=None):
    """Calculate icy air potential enthalpy.
    
    Calculate the potential enthalpy of icy air.
    
    :arg float wair: Total dry fraction in kg/kg.
    :arg float temp: In-situ temperature in K.
    :arg float pres: In-situ pressure in Pa.
    :arg float ppot: Potential pressure in Pa.
    :arg airf: In-situ humid air dry fraction in kg/kg. If unknown, pass
        None (default) and it will be calculated.
    :type airf: float or None
    :arg dhum: In-situ humid air density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dhum: float or None
    :arg apot: Potential humid air dry fraction in kg/kg. If unknown,
        pass None (default) and it will be calculated.
    :type apot: float or None
    :arg tpot: Potential temperature in K. If unknown, pass None
        (default) and it will be calculated.
    :type tpot: float or None
    :arg dhpot: Potential humid air density in kg/m3. If unknown, pass
        None (default) and it will be calculated.
    :type dhpot: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg airf0: Initial guess for the in-situ humid air dry fraction in
        kg/kg. If None (default) then `iceair4a._approx_tp` is used.
    :type airf0: float or None
    :arg dhum0: Initial guess for the in-situ humid air density in
        kg/m3. If None (default) then `iceair4a._approx_tp` is used.
    :type dhum0: float or None
    :arg apot0: Initial guess for the potential humid air dry fraction
        in kg/kg. If None (default) then `_approx_pot` is used.
    :type apot0: float or None
    :arg temp0: Initial guess for the potential temperature in K. If
        None (default) then `_approx_pot` is used.
    :type temp0: float or None
    :arg dhpot0: Initial guess for the potential humid air density in
        kg/m3. If None (default) then `_approx_pot` is used.
    :type dhpot0: float or None
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
    
    >>> potenthalpy(0.9,230.,5e4,1e5)
    -35781.25645
    """
    airf, dhum, apot, tpot, dhpot = eq_pot(wair,temp,pres,ppot,airf=airf,
        dhum=dhum,apot=apot,tpot=tpot,dhpot=dhpot,chkvals=chkvals,chktol=chktol,
        airf0=airf0,dhum0=dhum0,apot0=apot0,tpot0=tpot0,dhpot0=dhpot0,
        chkbnd=chkbnd,mathargs=mathargs)
    hp = iceair_h(0,0,0,wair,ppot,temp=tpot,airf=apot,dhum=dhpot)
    return hp

def pottemp(wair,temp,pres,ppot,airf=None,dhum=None,apot=None,tpot=None,
    dhpot=None,chkvals=False,chktol=_CHKTOL,airf0=None,dhum0=None,
    apot0=None,tpot0=None,dhpot0=None,chkbnd=False,mathargs=None):
    """Calculate icy air potential temperature.
    
    Calculate the potential temperature of icy air.
    
    :arg float wair: Total dry fraction in kg/kg.
    :arg float temp: In-situ temperature in K.
    :arg float pres: In-situ pressure in Pa.
    :arg float ppot: Potential pressure in Pa.
    :arg airf: In-situ humid air dry fraction in kg/kg. If unknown, pass
        None (default) and it will be calculated.
    :type airf: float or None
    :arg dhum: In-situ humid air density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dhum: float or None
    :arg apot: Potential humid air dry fraction in kg/kg. If unknown,
        pass None (default) and it will be calculated.
    :type apot: float or None
    :arg tpot: Potential temperature in K. If unknown, pass None
        (default) and it will be calculated.
    :type tpot: float or None
    :arg dhpot: Potential humid air density in kg/m3. If unknown, pass
        None (default) and it will be calculated.
    :type dhpot: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg airf0: Initial guess for the in-situ humid air dry fraction in
        kg/kg. If None (default) then `iceair4a._approx_tp` is used.
    :type airf0: float or None
    :arg dhum0: Initial guess for the in-situ humid air density in
        kg/m3. If None (default) then `iceair4a._approx_tp` is used.
    :type dhum0: float or None
    :arg apot0: Initial guess for the potential humid air dry fraction
        in kg/kg. If None (default) then `_approx_pot` is used.
    :type apot0: float or None
    :arg temp0: Initial guess for the potential temperature in K. If
        None (default) then `_approx_pot` is used.
    :type temp0: float or None
    :arg dhpot0: Initial guess for the potential humid air density in
        kg/m3. If None (default) then `_approx_pot` is used.
    :type dhpot0: float or None
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
    
    >>> pottemp(0.9,230.,5e4,1e5)
    266.1052089
    """
    airf, dhum, apot, tpot, dhpot = eq_pot(wair,temp,pres,ppot,airf=airf,
        dhum=dhum,apot=apot,tpot=tpot,dhpot=dhpot,chkvals=chkvals,chktol=chktol,
        airf0=airf0,dhum0=dhum0,apot0=apot0,tpot0=tpot0,dhpot0=dhpot0,
        chkbnd=chkbnd,mathargs=mathargs)
    return tpot

