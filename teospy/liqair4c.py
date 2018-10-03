"""Wet air enthalpy and potential temperature.

This module provides the specific enthalpy of liquid water-saturated
(wet) air. It also provides the potential temperature and other
adiabatic properties. The primary variables are the total dry air mass
fraction, entropy, and pressure. As in `liqair4b`, the total dry
fraction here is the ratio of dry air mass to total parcel mass,
including liquid.

:Examples:

>>> liqair_h(0,0,0,0.5,100.,1e5)
26898.5215492
>>> liqair_h(0,1,0,0.5,100.,1e5)
280.393544899
>>> liqair_h(0,1,1,0.5,100.,1e5)
1.55067379031e-4
>>> cp(0.5,1e5,entr=100.)
3144.21265404
>>> lapserate(0.5,1e5,entr=100.)
1.55067379031e-4
>>> temperature(0.5,1e5,entr=100.)
280.393544899
>>> potdensity(0.5,300.,1e4,1e5)
1.22550664945
>>> potenthalpy(0.5,300.,1e4,1e5)
655155.797982
>>> pottemp(0.5,300.,1e4,1e5)
348.222379217

:Functions:

* :func:`eq_wpte`: Calculate wet air equilibrium properties at total dry
  fraction, pressure, and either entropy or temperature.
* :func:`liqair_h`: Wet air specific enthalpy with derivatives.
* :func:`cp`: Wet air isobaric heat capacity.
* :func:`density`: Wet air density.
* :func:`kappa_s`: Wet air isentropic compressibility.
* :func:`lapserate`: Wet air adiabatic lapse rate.
* :func:`temperature`: Wet air temperature.
* :func:`eq_pot`: Calculate wet air equilibrium properties under
  adiabatic displacement.
* :func:`potdensity`: Wet air potential density.
* :func:`potenthalpy`: Wet air potential enthalpy.
* :func:`pottemp`: Wet air potential temperature.

"""

__all__ = ['eq_wpte','liqair_h',
    'cp','density','kappa_s','lapserate','temperature',
    'eq_pot','potdensity','potenthalpy','pottemp']

import numpy
import warnings
import constants0
import flu1
import air2
import flu2
import maths3
import air3c
import flu3a
import maths4
import liqair4a
import liqair4b

_CHKTOL = constants0.CHKTOL
_RDRY = constants0.RDRY
_RWAT = constants0.RWAT
_TTP = constants0.TTP
_PTPE = constants0.PTPE
_LLVTP = constants0.LLVTP
_PATM = constants0.PATM
_TCELS = constants0.TCELS
_CDRY = constants0.CDRY
_CVAP = constants0.CVAP
_CLIQ = constants0.CLIQ
_chkhumbnds = constants0.chkhumbnds
_chkflubnds = constants0.chkflubnds
_flu_f = flu1.flu_f
_air_f = air2.air_f
_air_eq_pressure = air2.eq_pressure
_air_eq_vappot = air2.eq_vappot
_flu_eq_chempot = flu2.eq_chempot
_flu_eq_pressure = flu2.eq_pressure
_newton = maths3.newton
_dliq_default = flu3a._dliq_default
_AVL = _LLVTP / (_RWAT*_TTP)
_BVL = (_CLIQ-_CVAP) / _RWAT
_EPSW = _RDRY/_RWAT


## Equilibrium functions
def _approx_wep(wair,entr,pres):
    """Approximate ATDhDl at WEP.
    
    Approximate the humid air dry fraction, temperature, humid air
    density, and liquid water density at the given total dry fraction,
    entropy, and pressure.
    
    :arg float wair: Total dry fraction in kg/kg.
    :arg float entr: Entropy in J/kg/K.
    :arg float pres: Pressure in Pa.
    :returns: Humid air dry fraction, temperature, humid air density,
        and liquid water density (all in SI units).
    """
    pvmax = pres * (1-wair) / (1-wair + _EPSW*wair)
    if pvmax >= _PTPE:
        # Saturation would start at _TTP; use saturated heat capacity at _TTP
        a_t = (pres - _PTPE)/(pres - _PTPE + _EPSW*_PTPE)
        s_t = (wair*_CDRY*numpy.log(_TTP/_TCELS)
            - wair*_RDRY*numpy.log((pres-_PTPE)/_PATM)
            + wair*(a_t**(-1)-1)*_RWAT*_AVL)
        c_t = (wair*_CDRY + wair*(1-a_t)/a_t*_CVAP + (1-wair/a_t)*_CLIQ
            + wair*_RWAT*(1-a_t)/a_t*pres/(pres-_PTPE)*_AVL**2)
        temp = _TTP * numpy.exp((entr-s_t)/c_t)
    else:
        # Get approximate saturation temperature
        v = numpy.log(pres*(1-wair)/(_PTPE*(_EPSW*wair + 1-wair)))/_BVL
        r = _AVL/BVL
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
                    * ((_AVL+_BVL)*_TTP/tsat - _BVL)**2)
            temp = tsat * numpy.exp(-(ssat-entr)/csat)
    pvap = _PTPE * numpy.exp((_AVL+_BVL)*(1 - _TTP/temp)
        - _BVL*numpy.log(temp/_TTP))
    airf = (pres - pvap) / (pres - pvap + _EPSW*pvap)
    dhum = pres/(_RDRY*temp) / (airf + (1-airf)/_EPSW)
    dliq = _dliq_default(temp,pres)
    return airf, temp, dhum, dliq

def _diff_wep(a,t,dh,dl,wair,entr,pres):
    """Calculate wet air disequilibrium at WEP.
    
    Calculate both sides of the equations
    
        given pressure = pressure of humid air
        given pressure = pressure of liquid water
        liquid water chemical potential = water vapour potential
        given entropy = entropy of wet air
    
    and their Jacobians with respect to the humid air dry fraction,
    temperature, humid air density, and liquid water density. Solving
    these equations gives equilibrium values at the given total dry
    fraction, entropy, and pressure.
    
    :arg float a: Humid air dry fraction in kg/kg.
    :arg float t: Temperature in K.
    :arg float dh: Humid air density in kg/m3.
    :arg float dl: Liquid water density in kg/m3.
    :arg float wair: Total dry fraction in kg/kg.
    :arg float entr: Entropy in J/kg/K.
    :arg float pres: Pressure in Pa.
    :returns: Left-hand side of the equation, right-hand side,
        Jacobian of LHS, and Jacobian of RHS.
    :rtype: tuple(array(float))
    """
    ph = _air_eq_pressure(0,0,0,a,t,dh)
    pl = _flu_eq_pressure(0,0,t,dl)
    gl = _flu_eq_chempot(0,0,t,dl)
    gv = _air_eq_vappot(0,0,0,a,t,dh)
    sh = -_air_f(0,1,0,a,t,dh)
    sl = -_flu_f(1,0,t,dl)
    s = wair/a*sh + (1-wair/a)*sl
    lhs = numpy.array([pres, pres, gl, entr])
    rhs = numpy.array([ph, pl, gv, s])
    
    ph_a = _air_eq_pressure(1,0,0,a,t,dh)
    ph_t = _air_eq_pressure(0,1,0,a,t,dh)
    ph_d = _air_eq_pressure(0,0,1,a,t,dh)
    pl_t = _flu_eq_pressure(1,0,t,dl)
    pl_d = _flu_eq_pressure(0,1,t,dl)
    gl_t = _flu_eq_chempot(1,0,t,dl)
    gl_d = _flu_eq_chempot(0,1,t,dl)
    gv_a = _air_eq_vappot(1,0,0,a,t,dh)
    gv_t = _air_eq_vappot(0,1,0,a,t,dh)
    gv_d = _air_eq_vappot(0,0,1,a,t,dh)
    sh_a = -_air_f(1,1,0,a,t,dh)
    sh_t = -_air_f(0,2,0,a,t,dh)
    sh_d = -_air_f(0,1,1,a,t,dh)
    sl_t = -_flu_f(2,0,t,dl)
    sl_d = -_flu_f(1,1,t,dl)
    s_a = -wair/a**2 * (sh - a*sh_a - sl)
    s_t = wair/a*sh_t + (1-wair/a)*sl_t
    s_dh = wair/a*sh_d
    s_dl = (1-wair/a)*sl_d
    dlhs = numpy.array([[0.,0.,0.,0.], [0.,0.,0.,0.], [0.,gl_t,0.,gl_d],
        [0.,0.,0.,0.]])
    drhs = numpy.array([[ph_a,ph_t,ph_d,0.], [0.,pl_t,0.,pl_d],
        [gv_a,gv_t,gv_d,0.], [s_a,s_t,s_dh,s_dl]])
    return lhs, rhs, dlhs, drhs

def eq_wpte(wair,pres,entr=None,airf=None,temp=None,dhum=None,dliq=None,
    chkvals=False,chktol=_CHKTOL,airf0=None,temp0=None,dhum0=None,
    dliq0=None,chkbnd=False,mathargs=None):
    """Get primary wet air variables at WP and T/E.
    
    Get the values of all primary variables for wet air at the given
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
    :arg dliq: Liquid water density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dliq: float or None
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
    :arg dliq0: Initial guess for the liquid water density in kg/m3. If
        None (default) then `_approx_wep` is used.
    :type dliq0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Humid air dry fraction, temperature, humid air density,
        and liquid water density (all in SI units).
    :raises ValueError: If both entr and temp are None.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    """
    if entr is None and temp is None:
        errmsg = 'Must provide at least one of entr or temp'
        raise ValueError(errmsg)
    if temp is None:
        x0 = (airf0,temp0,dhum0,dliq0)
        fargs = (wair,entr,pres)
        if mathargs is None:
            mathargs = dict()
        x1 = _newton(_diff_wep,x0,_approx_wep,fargs=fargs,**mathargs)
        airf, temp, dhum, dliq = x1
    elif any(val is None for val in (airf,dhum,dliq)):
        airf, __, __, dhum, dliq = liqair4a.eq_atpe(temp=temp,pres=pres,
            chkvals=chkvals,chktol=chktol,airf0=airf0,dhum0=dhum0,dliq0=dliq0,
            mathargs=mathargs)
    
    _chkhumbnds(airf,temp,dhum,chkbnd=chkbnd)
    _chkflubnds(temp,dliq,chkbnd=chkbnd)
    if not chkvals:
        return airf, temp, dhum, dliq
    
    if entr is None:
        entr = liqair4b.entropy(wair,temp,pres,airf=airf,dhum=dhum,dliq=dliq)
    lhs, rhs, __, __ = _diff_wep(airf,temp,dhum,dliq,wair,entr,pres)
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
    return airf, temp, dhum, dliq


## Enthalpy
def liqair_h(drvw,drve,drvp,wair,pres,entr=None,temp=None,airf=None,
    dhum=None,dliq=None,chkvals=False,chktol=_CHKTOL,airf0=None,
    temp0=None,dhum0=None,dliq0=None,chkbnd=False,mathargs=None):
    """Calculate wet air enthalpy with derivatives.
    
    Calculate the specific enthalpy of wet air or its derivatives with
    respect to total dry air fraction, entropy, and pressure.
    
    :arg float drvw: Number of total dry fraction derivatives.
    :arg float drve: Number of entropy derivatives.
    :arg float drvp: Number of pressure derivatives.
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
    :arg dliq: Liquid water density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dliq: float or None
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
    :arg dliq0: Initial guess for the liquid water density in kg/m3. If
        None (default) then `_approx_wep` is used.
    :type dliq0: float or None
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
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    :raises ValueError: If any of (drvw,drve,drvp) are negative, or if
        (drvw+drve+drvp)>2.
    
    :Examples:
    
    >>> liqair_h(0,0,0,0.5,100.,1e5)
    26898.5215492
    >>> liqair_h(1,0,0,0.5,100.,1e5)
    -1681.79366113
    >>> liqair_h(0,1,0,0.5,100.,1e5)
    280.393544899
    >>> liqair_h(0,0,1,0.5,100.,1e5)
    0.40687293002
    >>> liqair_h(2,0,0,0.5,100.,1e5)
    35.7689708915
    >>> liqair_h(1,1,0,0.5,100.,1e5)
    1.78599925196
    >>> liqair_h(1,0,1,0.5,100.,1e5)
    0.811745643965
    >>> liqair_h(0,2,0,0.5,100.,1e5)
    8.91776656830e-2
    >>> liqair_h(0,1,1,0.5,100.,1e5)
    1.55067379031e-4
    >>> liqair_h(0,0,2,0.5,100.,1e5)
    -3.83770118470e-6
    """
    if any(drv < 0 for drv in (drvw,drve,drvp)) or (drvw+drve+drvp) > 2:
        errmsg = 'Derivatives {0} not recognized'.format((drvw,drve,drvp))
        raise ValueError(errmsg)
    airf, temp, dhum, dliq = eq_wpte(wair,pres,entr=entr,temp=temp,
        airf=airf,dhum=dhum,dliq=dliq,chkvals=chkvals,chktol=chktol,
        airf0=airf0,temp0=temp0,dhum0=dhum0,dliq0=dliq0,chkbnd=chkbnd,
        mathargs=mathargs)
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
        fh = _air_f(0,0,0,airf,temp,dhum)
        fh_t = _air_f(0,1,0,airf,temp,dhum)
        fh_d = _air_f(0,0,1,airf,temp,dhum)
        fl = _flu_f(0,0,temp,dliq)
        fl_t = _flu_f(1,0,temp,dliq)
        fl_d = _flu_f(0,1,temp,dliq)
        hh = fh - temp*fh_t + dhum*fh_d
        hl = fl - temp*fl_t + dliq*fl_d
        h = w*hh + (1-w)*hl
        return h
    elif (drvw,drve,drvp) == (1,0,0):
        fh_a = _air_f(1,0,0,airf,temp,dhum)
        h_w = fh_a
        return h_w
    elif (drvw,drve,drvp) == (0,1,0):
        h_e = temp
        return h_e
    elif (drvw,drve,drvp) == (0,0,1):
        h_p = w/dhum + (1-w)/dliq
        return h_p
    
    # Higher-order derivatives require inversion
    if entr is None:
        entr = liqair4b.entropy(wair,temp,pres,airf=airf,dhum=dhum,dliq=dliq)
    __, __, dlhs, drhs = _diff_wep(airf,temp,dhum,dliq,wair,entr,pres)
    ppge_x = drhs - dlhs
    if (drvw,drve,drvp) == (2,0,0):
        fh_t = _air_f(0,1,0,airf,temp,dhum)
        fl_t = _flu_f(1,0,temp,dliq)
        ppge_w = numpy.array([0.,0.,0.,(-fh_t + fl_t)/airf])
        x_w = numpy.linalg.solve(ppge_x,-ppge_w)
        
        h_wa = _air_f(2,0,0,airf,temp,dhum)
        h_wt = _air_f(1,1,0,airf,temp,dhum)
        h_wh = _air_f(1,0,1,airf,temp,dhum)
        h_wx = numpy.array([h_wa,h_wt,h_wh,0.])
        h_ww = h_wx.dot(x_w)
        return h_ww
    elif (drvw,drve,drvp) == (1,1,0):
        fh_t = _air_f(0,1,0,airf,temp,dhum)
        fl_t = _flu_f(1,0,temp,dliq)
        ppge_w = numpy.array([0.,0.,0.,(-fh_t + fl_t)/airf])
        x_w = numpy.linalg.solve(ppge_x,-ppge_w)
        t_w = x_w[1]
        h_we = t_w
        return h_we
    elif (drvw,drve,drvp) == (1,0,1):
        fh_t = _air_f(0,1,0,airf,temp,dhum)
        fl_t = _flu_f(1,0,temp,dliq)
        ppge_w = numpy.array([0.,0.,0.,(-fh_t + fl_t)/airf])
        x_w = numpy.linalg.solve(ppge_x,-ppge_w)
        
        h_pa = -w/airf*(dhum**(-1) - dliq**(-1))
        h_ph = -w/dhum**2
        h_pl = -(1-w)/dliq**2
        h_px = numpy.array([h_pa,0.,h_ph,h_pl])
        h_wp = (dhum**(-1) - dliq**(-1))/airf + h_px.dot(x_w)
        return h_wp
    elif (drvw,drve,drvp) == (0,2,0):
        ppge_e = numpy.array([0.,0.,0.,1.])
        x_e = numpy.linalg.solve(ppge_x,ppge_e)
        t_e = x_e[1]
        h_ee = t_e
        return h_ee
    elif (drvw,drve,drvp) == (0,1,1):
        ppge_p = numpy.array([1.,1.,0.,0.])
        x_p = numpy.linalg.solve(ppge_x,ppge_p)
        t_p = x_p[1]
        h_ep = t_p
        return h_ep
    elif (drvw,drve,drvp) == (0,0,2):
        ppge_p = numpy.array([1.,1.,0.,0.])
        x_p = numpy.linalg.solve(ppge_x,ppge_p)
        
        h_pa = -w/airf*(dhum**(-1) - dliq**(-1))
        h_ph = -w/dhum**2
        h_pl = -(1-w)/dliq**2
        h_px = numpy.array([h_pa,0.,h_ph,h_pl])
        h_pp = h_px.dot(x_p)
        return h_pp
    
    # Should not have made it this far!
    errmsg = 'Derivatives {0} not recognized'.format((drvw,drve,drvp))
    raise ValueError(errmsg)


## Thermodynamic properties
def cp(wair,pres,entr=None,temp=None,airf=None,dhum=None,dliq=None,
    chkvals=False,chktol=_CHKTOL,airf0=None,temp0=None,dhum0=None,
    dliq0=None,chkbnd=False,mathargs=None):
    """Calculate wet air isobaric heat capacity.
    
    Calculate the isobaric heat capacity of wet air.
    
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
    :arg dliq: Liquid water density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dliq: float or None
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
    :arg dliq0: Initial guess for the liquid water density in kg/m3. If
        None (default) then `_approx_wep` is used.
    :type dliq0: float or None
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
    
    >>> cp(0.5,1e5,entr=100.)
    3144.21265404
    """
    airf, temp, dhum, dliq = eq_wpte(wair,pres,entr=entr,temp=temp,airf=airf,
        dhum=dhum,dliq=dliq,chkvals=chkvals,chktol=chktol,airf0=airf0,
        temp0=temp0,dhum0=dhum0,dliq0=dliq0,chkbnd=chkbnd,mathargs=mathargs)
    h_ee = liqair_h(0,2,0,wair,pres,airf=airf,temp=temp,dhum=dhum,dliq=dliq)
    cp = temp/h_ee
    return cp

def density(wair,pres,entr=None,temp=None,airf=None,dhum=None,dliq=None,
    chkvals=False,chktol=_CHKTOL,airf0=None,temp0=None,dhum0=None,
    dliq0=None,chkbnd=False,mathargs=None):
    """Calculate wet air density.
    
    Calculate the density of wet air.
    
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
    :arg dliq: Liquid water density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dliq: float or None
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
    :arg dliq0: Initial guess for the liquid water density in kg/m3. If
        None (default) then `_approx_wep` is used.
    :type dliq0: float or None
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
    
    >>> density(0.5,1e5,entr=100.)
    2.45776980040
    """
    h_p = liqair_h(0,0,1,wair,pres,entr=entr,airf=airf,temp=temp,dhum=dhum,
        dliq=dliq,chkvals=chkvals,chktol=chktol,airf0=airf0,temp0=temp0,
        dhum0=dhum0,dliq0=dliq0,chkbnd=chkbnd,mathargs=mathargs)
    dtot = h_p**(-1)
    return dtot

def kappa_s(wair,pres,entr=None,temp=None,airf=None,dhum=None,dliq=None,
    chkvals=False,chktol=_CHKTOL,airf0=None,temp0=None,dhum0=None,
    dliq0=None,chkbnd=False,mathargs=None):
    """Calculate wet air isentropic compressibility.
    
    Calculate the isentropic compressibility of wet air.
    
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
    :arg dliq: Liquid water density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dliq: float or None
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
    :arg dliq0: Initial guess for the liquid water density in kg/m3. If
        None (default) then `_approx_wep` is used.
    :type dliq0: float or None
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
    
    >>> kappa_s(0.5,1e5,entr=100.)
    9.43218607469e-6
    """
    airf, temp, dhum, dliq = eq_wpte(wair,pres,entr=entr,temp=temp,airf=airf,
        dhum=dhum,dliq=dliq,chkvals=chkvals,chktol=chktol,airf0=airf0,
        temp0=temp0,dhum0=dhum0,dliq0=dliq0,chkbnd=chkbnd,mathargs=mathargs)
    h_p = liqair_h(0,0,1,wair,pres,airf=airf,temp=temp,dhum=dhum,dliq=dliq)
    h_pp = liqair_h(0,0,2,wair,pres,airf=airf,temp=temp,dhum=dhum,dliq=dliq)
    kappa = -h_pp/h_p
    return kappa

def lapserate(wair,pres,entr=None,temp=None,airf=None,dhum=None,
    dliq=None,chkvals=False,chktol=_CHKTOL,airf0=None,temp0=None,
    dhum0=None,dliq0=None,chkbnd=False,mathargs=None):
    """Calculate wet air adiabatic lapse rate.
    
    Calculate the adiabatic lapse rate of wet air.
    
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
    :arg dliq: Liquid water density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dliq: float or None
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
    :arg dliq0: Initial guess for the liquid water density in kg/m3. If
        None (default) then `_approx_wep` is used.
    :type dliq0: float or None
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
    
    >>> lapserate(0.5,1e5,entr=100.)
    1.55067379031e-4
    """
    h_ep = liqair_h(0,1,1,wair,pres,entr=entr,airf=airf,temp=temp,dhum=dhum,
        dliq=dliq,chkvals=chkvals,chktol=chktol,airf0=airf0,temp0=temp0,
        dhum0=dhum0,dliq0=dliq0,chkbnd=chkbnd,mathargs=mathargs)
    gamma = h_ep
    return gamma

def temperature(wair,pres,entr=None,temp=None,airf=None,dhum=None,
    dliq=None,chkvals=False,chktol=_CHKTOL,airf0=None,temp0=None,
    dhum0=None,dliq0=None,chkbnd=False,mathargs=None):
    """Calculate wet air temperature.
    
    Calculate the temperature of wet air.
    
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
    :arg dliq: Liquid water density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dliq: float or None
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
    :arg dliq0: Initial guess for the liquid water density in kg/m3. If
        None (default) then `_approx_wep` is used.
    :type dliq0: float or None
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
    
    >>> temperature(0.5,1e5,entr=100.)
    280.393544899
    """
    h_e = liqair_h(0,1,0,wair,pres,entr=entr,airf=airf,temp=temp,dhum=dhum,
        dliq=dliq,chkvals=chkvals,chktol=chktol,airf0=airf0,temp0=temp0,
        dhum0=dhum0,dliq0=dliq0,chkbnd=chkbnd,mathargs=mathargs)
    temp = h_e
    return temp


## Potential thermodynamic functions
def _approx_pot(wair,temp,pres,ppot,airf,dhum,dliq):
    """Approximate ATDhDl2 at WTP1P2.
    
    Approximate the humid air dry fraction, temperature, humid air
    density, and liquid water density after adiabatic displacement.
    
    :arg float wair: Total dry fraction in kg/kg.
    :arg float temp: In-situ temperature in K.
    :arg float pres: In-situ pressure in Pa.
    :arg float ppot: Potential pressure in Pa.
    :arg float airf: In-situ humid air dry fraction in kg/kg (unused).
    :arg float dhum: In-situ humid air density in kg/m3 (unused).
    :arg float dliq: In-situ liquid water density in kg/m3 (unused).
    :returns: Potential humid air dry fraction, potential temperature,
        potential humid air density, and potential liquid water density
        (all in SI units).
    """
    pvsat0 = _PTPE*numpy.exp((_AVL+_BVL)*(1 - _TTP/temp)
        + _BVL*numpy.log(_TTP/temp))
    pvmax0 = pres * (1-wair)/(_EPSW*wair + 1-wair)
    if pvmax0 > pvsat0:
        # Parcel starts saturated
        pv0 = pvsat0
        a0 = (pres-pv0) / (pres-pv0 + _EPSW*pv0)
        ceff0 = (wair*_CDRY + (1-wair/a0)*_CLIQ + wair*(1-a0)/a0*_CVAP
            + wair*_RWAT*(1-a0)/a0 * pres/(pres-pv0)
                * (_AVL*_TTP/temp + _BVL*(_TTP/temp-1))**2)
        reff0 = (wair*_RDRY + wair*(1-a0)/a0*_RWAT
            + wair*(1-a0)/a0*_RWAT*pres/(pres-pv0)
                * (_AVL*_TTP/temp + _BVL*(_TTP/temp-1)))
        ginv0 = ceff0/reff0
        r = (_AVL+_BVL)/(ginv0+_BVL) - 1
        v = numpy.log((_TTP/temp)**ginv0 * pvmax0/_PTPE)/(ginv0+_BVL)
        if v <= r:
            # Parcel is unlikely to ever be unsaturated
            tpot = temp * (ppot/pres)**(reff0/ceff0)
            pv2 = _PTPE*numpy.exp((_AVL+_BVL)*(1 - _TTP/tpot)
                + _BVL*numpy.log(_TTP/tpot))
            apot = (ppot-pv2) / (ppot-pv2 + _EPSW*pv2)
        else:
            # Find where parcel de-saturates
            x = maths4.lamb2(v,r)
            ticl = _TTP/x
            picl = pres * (ticl/picl)**ginv0
            if ppot < picl:
                # Parcel ends saturated
                tpot = temp * (ppot/pres)**(reff0/ceff0)
                pv2 = _PTPE*numpy.exp((_AVL+_BVL)*(1 - _TTP/tpot)
                    + _BVL*numpy.log(_TTP/tpot))
                apot = (ppot-pv2) / (ppot-pv2 + _EPSW*pv2)
            else:
                #Parcel ends unsaturated
                p1 = picl
                t1 = ticl
                ceff1 = wair*_CDRY + (1-wair)*_CVAP
                reff1 = wair*_RDRY + (1-wair)*_RWAT
                tpot = t1 * (ppot/p1)**(reff1/ceff1)
                apot = wair
    else:
        # Parcel starts unsaturated
        ticl, __, __ = liqair4a._approx_icl(wair,temp,pres,dhum)
        ceff0 = wair*_CDRY + (1-wair)*_CVAP
        reff0 = wair*_RDRY + (1-wair)*_RWAT
        picl = pres * (ticl/temp)**(ceff0/reff0)
        if ppot < picl:
            # Parcel ends saturated
            t1 = ticl
            p1 = picl
            ceff1 = (wair*_CDRY + (1-wair)*_CVAP
                + (1-wair)*_RWAT*p1/(p1-pv1) * ((_AVL+_BVL)*_TTP/t1 - _BVL)**2)
            reff1 = (wair*_RDRY + (1-wair)*_RWAT
                + (1-wair)*_RWAT*p1/(p1-pv1) * ((_AVL+_BVL)*_TTP/t1 - _BVL))
            tpot = t1 * (ppot/p1)**(reff1/ceff1)
            pv2 = _PTPE*numpy.exp((_AVL+_BVL)*(1 - _TTP/tpot)
                + _BVL*numpy.log(_TTP/tpot))
            apot = (ppot-pv2) / (ppot-pv2 + _EPSW*pv2)
        else:
            # Parcel ends unsaturated
            tpot = temp * (ppot/pres)**(reff0/ceff0)
            apot = wair
    dhpot = ppot/(_RDRY*tpot) / (apot + (1-apot)/_EPSW)
    dlpot = _dliq_default(tpot,ppot)
    return apot, tpot, dhpot, dlpot

def _diff_pot(a2,t2,dh2,dl2,wair,temp,pres,ppot,airf,dhum,dliq):
    """Calculate wet air disequilibrium at WTP1P2.
    
    Calculate both sides of the equations
    
        given potential pressure = potential pressure of humid air
        given potential pressure = potential pressure of liquid water
        potential liquid water potential = potential vapour potential
        initial entropy of wet air = potential entropy of wet air
    
    and their Jacobians with respect to potential humid air dry
    fraction, temperature, humid air density, and liquid water density.
    Solving these equations gives equilibrium values at the given total
    dry air fraction, in-situ temperature, in-situ pressure, and
    potential pressure.
    
    :arg float a2: Potential humid air dry fraction in kg/kg.
    :arg float t2: Potential temperature in K.
    :arg float dh2: Potential humid air density in kg/m3.
    :arg float dl2: Potential liquid water density in kg/m3.
    :arg float wair: Total dry fraction in kg/kg.
    :arg float temp: In-situ temperature in K.
    :arg float pres: In-situ pressure in Pa.
    :arg float ppot: Potential pressure in Pa.
    :arg float airf: In-situ humid air dry fraction in kg/kg.
    :arg float dhum: In-situ humid air density in kg/m3.
    :arg float dliq: In-situ liquid water density in kg/m3.
    :returns: Left-hand side of the equations, right-hand side,
        Jacobian of LHS, and Jacobian of RHS.
    :rtype: tuple(array(float))
    """
    phpot = _air_eq_pressure(0,0,0,a2,t2,dh2)
    plpot = _flu_eq_pressure(0,0,t2,dl2)
    glpot = _flu_eq_chempot(0,0,t2,dl2)
    gvpot = _air_eq_vappot(0,0,0,a2,t2,dh2)
    sh1 = -_air_f(0,1,0,airf,temp,dhum)
    sl1 = -_flu_f(1,0,temp,dliq)
    entr = wair/airf*sh1 + (1-wair/airf)*sl1
    sh2 = -_air_f(0,1,0,a2,t2,dh2)
    sl2 = -_flu_f(1,0,t2,dl2)
    epot = wair/a2*sh2 + (1-wair/a2)*sl2
    lhs = numpy.array([ppot, ppot, glpot, entr])
    rhs = numpy.array([phpot, plpot, gvpot, epot])
    
    ph2_a = _air_eq_pressure(1,0,0,a2,t2,dh2)
    ph2_t = _air_eq_pressure(0,1,0,a2,t2,dh2)
    ph2_d = _air_eq_pressure(0,0,1,a2,t2,dh2)
    pl2_t = _flu_eq_pressure(1,0,t2,dl2)
    pl2_d = _flu_eq_pressure(0,1,t2,dl2)
    gl2_t = _flu_eq_chempot(1,0,t2,dl2)
    gl2_d = _flu_eq_chempot(0,1,t2,dl2)
    gv2_a = _air_eq_vappot(1,0,0,a2,t2,dh2)
    gv2_t = _air_eq_vappot(0,1,0,a2,t2,dh2)
    gv2_d = _air_eq_vappot(0,0,1,a2,t2,dh2)
    sh2_a = -_air_f(1,1,0,a2,t2,dh2)
    sh2_t = -_air_f(0,2,0,a2,t2,dh2)
    sh2_d = -_air_f(0,1,1,a2,t2,dh2)
    sl2_t = -_flu_f(2,0,t2,dl2)
    sl2_d = -_flu_f(1,1,t2,dl2)
    s2_a = -wair/a2**2 * (sh2 - a2*sh2_a - sl2)
    s2_t = wair/a2*sh2_t + (1-wair/a2)*sl2_t
    s2_dh = wair/a2*sh2_d
    s2_dl = (1-wair/a2)*sl2_d
    dlhs = numpy.array([[0.,0.,0.,0.], [0.,0.,0.,0.], [0.,gl2_t,0.,gl2_d],
        [0.,0.,0.,0.]])
    drhs = numpy.array([[ph2_a,ph2_t,ph2_d,0.], [0.,pl2_t,0.,pl2_d],
        [gv2_a,gv2_t,gv2_d,0.], [s2_a,s2_t,s2_dh,s2_dl]])
    return lhs, rhs, dlhs, drhs

def eq_pot(wair,temp,pres,ppot,airf=None,dhum=None,dliq=None,apot=None,
    tpot=None,dhpot=None,dlpot=None,chkvals=False,chktol=_CHKTOL,
    airf0=None,dhum0=None,dliq0=None,apot0=None,tpot0=None,dhpot0=None,
    dlpot0=None,chkbnd=False,mathargs=None):
    """Get primary variables at WTP1P2.
    
    Get the values of the equilibrium potential humid air dry fraction,
    temperature, humid air density, and liquid water density for the
    given total dry fraction, in-situ temperature, in-situ pressure, and
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
    :arg dliq: In-situ liquid water density in kg/m3. If unknown, pass
        None (default) and it will be calculated.
    :type dliq: float or None
    :arg apot: Potential humid air dry fraction in kg/kg. If unknown,
        pass None (default) and it will be calculated.
    :type apot: float or None
    :arg tpot: Potential temperature in K. If unknown, pass None
        (default) and it will be calculated.
    :type tpot: float or None
    :arg dhpot: Potential humid air density in kg/m3. If unknown, pass
        None (default) and it will be calculated.
    :type dhpot: float or None
    :arg dlpot: Potential liquid water density in kg/m3. If unknown,
        pass None (default) and it will be calculated.
    :type dlpot: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg airf0: Initial guess for the in-situ humid air dry fraction in
        kg/kg. If None (default) then `liqair4a._approx_tp` is used.
    :type airf0: float or None
    :arg dhum0: Initial guess for the in-situ humid air density in
        kg/m3. If None (default) then `liqair4a._approx_tp` is used.
    :type dhum0: float or None
    :arg dliq0: Initial guess for the in-situ liquid water density in
        kg/m3. If None (default) then `liqair4a._approx_tp` is used.
    :type dliq0: float or None
    :arg apot0: Initial guess for the potential humid air dry fraction
        in kg/kg. If None (default) then `_approx_pot` is used.
    :type apot0: float or None
    :arg temp0: Initial guess for the potential temperature in K. If
        None (default) then `_approx_pot` is used.
    :type temp0: float or None
    :arg dhpot0: Initial guess for the potential humid air density in
        kg/m3. If None (default) then `_approx_pot` is used.
    :type dhpot0: float or None
    :arg dlpot0: Initial guess for the potential liquid water density in
        kg/m3. If None (default) then `_approx_pot` is used.
    :type dlpot0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: In-situ humid air dry fraction, in-situ humid air density,
        in-situ liquid water density, potential humid air dry fraction,
        potential temperature, potential humid air density, and
        potential liquid water density (all in SI units).
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    """
    if any(val is None for val in (airf,dhum,dliq)):
        airf, __, __, dhum, dliq = liqair4a.eq_atpe(temp=temp,pres=pres,
            airf0=airf0,dhum0=dhum0,dliq0=dliq0,mathargs=mathargs)
    if any(val is None for val in (apot,tpot,dhpot,dlpot)):
        x0 = (apot,tpot,dhpot,dlpot)
        fargs = (wair,temp,pres,ppot,airf,dhum,dliq)
        if mathargs is None:
            mathargs = dict()
        x1 = _newton(_diff_pot,x0,_approx_pot,fargs=fargs,**mathargs)
        apot, tpot, dhpot, dlpot = x1
    
    _chkhumbnds(airf,temp,dhum,chkbnd=chkbnd)
    _chkhumbnds(apot,tpot,dhpot,chkbnd=chkbnd)
    _chkflubnds(tpot,dlpot,chkbnd=chkbnd)
    if not chkvals:
        return airf, dhum, dliq, apot, tpot, dhpot, dlpot
    
    lhs1, rhs1, __, __ = liqair4a._diff_tp(airf,dhum,dliq,temp,pres)
    lhs2, rhs2, __, __ = _diff_pot(apot,tpot,dhpot,dlpot,wair,temp,pres,ppot,
        airf,dhum,dliq)
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
    return airf, dhum, dliq, apot, tpot, dhpot, dlpot

def potdensity(wair,temp,pres,ppot,airf=None,dhum=None,dliq=None,
    apot=None,tpot=None,dhpot=None,dlpot=None,chkvals=False,
    chktol=_CHKTOL,airf0=None,dhum0=None,dliq0=None,apot0=None,
    tpot0=None,dhpot0=None,dlpot0=None,chkbnd=False,mathargs=None):
    """Calculate wet air potential density.
    
    Calculate the potential density of wet air.
    
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
    :arg dliq: In-situ liquid water density in kg/m3. If unknown, pass
        None (default) and it will be calculated.
    :type dliq: float or None
    :arg apot: Potential humid air dry fraction in kg/kg. If unknown,
        pass None (default) and it will be calculated.
    :type apot: float or None
    :arg tpot: Potential temperature in K. If unknown, pass None
        (default) and it will be calculated.
    :type tpot: float or None
    :arg dhpot: Potential humid air density in kg/m3. If unknown, pass
        None (default) and it will be calculated.
    :type dhpot: float or None
    :arg dlpot: Potential liquid water density in kg/m3. If unknown,
        pass None (default) and it will be calculated.
    :type dlpot: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg airf0: Initial guess for the in-situ humid air dry fraction in
        kg/kg. If None (default) then `liqair4a._approx_tp` is used.
    :type airf0: float or None
    :arg dhum0: Initial guess for the in-situ humid air density in
        kg/m3. If None (default) then `liqair4a._approx_tp` is used.
    :type dhum0: float or None
    :arg dliq0: Initial guess for the in-situ liquid water density in
        kg/m3. If None (default) then `liqair4a._approx_tp` is used.
    :type dliq0: float or None
    :arg apot0: Initial guess for the potential humid air dry fraction
        in kg/kg. If None (default) then `_approx_pot` is used.
    :type apot0: float or None
    :arg temp0: Initial guess for the potential temperature in K. If
        None (default) then `_approx_pot` is used.
    :type temp0: float or None
    :arg dhpot0: Initial guess for the potential humid air density in
        kg/m3. If None (default) then `_approx_pot` is used.
    :type dhpot0: float or None
    :arg dlpot0: Initial guess for the potential liquid water density in
        kg/m3. If None (default) then `_approx_pot` is used.
    :type dlpot0: float or None
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
    
    >>> potdensity(0.5,300.,1e4,1e5)
    1.22550664945
    """
    airf, dhum, dliq, apot, tpot, dhpot, dlpot = eq_pot(wair,temp,pres,ppot,
        airf=airf,dhum=dhum,dliq=dliq,apot=apot,tpot=tpot,dhpot=dhpot,
        dlpot=dlpot,chkvals=chkvals,chktol=chktol,airf0=airf0,dhum0=dhum0,
        dliq0=dliq0,apot0=apot0,tpot0=tpot0,dhpot0=dhpot0,dlpot0=dlpot0,
        chkbnd=chkbnd,mathargs=mathargs)
    hpot_p = liqair_h(0,0,1,wair,ppot,temp=tpot,airf=apot,dhum=dhpot,dliq=dlpot)
    dpot = hpot_p**(-1)
    return dpot

def potenthalpy(wair,temp,pres,ppot,airf=None,dhum=None,dliq=None,
    apot=None,tpot=None,dhpot=None,dlpot=None,chkvals=False,
    chktol=_CHKTOL,airf0=None,dhum0=None,dliq0=None,apot0=None,
    tpot0=None,dhpot0=None,dlpot0=None,chkbnd=False,mathargs=None):
    """Calculate wet air potential enthalpy.
    
    Calculate the potential enthalpy of wet air.
    
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
    :arg dliq: In-situ liquid water density in kg/m3. If unknown, pass
        None (default) and it will be calculated.
    :type dliq: float or None
    :arg apot: Potential humid air dry fraction in kg/kg. If unknown,
        pass None (default) and it will be calculated.
    :type apot: float or None
    :arg tpot: Potential temperature in K. If unknown, pass None
        (default) and it will be calculated.
    :type tpot: float or None
    :arg dhpot: Potential humid air density in kg/m3. If unknown, pass
        None (default) and it will be calculated.
    :type dhpot: float or None
    :arg dlpot: Potential liquid water density in kg/m3. If unknown,
        pass None (default) and it will be calculated.
    :type dlpot: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg airf0: Initial guess for the in-situ humid air dry fraction in
        kg/kg. If None (default) then `liqair4a._approx_tp` is used.
    :type airf0: float or None
    :arg dhum0: Initial guess for the in-situ humid air density in
        kg/m3. If None (default) then `liqair4a._approx_tp` is used.
    :type dhum0: float or None
    :arg dliq0: Initial guess for the in-situ liquid water density in
        kg/m3. If None (default) then `liqair4a._approx_tp` is used.
    :type dliq0: float or None
    :arg apot0: Initial guess for the potential humid air dry fraction
        in kg/kg. If None (default) then `_approx_pot` is used.
    :type apot0: float or None
    :arg temp0: Initial guess for the potential temperature in K. If
        None (default) then `_approx_pot` is used.
    :type temp0: float or None
    :arg dhpot0: Initial guess for the potential humid air density in
        kg/m3. If None (default) then `_approx_pot` is used.
    :type dhpot0: float or None
    :arg dlpot0: Initial guess for the potential liquid water density in
        kg/m3. If None (default) then `_approx_pot` is used.
    :type dlpot0: float or None
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
    
    >>> potenthalpy(0.5,300.,1e4,1e5)
    655155.797982
    """
    airf, dhum, dliq, apot, tpot, dhpot, dlpot = eq_pot(wair,temp,pres,ppot,
        airf=airf,dhum=dhum,dliq=dliq,apot=apot,tpot=tpot,dhpot=dhpot,
        dlpot=dlpot,chkvals=chkvals,chktol=chktol,airf0=airf0,dhum0=dhum0,
        dliq0=dliq0,apot0=apot0,tpot0=tpot0,dhpot0=dhpot0,dlpot0=dlpot0,
        chkbnd=chkbnd,mathargs=mathargs)
    hpot = liqair_h(0,0,0,wair,ppot,temp=tpot,airf=apot,dhum=dhpot,dliq=dlpot)
    return hpot

def pottemp(wair,temp,pres,ppot,airf=None,dhum=None,dliq=None,apot=None,
    tpot=None,dhpot=None,dlpot=None,chkvals=False,chktol=_CHKTOL,
    airf0=None,dhum0=None,dliq0=None,apot0=None,tpot0=None,dhpot0=None,
    dlpot0=None,chkbnd=False,mathargs=None):
    """Calculate wet air potential temperature.
    
    Calculate the potential temperature of wet air.
    
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
    :arg dliq: In-situ liquid water density in kg/m3. If unknown, pass
        None (default) and it will be calculated.
    :type dliq: float or None
    :arg apot: Potential humid air dry fraction in kg/kg. If unknown,
        pass None (default) and it will be calculated.
    :type apot: float or None
    :arg tpot: Potential temperature in K. If unknown, pass None
        (default) and it will be calculated.
    :type tpot: float or None
    :arg dhpot: Potential humid air density in kg/m3. If unknown, pass
        None (default) and it will be calculated.
    :type dhpot: float or None
    :arg dlpot: Potential liquid water density in kg/m3. If unknown,
        pass None (default) and it will be calculated.
    :type dlpot: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg airf0: Initial guess for the in-situ humid air dry fraction in
        kg/kg. If None (default) then `liqair4a._approx_tp` is used.
    :type airf0: float or None
    :arg dhum0: Initial guess for the in-situ humid air density in
        kg/m3. If None (default) then `liqair4a._approx_tp` is used.
    :type dhum0: float or None
    :arg dliq0: Initial guess for the in-situ liquid water density in
        kg/m3. If None (default) then `liqair4a._approx_tp` is used.
    :type dliq0: float or None
    :arg apot0: Initial guess for the potential humid air dry fraction
        in kg/kg. If None (default) then `_approx_pot` is used.
    :type apot0: float or None
    :arg temp0: Initial guess for the potential temperature in K. If
        None (default) then `_approx_pot` is used.
    :type temp0: float or None
    :arg dhpot0: Initial guess for the potential humid air density in
        kg/m3. If None (default) then `_approx_pot` is used.
    :type dhpot0: float or None
    :arg dlpot0: Initial guess for the potential liquid water density in
        kg/m3. If None (default) then `_approx_pot` is used.
    :type dlpot0: float or None
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
    
    >>> pottemp(0.5,300.,1e4,1e5)
    348.222379217
    """
    airf, dhum, dliq, apot, tpot, dhpot, dlpot = eq_pot(wair,temp,pres,ppot,
        airf=airf,dhum=dhum,dliq=dliq,apot=apot,tpot=tpot,dhpot=dhpot,
        dlpot=dlpot,chkvals=chkvals,chktol=chktol,airf0=airf0,dhum0=dhum0,
        dliq0=dliq0,apot0=apot0,tpot0=tpot0,dhpot0=dhpot0,dlpot0=dlpot0,
        chkbnd=chkbnd,mathargs=mathargs)
    return tpot

