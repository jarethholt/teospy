"""Seawater enthalpy and potential temperature.

This module provides the enthalpy of seawater (liquid water and salt) as
a function of salinity, entropy, and pressure. It also provides the
potential temperature and potential enthalpy and the related expansion
and contraction coefficients.

:Examples:

>>> sea_h(0,0,0,0.035,500.,1e5)
145481.970750
>>> sea_h(0,0,1,0.035,500.,1e5)
9.81092930969e-4
>>> sea_h(1,0,1,0.035,500.,1e5)
-6.84629317367e-4
>>> sea_h(0,2,0,0.035,500.,1e5)
7.72873234085e-2
>>> contraction_t(0.035,500.,1e5)
0.725209049049
>>> pottemp(0.035,1e7,1e5,temp=300.)
299.771869405
>>> potenthalpy(0.035,1e7,1e5,temp=300.)
106307.996083
>>> contraction_h(0.035,1e7,1e5,entr=500.)
0.697779873590
>>> contraction_theta(0.035,1e7,1e5,entr=500.)
0.717342103505

:Functions:

* :func:`eq_sep`: Calculate equilibrium quantities at salinity, entropy,
  and pressure.
* :func:`sea_h`: Seawater enthalpy with derivatives.
* :func:`temperature`: Seawater temperature.
* :func:`contraction_t`: Seawater haline contraction coefficient at
  constant in-situ temperature.
* :func:`expansion_t`: Seawater expansion coefficient with respect to
  in-situ temperature.
* :func:`eq_pot`: Calculate equilibrium quantities at salinity, in-situ
  pressure, potential pressure, and either entropy or in-situ
  temperature.
* :func:`pottemp`: Seawater potential temperature.
* :func:`potdensity`: Seawater potential density.
* :func:`potenthalpy`: Seawater potential enthalpy.
* :func:`contraction_h`: Seawater haline contraction coefficient at
  constant potential enthalpy.
* :func:`contraction_theta`: Seawater haline contraction coefficient at
  constant potential temperature (entropy).
* :func:`expansion_h`: Seawater expansion coefficient with respect to
  potential enthalpy.
* :func:`expansion_theta`: Seawater expansion coefficient with respect
  to potential temperature.

"""

__all__ = ['eq_sep','sea_h','temperature','contraction_t','expansion_t',
    'eq_pot','pottemp','potdensity','potenthalpy','contraction_h',
    'contraction_theta','expansion_h','expansion_theta']

import numpy
import constants0
import flu1
import sal1
import flu2
import sal2
import maths3
import flu3a
import sea3a

_CHKTOL = constants0.CHKTOL
_RUNIV = constants0.RUNIV
_MSAL = constants0.MSAL
_TTP = constants0.TTP
_DLTP = constants0.DLTP
_CLIQ = constants0.CLIQ
_chkflubnds = constants0.chkflubnds
_chksalbnds = constants0.chksalbnds
_flu_f = flu1.flu_f
_flu_eq_pressure = flu2.eq_pressure
_flu_eq_entropy = flu2.eq_entropy
_sal_g = sal2.sal_g
_sal_eq_entropy = sal2.eq_entropy
_newton = maths3.newton


## Equilibrium functions for enthalpy
def _approx_sep(salt,entr,pres):
    """Approximate seawater T/DL at SEP.
    
    Approximate the temperature and liquid water density for seawater at
    the given salinity, entropy, and pressure. Uses a first-order Gibbs
    energy of salt and the default liquid density approximation from
    flu3a.
    
    :arg float salt: Salinity in kg/kg.
    :arg float entr: Entropy in J/kg/K.
    :arg float pres: Pressure in Pa.
    :returns: Temperature in K and liquid water density in kg/m3.
    """
    if salt == 0:
        s0 = 0.
    else:
        s0 = -_RUNIV/_MSAL * salt*numpy.log(salt)
    temp = _TTP*numpy.exp((entr - s0)/_CLIQ)
    dliq = flu3a._dliq_default(temp,pres)
    return temp, dliq

def _diff_sep(t,d,salt,entr,pres,useext=False):
    """Calculate seawater disequilibrium at SEP.
    
    Calculate both sides of the equations
    
        given entropy = entropy of seawater
        given pressure = pressure of liquid water
    
    and their Jacobians with respect to temperature and liquid water
    density. Solving these equations gives the equilibrium temperature
    and fluid water density for the given salinity, entropy, and
    pressure.
    
    :arg float t: Temperature in K.
    :arg float d: Liquid water density in kg/m3.
    :arg float salt: Salinity in kg/kg.
    :arg float entr: Entropy in J/kg/K.
    :arg float pres: Pressure in Pa.
    :arg bool useext: If False (default) then the salt contribution is
        calculated from _GSCOEFFS in sal1; if True, from _GSCOEFFS_EXT.
    :returns: Left-hand side of the equation, right-hand side,
        Jacobian of LHS, and Jacobian of RHS.
    :rtype: tuple(array(float))
    """
    sl = _flu_eq_entropy(0,0,t,d)
    ss = _sal_eq_entropy(0,0,0,salt,t,pres,useext=useext)
    ssea = sl + ss
    pliq = _flu_eq_pressure(0,0,t,d)
    lhs = numpy.array([entr, pres])
    rhs = numpy.array([ssea, pliq])
    
    sl_t = _flu_eq_entropy(1,0,t,d)
    sl_d = _flu_eq_entropy(0,1,t,d)
    ss_t = _sal_eq_entropy(0,1,0,salt,t,pres,useext=useext)
    ssea_t = sl_t + ss_t
    ssea_d = sl_d
    pliq_t = _flu_eq_pressure(1,0,t,d)
    pliq_d = _flu_eq_pressure(0,1,t,d)
    dlhs = numpy.array([[0.,0.], [0.,0.]])
    drhs = numpy.array([[ssea_t,ssea_d], [pliq_t,pliq_d]])
    return lhs, rhs, dlhs, drhs

def eq_sep(salt,entr,pres,temp=None,dliq=None,chkvals=False,
    chktol=_CHKTOL,temp0=None,dliq0=None,chkbnd=False,useext=False,
    mathargs=None):
    """Get primary seawater variables at SEP.
    
    Get the values of the equilibrium temperature and liquid water
    density for seawater at the given salinity, entropy, and pressure.
    
    If the calculation has already been done, the results can be passed
    to avoid unnecessary repeat calculations. If enough values are
    passed, they will be checked for consistency if chkvals is True.
    
    :arg float salt: Salinity in kg/kg.
    :arg float entr: Entropy in J/kg/K.
    :arg float pres: Pressure in Pa.
    :arg temp: Temperature in K. If unknown, pass None (default) and it
        will be calculated.
    :type dliq: float or None
    :arg dliq: Liquid water density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dliq: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg temp0: Initial guess for the temperature in K. If None
        (default) then `_approx_sep` is used.
    :type temp0: float or None
    :arg dliq0: Initial guess for the liquid water density in kg/m3. If
        None (default) then `_approx_sep` is used.
    :type dliq0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Temperature and liquid water density (in SI units).
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    """
    if any(val is None for val in (temp,dliq)):
        x0 = (temp0,dliq0)
        fargs = (salt,entr,pres)
        fkwargs = {'useext': useext}
        if mathargs is None:
            mathargs = dict()
        x1 = _newton(_diff_sep,x0,_approx_sep,fargs=fargs,fkwargs=fkwargs,
            **mathargs)
        temp, dliq = x1
    
    _chkflubnds(temp,dliq,chkbnd=chkbnd)
    _chksalbnds(salt,temp,pres,chkbnd=chkbnd)
    if not chkvals:
        return temp, dliq
    
    lhs, rhs, __, __ = _diff_sep(temp,dliq,salt,entr,pres,useext=useext)
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
    return temp, dliq


## Enthalpy functions
def sea_h(drvs,drve,drvp,salt,entr,pres,temp=None,dliq=None,
    chkvals=False,chktol=_CHKTOL,temp0=None,dliq0=None,chkbnd=False,
    useext=False,mathargs=None):
    """Calculate seawater enthalpy with derivatives.
    
    Calculate the specific enthalpy of seawater or its derivatives with
    respect to salinity, entropy, and pressure.
    
    :arg int drvs: Number of salinity derivatives.
    :arg int drve: Number of entropy derivatives.
    :arg int drvp: Number of pressure derivatives.
    :arg float salt: Salinity in kg/kg.
    :arg float entr: Entropy in J/kg/K.
    :arg float pres: Pressure in Pa.
    :arg temp: Temperature in K. If unknown, pass None (default) and it
        will be calculated.
    :type dliq: float or None
    :arg dliq: Liquid water density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dliq: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg temp0: Initial guess for the temperature in K. If None
        (default) then `_approx_sep` is used.
    :type temp0: float or None
    :arg dliq0: Initial guess for the liquid water density in kg/m3. If
        None (default) then `_approx_sep` is used.
    :type dliq0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Enthalpy in units of
        (J/kg) / (kg/kg)^drvs / (J/kg/K)^drve / Pa^drvp.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> sea_h(0,0,0,0.035,500.,1e5)
    145481.970750
    >>> sea_h(0,0,1,0.035,500.,1e5)
    9.81092930969e-4
    >>> sea_h(1,0,1,0.035,500.,1e5)
    -6.84629317367e-4
    >>> sea_h(0,2,0,0.035,500.,1e5)
    7.72873234085e-2
    """
    temp, dliq = eq_sep(salt,entr,pres,temp=temp,dliq=dliq,chkvals=chkvals,
        chktol=chktol,temp0=temp0,dliq0=dliq0,chkbnd=chkbnd,useext=useext,
        mathargs=mathargs)
    
    # Simple cases: First-order derivatives
    if (drvs,drve,drvp) == (0,0,0):
        fl = _flu_f(0,0,temp,dliq)
        fl_t = _flu_f(1,0,temp,dliq)
        fl_d = _flu_f(0,1,temp,dliq)
        gs = _sal_g(0,0,0,salt,temp,pres,useext=useext)
        gs_t = _sal_g(0,1,0,salt,temp,pres,useext=useext)
        hl = fl - temp*fl_t + dliq*fl_d
        hs = gs - temp*gs_t
        h = hl + hs
        return h
    elif (drvs,drve,drvp) == (1,0,0):
        gs_s = _sal_g(1,0,0,salt,temp,pres,useext=useext)
        h_s = gs_s
        return h_s
    elif (drvs,drve,drvp) == (0,1,0):
        h_e = temp
        return h_e
    elif (drvs,drve,drvp) == (0,0,1):
        gs_p = _sal_g(0,0,1,salt,temp,pres,useext=useext)
        hl_p = dliq**(-1)
        hs_p = gs_p
        h_p = hl_p + hs_p
        return h_p
    
    # Higher-order derivatives require the Jacobian
    __, __, __, ep_td = _diff_sep(temp,dliq,salt,entr,pres,useext=useext)
    if (drvs,drve,drvp) == (2,0,0):
        gs_ss = _sal_g(2,0,0,salt,temp,pres,useext=useext)
        gs_st = _sal_g(1,1,0,salt,temp,pres,useext=useext)
        e_s = -gs_st
        p_s = 0.
        ep_s = numpy.array([e_s,p_s])
        h_sx = numpy.array([gs_st,0.])
        td_s = numpy.linalg.solve(ep_td,-ep_s)
        h_ss = gs_ss + h_sx.dot(td_s)
        return h_ss
    elif (drvs,drve,drvp) == (1,1,0):
        gs_st = _sal_g(1,1,0,salt,temp,pres,useext=useext)
        ep_e = numpy.array([1.,0.])
        h_sx = numpy.array([gs_st,0.])
        td_e = numpy.linalg.solve(ep_td,ep_e)
        h_se = h_sx.dot(td_e)
        return h_se
    elif (drvs,drve,drvp) == (1,0,1):
        gs_st = _sal_g(1,1,0,salt,temp,pres,useext=useext)
        gs_sp = _sal_g(1,0,1,salt,temp,pres,useext=useext)
        gs_tp = _sal_g(0,1,1,salt,temp,pres,useext=useext)
        ep_p = numpy.array([gs_tp,1.])
        h_sx = numpy.array([gs_st,0.])
        td_p = numpy.linalg.solve(ep_td,ep_p)
        h_sp = gs_sp + h_sx.dot(td_p)
        return h_sp
    elif (drvs,drve,drvp) == (0,2,0):
        ep_e = numpy.array([1.,0.])
        t_e, d_e = numpy.linalg.solve(ep_td,ep_e)
        h_ee = t_e
        return h_ee
    elif (drvs,drve,drvp) == (0,1,1):
        gs_tp = _sal_g(0,1,1,salt,temp,pres,useext=useext)
        ep_p = numpy.array([gs_tp,1.])
        t_p, d_p = numpy.linalg.solve(ep_td,ep_p)
        h_ep = t_p
        return h_ep
    elif (drvs,drve,drvp) == (0,0,2):
        gs_tp = _sal_g(0,1,1,salt,temp,pres,useext=useext)
        gs_pp = _sal_g(0,0,2,salt,temp,pres,useext=useext)
        ep_p = numpy.array([gs_tp,1.])
        h_px = numpy.array([gs_tp,-dliq**(-2)])
        td_p = numpy.linalg.solve(ep_td,ep_p)
        h_pp = gs_pp + h_px.dot(td_p)
        return h_pp
    
    # Should not have made it this far!
    errmsg = 'Derivatives {0} not recognized'.format((drvs,drve,drvp))
    raise ValueError(errmsg)

def temperature(salt,entr,pres,temp=None,dliq=None,chkvals=False,
    chktol=_CHKTOL,temp0=None,dliq0=None,chkbnd=False,useext=False,
    mathargs=None):
    """Calculate seawater temperature.
    
    Calculate the temperature of seawater from salinity, entropy, and
    pressure.
    
    :arg float salt: Salinity in kg/kg.
    :arg float entr: Entropy in J/kg/K.
    :arg float pres: Pressure in Pa.
    :arg temp: Temperature in K. If unknown, pass None (default) and it
        will be calculated.
    :type dliq: float or None
    :arg dliq: Liquid water density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dliq: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg temp0: Initial guess for the temperature in K. If None
        (default) then `_approx_sep` is used.
    :type temp0: float or None
    :arg dliq0: Initial guess for the liquid water density in kg/m3. If
        None (default) then `_approx_sep` is used.
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
    
    >>> temperature(0.035,500.,1e5)
    309.557955853
    """
    temp, dliq = eq_sep(salt,entr,pres,temp=temp,dliq=dliq,chkvals=chkvals,
        chktol=chktol,temp0=temp0,dliq0=dliq0,chkbnd=chkbnd,useext=useext,
        mathargs=mathargs)
    return temp

def contraction_t(salt,entr,pres,temp=None,dliq=None,chkvals=False,
    chktol=_CHKTOL,temp0=None,dliq0=None,chkbnd=False,useext=False,
    mathargs=None):
    """Calculate seawater haline contraction coefficient.
    
    Calculate the haline contraction coefficient of seawater at constant
    in-situ temperature from salinity, entropy, and pressure.
    
    :arg float salt: Salinity in kg/kg.
    :arg float entr: Entropy in J/kg/K.
    :arg float pres: Pressure in Pa.
    :arg temp: Temperature in K. If unknown, pass None (default) and it
        will be calculated.
    :type dliq: float or None
    :arg dliq: Liquid water density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dliq: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg temp0: Initial guess for the temperature in K. If None
        (default) then `_approx_sep` is used.
    :type temp0: float or None
    :arg dliq0: Initial guess for the liquid water density in kg/m3. If
        None (default) then `_approx_sep` is used.
    :type dliq0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Contraction coefficient in 1/(kg/kg).
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> contraction_t(0.035,500.,1e5)
    0.725209049049
    """
    temp, dliq = eq_sep(salt,entr,pres,temp=temp,dliq=dliq,chkvals=chkvals,
        chktol=chktol,temp0=temp0,dliq0=dliq0,chkbnd=chkbnd,useext=useext,
        mathargs=mathargs)
    h_p = sea_h(0,0,1,salt,entr,pres,temp=temp,dliq=dliq,useext=useext)
    h_se = sea_h(1,1,0,salt,entr,pres,temp=temp,dliq=dliq,useext=useext)
    h_sp = sea_h(1,0,1,salt,entr,pres,temp=temp,dliq=dliq,useext=useext)
    h_ee = sea_h(0,2,0,salt,entr,pres,temp=temp,dliq=dliq,useext=useext)
    h_ep = sea_h(0,1,1,salt,entr,pres,temp=temp,dliq=dliq,useext=useext)
    beta = (h_se*h_ep/h_ee - h_sp) / h_p
    return beta

def expansion_t(salt,entr,pres,temp=None,dliq=None,chkvals=False,
    chktol=_CHKTOL,temp0=None,dliq0=None,chkbnd=False,useext=False,
    mathargs=None):
    """Calculate seawater thermal expansion coefficient.
    
    Calculate the thermal expansion coefficient of seawater with respect
    to in-situ temperature from salinity, entropy, and pressure.
    
    :arg float salt: Salinity in kg/kg.
    :arg float entr: Entropy in J/kg/K.
    :arg float pres: Pressure in Pa.
    :arg temp: Temperature in K. If unknown, pass None (default) and it
        will be calculated.
    :type dliq: float or None
    :arg dliq: Liquid water density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dliq: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg temp0: Initial guess for the temperature in K. If None
        (default) then `_approx_sep` is used.
    :type temp0: float or None
    :arg dliq0: Initial guess for the liquid water density in kg/m3. If
        None (default) then `_approx_sep` is used.
    :type dliq0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Expansion coefficient in 1/K.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> expansion_t(0.035,500.,1e5)
    3.77581809091e-4
    """
    temp, dliq = eq_sep(salt,entr,pres,temp=temp,dliq=dliq,chkvals=chkvals,
        chktol=chktol,temp0=temp0,dliq0=dliq0,chkbnd=chkbnd,useext=useext,
        mathargs=mathargs)
    h_p = sea_h(0,0,1,salt,entr,pres,temp=temp,dliq=dliq,useext=useext)
    h_ee = sea_h(0,2,0,salt,entr,pres,temp=temp,dliq=dliq,useext=useext)
    h_ep = sea_h(0,1,1,salt,entr,pres,temp=temp,dliq=dliq,useext=useext)
    alpha = h_ep / (h_p * h_ee)
    return alpha


## Potential temperature functions
def eq_pot(salt,pres,ppot,entr=None,temp=None,dliq=None,tpot=None,
    dlpot=None,chkvals=False,chktol=_CHKTOL,temp0=None,dliq0=None,
    tpot0=None,dlpot0=None,chkbnd=False,useext=False,mathargs=None):
    """Get primary values at SP1P2 and E1/T1.
    
    Get the values of the equilibrium in-situ liquid water density,
    potential temperature, and potential liquid density for the given
    salinity, in-situ pressure, potential pressure, and either the
    entropy or in-situ temperature.
    
    If the calculation has already been done, the results can be passed
    to avoid unnecessary repeat calculations. If enough values are
    passed, they will be checked for consistency if chkvals is True.
    
    :arg float salt: Salinity in kg/kg.
    :arg float pres: In-situ pressure in Pa.
    :arg float ppot: Potential pressure in Pa.
    :arg entr: Entropy in J/kg/K. Either this or `temp` must be
        provided.
    :type entr: float or None
    :arg temp: In-situ temperature in K. If `entr` is given, pass None
        (default) and it will be calculated.
    :type temp: float or None
    :arg dliq: In-situ liquid water density in kg/m3. If unknown, pass
        None (default) and it will be calculated.
    :type dliq: float or None
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
        (default) then `_approx_sep` is used.
    :type temp0: float or None
    :arg dliq0: Initial guess for the in-situ liquid water density in
        kg/m3. If None (default) then `_approx_sep` or
        `flu3a._dliq_default` are used.
    :type dliq0: float or None
    :arg tpot0: Initial guess for the potential temperature in K. If
        None (default) then `_approx_sep` is used.
    :type tpot0: float or None
    :arg dlpot0: Initial guess for the potential liquid water density in
        kg/m3. If None (default) then `_approx_sep` is used.
    :type dlpot0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Entropy, in-situ temperature, in-situ density, potential
        temperature, and potential density (all in SI units).
    :raises ValueError: If both `entr` and `temp` are None.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    """
    if temp is None:
        if entr is None:
            errmsg = 'One of entr or temp must be provided'
            raise ValueError(errmsg)
        temp, dliq = eq_sep(salt,entr,pres,temp0=temp0,dliq0=dliq0,
            chkbnd=chkbnd,useext=useext,mathargs=mathargs)
    else:
        if dliq is None:
            dliq = flu3a.eq_tp_liq(temp,pres,dliq0=dliq0,chkbnd=chkbnd,
                mathargs=mathargs)
        if entr is None:
            entr = sea3a.entropy(salt,temp,pres,dliq=dliq)
    
    if any(val is None for val in (tpot,dlpot)):
        tpot, dlpot = eq_sep(salt,entr,ppot,temp0=tpot0,dliq0=dlpot0,
            chkbnd=chkbnd,useext=useext,mathargs=mathargs)
    
    if not chkvals:
        return entr, temp, dliq, tpot, dlpot
    
    # Check for consistency
    lhs1, rhs1, __, __ = _diff_sep(temp,dliq,salt,entr,pres,useext=useext)
    lhs2, rhs2, __, __ = _diff_sep(tpot,dlpot,salt,entr,ppot,useext=useext)
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
    return entr, temp, dliq, tpot, dlpot

def pottemp(salt,pres,ppot,entr=None,temp=None,dliq=None,tpot=None,
    dlpot=None,chkvals=False,chktol=_CHKTOL,temp0=None,dliq0=None,
    tpot0=None,dlpot0=None,chkbnd=False,useext=False,mathargs=None):
    """Calculate seawater potential temperature.
    
    Calculate the potential temperature of seawater from the salinity,
    in-situ pressure, potential pressure, and entropy or in-situ
    temperature.
    
    :arg float salt: Salinity in kg/kg.
    :arg float pres: In-situ pressure in Pa.
    :arg float ppot: Potential pressure in Pa.
    :arg entr: Entropy in J/kg/K. Either this or `temp` must be
        provided.
    :type entr: float or None
    :arg temp: In-situ temperature in K. If `entr` is given, pass None
        (default) and it will be calculated.
    :type temp: float or None
    :arg dliq: In-situ liquid water density in kg/m3. If unknown, pass
        None (default) and it will be calculated.
    :type dliq: float or None
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
        (default) then `_approx_sep` is used.
    :type temp0: float or None
    :arg dliq0: Initial guess for the in-situ liquid water density in
        kg/m3. If None (default) then `_approx_sep` or
        `flu3a._dliq_default` are used.
    :type dliq0: float or None
    :arg tpot0: Initial guess for the potential temperature in K. If
        None (default) then `_approx_sep` is used.
    :type tpot0: float or None
    :arg dlpot0: Initial guess for the potential liquid water density in
        kg/m3. If None (default) then `_approx_sep` is used.
    :type dlpot0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Potential temperature in K.
    :raises ValueError: If both `entr` and `temp` are None.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> pottemp(0.035,1e7,1e5,temp=300.)
    299.771869405
    """
    entr, temp, dliq, tpot, dlpot = eq_pot(salt,pres,ppot,entr=entr,temp=temp,
        dliq=dliq,tpot=tpot,dlpot=dlpot,chkvals=chkvals,chktol=chktol,
        temp0=temp0,dliq0=dliq0,tpot0=tpot0,dlpot0=dlpot0,chkbnd=chkbnd,
        useext=useext,mathargs=mathargs)
    return tpot

def potdensity(salt,pres,ppot,entr=None,temp=None,dliq=None,tpot=None,
    dlpot=None,chkvals=False,chktol=_CHKTOL,temp0=None,dliq0=None,
    tpot0=None,dlpot0=None,chkbnd=False,useext=False,mathargs=None):
    """Calculate seawater potential density.
    
    Calculate the potential density of seawater from the salinity,
    in-situ pressure, potential pressure, and entropy or in-situ
    temperature.
    
    :arg float salt: Salinity in kg/kg.
    :arg float pres: In-situ pressure in Pa.
    :arg float ppot: Potential pressure in Pa.
    :arg entr: Entropy in J/kg/K. Either this or `temp` must be
        provided.
    :type entr: float or None
    :arg temp: In-situ temperature in K. If `entr` is given, pass None
        (default) and it will be calculated.
    :type temp: float or None
    :arg dliq: In-situ liquid water density in kg/m3. If unknown, pass
        None (default) and it will be calculated.
    :type dliq: float or None
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
        (default) then `_approx_sep` is used.
    :type temp0: float or None
    :arg dliq0: Initial guess for the in-situ liquid water density in
        kg/m3. If None (default) then `_approx_sep` or
        `flu3a._dliq_default` are used.
    :type dliq0: float or None
    :arg tpot0: Initial guess for the potential temperature in K. If
        None (default) then `_approx_sep` is used.
    :type tpot0: float or None
    :arg dlpot0: Initial guess for the potential liquid water density in
        kg/m3. If None (default) then `_approx_sep` is used.
    :type dlpot0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Potential density in kg/m3.
    :raises ValueError: If both `entr` and `temp` are None.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> potdensity(0.035,1e7,1e5,temp=300.)
    1022.71520130
    """
    entr, temp, dliq, tpot, dlpot = eq_pot(salt,pres,ppot,entr=entr,temp=temp,
        dliq=dliq,tpot=tpot,dlpot=dlpot,chkvals=chkvals,chktol=chktol,
        temp0=temp0,dliq0=dliq0,tpot0=tpot0,dlpot0=dlpot0,chkbnd=chkbnd,
        useext=useext,mathargs=mathargs)
    dpot = sea3a.density(salt,tpot,ppot,dliq=dlpot,useext=useext)
    return dpot

def potenthalpy(salt,pres,ppot,entr=None,temp=None,dliq=None,tpot=None,
    dlpot=None,chkvals=False,chktol=_CHKTOL,temp0=None,dliq0=None,
    tpot0=None,dlpot0=None,chkbnd=False,useext=False,mathargs=None):
    """Calculate seawater potential enthalpy.
    
    Calculate the potential enthalpy of seawater from the salinity,
    in-situ pressure, potential pressure, and entropy or in-situ
    temperature.
    
    :arg float salt: Salinity in kg/kg.
    :arg float pres: In-situ pressure in Pa.
    :arg float ppot: Potential pressure in Pa.
    :arg entr: Entropy in J/kg/K. Either this or `temp` must be
        provided.
    :type entr: float or None
    :arg temp: In-situ temperature in K. If `entr` is given, pass None
        (default) and it will be calculated.
    :type temp: float or None
    :arg dliq: In-situ liquid water density in kg/m3. If unknown, pass
        None (default) and it will be calculated.
    :type dliq: float or None
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
        (default) then `_approx_sep` is used.
    :type temp0: float or None
    :arg dliq0: Initial guess for the in-situ liquid water density in
        kg/m3. If None (default) then `_approx_sep` or
        `flu3a._dliq_default` are used.
    :type dliq0: float or None
    :arg tpot0: Initial guess for the potential temperature in K. If
        None (default) then `_approx_sep` is used.
    :type tpot0: float or None
    :arg dlpot0: Initial guess for the potential liquid water density in
        kg/m3. If None (default) then `_approx_sep` is used.
    :type dlpot0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Potential enthalpy in J/kg.
    :raises ValueError: If both `entr` and `temp` are None.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> potenthalpy(0.035,1e7,1e5,temp=300.)
    106307.996083
    """
    entr, temp, dliq, tpot, dlpot = eq_pot(salt,pres,ppot,entr=entr,temp=temp,
        dliq=dliq,tpot=tpot,dlpot=dlpot,chkvals=chkvals,chktol=chktol,
        temp0=temp0,dliq0=dliq0,tpot0=tpot0,dlpot0=dlpot0,chkbnd=chkbnd,
        useext=useext,mathargs=mathargs)
    hpot = sea3a.enthalpy(salt,tpot,ppot,dliq=dlpot,useext=useext)
    return hpot


## Potential expansion and contraction coefficients
def contraction_h(salt,pres,ppot,entr=None,temp=None,dliq=None,
    tpot=None,dlpot=None,chkvals=False,chktol=_CHKTOL,temp0=None,
    dliq0=None,tpot0=None,dlpot0=None,chkbnd=False,useext=False,
    mathargs=None):
    """Calculate seawater isenthalpic haline contraction coefficient.
    
    Calculate the haline contraction coefficient of seawater at constant
    potential enthalpy from salinity, in-situ pressure, potential
    pressure, and entropy or in-situ temperature.
    
    :arg float salt: Salinity in kg/kg.
    :arg float pres: In-situ pressure in Pa.
    :arg float ppot: Potential pressure in Pa.
    :arg entr: Entropy in J/kg/K. Either this or `temp` must be
        provided.
    :type entr: float or None
    :arg temp: In-situ temperature in K. If `entr` is given, pass None
        (default) and it will be calculated.
    :type temp: float or None
    :arg dliq: In-situ liquid water density in kg/m3. If unknown, pass
        None (default) and it will be calculated.
    :type dliq: float or None
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
        (default) then `_approx_sep` is used.
    :type temp0: float or None
    :arg dliq0: Initial guess for the in-situ liquid water density in
        kg/m3. If None (default) then `_approx_sep` or
        `flu3a._dliq_default` are used.
    :type dliq0: float or None
    :arg tpot0: Initial guess for the potential temperature in K. If
        None (default) then `_approx_sep` is used.
    :type tpot0: float or None
    :arg dlpot0: Initial guess for the potential liquid water density in
        kg/m3. If None (default) then `_approx_sep` is used.
    :type dlpot0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Contraction coefficient in 1/(kg/kg).
    :raises ValueError: If both `entr` and `temp` are None.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> contraction_h(0.035,1e7,1e5,entr=500.)
    0.697779873590
    """
    entr, temp, dliq, tpot, dlpot = eq_pot(salt,pres,ppot,entr=entr,temp=temp,
        dliq=dliq,tpot=tpot,dlpot=dlpot,chkvals=chkvals,chktol=chktol,
        temp0=temp0,dliq0=dliq0,tpot0=tpot0,dlpot0=dlpot0,chkbnd=chkbnd,
        useext=useext,mathargs=mathargs)
    h_p = sea_h(0,0,1,salt,entr,pres,temp=temp,dliq=dliq,useext=useext)
    h_sp = sea_h(1,0,1,salt,entr,pres,temp=temp,dliq=dliq,useext=useext)
    h_ep = sea_h(0,1,1,salt,entr,pres,temp=temp,dliq=dliq,useext=useext)
    hp_s = sea_h(1,0,0,salt,entr,ppot,temp=tpot,dliq=dlpot,useext=useext)
    hp_e = sea_h(0,1,0,salt,entr,ppot,temp=tpot,dliq=dlpot,useext=useext)
    beta = (hp_s*h_ep/hp_e - h_sp) / h_p
    return beta

def contraction_theta(salt,pres,ppot,entr=None,temp=None,dliq=None,
    tpot=None,dlpot=None,chkvals=False,chktol=_CHKTOL,temp0=None,
    dliq0=None,tpot0=None,dlpot0=None,chkbnd=False,useext=False,
    mathargs=None):
    """Calculate seawater isentropic haline contraction coefficient.
    
    Calculate the haline contraction coefficient of seawater at constant
    potential temperature (entropy) from salinity, in-situ pressure,
    potential pressure, and entropy or in-situ temperature.
    
    :arg float salt: Salinity in kg/kg.
    :arg float pres: In-situ pressure in Pa.
    :arg float ppot: Potential pressure in Pa.
    :arg entr: Entropy in J/kg/K. Either this or `temp` must be
        provided.
    :type entr: float or None
    :arg temp: In-situ temperature in K. If `entr` is given, pass None
        (default) and it will be calculated.
    :type temp: float or None
    :arg dliq: In-situ liquid water density in kg/m3. If unknown, pass
        None (default) and it will be calculated.
    :type dliq: float or None
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
        (default) then `_approx_sep` is used.
    :type temp0: float or None
    :arg dliq0: Initial guess for the in-situ liquid water density in
        kg/m3. If None (default) then `_approx_sep` or
        `flu3a._dliq_default` are used.
    :type dliq0: float or None
    :arg tpot0: Initial guess for the potential temperature in K. If
        None (default) then `_approx_sep` is used.
    :type tpot0: float or None
    :arg dlpot0: Initial guess for the potential liquid water density in
        kg/m3. If None (default) then `_approx_sep` is used.
    :type dlpot0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Contraction coefficient in 1/(kg/kg).
    :raises ValueError: If both `entr` and `temp` are None.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> contraction_theta(0.035,1e7,1e5,entr=500.)
    0.717342103505
    """
    entr, temp, dliq, tpot, dlpot = eq_pot(salt,pres,ppot,entr=entr,temp=temp,
        dliq=dliq,tpot=tpot,dlpot=dlpot,chkvals=chkvals,chktol=chktol,
        temp0=temp0,dliq0=dliq0,tpot0=tpot0,dlpot0=dlpot0,chkbnd=chkbnd,
        useext=useext,mathargs=mathargs)
    h_p = sea_h(0,0,1,salt,entr,pres,temp=temp,dliq=dliq,useext=useext)
    h_sp = sea_h(1,0,1,salt,entr,pres,temp=temp,dliq=dliq,useext=useext)
    h_ep = sea_h(0,1,1,salt,entr,pres,temp=temp,dliq=dliq,useext=useext)
    hp_se = sea_h(1,1,0,salt,entr,ppot,temp=tpot,dliq=dlpot,useext=useext)
    hp_ee = sea_h(0,2,0,salt,entr,ppot,temp=tpot,dliq=dlpot,useext=useext)
    beta = (hp_se*h_ep/hp_ee - h_sp) / h_p
    return beta

def expansion_h(salt,pres,ppot,entr=None,temp=None,dliq=None,tpot=None,
    dlpot=None,chkvals=False,chktol=_CHKTOL,temp0=None,dliq0=None,
    tpot0=None,dlpot0=None,chkbnd=False,useext=False,mathargs=None):
    """Calculate seawater enthalpy expansion coefficient.
    
    Calculate the expansion coefficient of seawater with respect to
    enthalpy from salinity, in-situ pressure, potential pressure, and
    entropy or in-situ temperature.
    
    :arg float salt: Salinity in kg/kg.
    :arg float pres: In-situ pressure in Pa.
    :arg float ppot: Potential pressure in Pa.
    :arg entr: Entropy in J/kg/K. Either this or `temp` must be
        provided.
    :type entr: float or None
    :arg temp: In-situ temperature in K. If `entr` is given, pass None
        (default) and it will be calculated.
    :type temp: float or None
    :arg dliq: In-situ liquid water density in kg/m3. If unknown, pass
        None (default) and it will be calculated.
    :type dliq: float or None
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
        (default) then `_approx_sep` is used.
    :type temp0: float or None
    :arg dliq0: Initial guess for the in-situ liquid water density in
        kg/m3. If None (default) then `_approx_sep` or
        `flu3a._dliq_default` are used.
    :type dliq0: float or None
    :arg tpot0: Initial guess for the potential temperature in K. If
        None (default) then `_approx_sep` is used.
    :type tpot0: float or None
    :arg dlpot0: Initial guess for the potential liquid water density in
        kg/m3. If None (default) then `_approx_sep` is used.
    :type dlpot0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Expansion coefficient in 1/(J/kg).
    :raises ValueError: If both `entr` and `temp` are None.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> expansion_h(0.035,1e7,1e5,entr=500.)
    9.60618615640e-08
    """
    entr, temp, dliq, tpot, dlpot = eq_pot(salt,pres,ppot,entr=entr,temp=temp,
        dliq=dliq,tpot=tpot,dlpot=dlpot,chkvals=chkvals,chktol=chktol,
        temp0=temp0,dliq0=dliq0,tpot0=tpot0,dlpot0=dlpot0,chkbnd=chkbnd,
        useext=useext,mathargs=mathargs)
    h_p = sea_h(0,0,1,salt,entr,pres,temp=temp,dliq=dliq,useext=useext)
    h_ep = sea_h(0,1,1,salt,entr,pres,temp=temp,dliq=dliq,useext=useext)
    hp_e = sea_h(0,1,0,salt,entr,ppot,temp=tpot,dliq=dlpot,useext=useext)
    alpha = h_ep / (h_p * hp_e)
    return alpha

def expansion_theta(salt,pres,ppot,entr=None,temp=None,dliq=None,
    tpot=None,dlpot=None,chkvals=False,chktol=_CHKTOL,temp0=None,
    dliq0=None,tpot0=None,dlpot0=None,chkbnd=False,useext=False,
    mathargs=None):
    """Calculate seawater potential temperature expansion coefficient.
    
    Calculate the expansion coefficient of seawater with respect to
    potential temperature from salinity, in-situ pressure, potential
    pressure, and entropy or in-situ temperature.
    
    :arg float salt: Salinity in kg/kg.
    :arg float pres: In-situ pressure in Pa.
    :arg float ppot: Potential pressure in Pa.
    :arg entr: Entropy in J/kg/K. Either this or `temp` must be
        provided.
    :type entr: float or None
    :arg temp: In-situ temperature in K. If `entr` is given, pass None
        (default) and it will be calculated.
    :type temp: float or None
    :arg dliq: In-situ liquid water density in kg/m3. If unknown, pass
        None (default) and it will be calculated.
    :type dliq: float or None
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
        (default) then `_approx_sep` is used.
    :type temp0: float or None
    :arg dliq0: Initial guess for the in-situ liquid water density in
        kg/m3. If None (default) then `_approx_sep` or
        `flu3a._dliq_default` are used.
    :type dliq0: float or None
    :arg tpot0: Initial guess for the potential temperature in K. If
        None (default) then `_approx_sep` is used.
    :type tpot0: float or None
    :arg dlpot0: Initial guess for the potential liquid water density in
        kg/m3. If None (default) then `_approx_sep` is used.
    :type dlpot0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Expansion coefficient in 1/K.
    :raises ValueError: If both `entr` and `temp` are None.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> expansion_theta(0.035,1e7,1e5,entr=500.)
    3.84755380181e-04
    """
    entr, temp, dliq, tpot, dlpot = eq_pot(salt,pres,ppot,entr=entr,temp=temp,
        dliq=dliq,tpot=tpot,dlpot=dlpot,chkvals=chkvals,chktol=chktol,
        temp0=temp0,dliq0=dliq0,tpot0=tpot0,dlpot0=dlpot0,chkbnd=chkbnd,
        useext=useext,mathargs=mathargs)
    h_p = sea_h(0,0,1,salt,entr,pres,temp=temp,dliq=dliq,useext=useext)
    h_ep = sea_h(0,1,1,salt,entr,pres,temp=temp,dliq=dliq,useext=useext)
    hp_ee = sea_h(0,2,0,salt,entr,ppot,temp=tpot,dliq=dlpot,useext=useext)
    alpha = h_ep / (h_p * hp_ee)
    return alpha

