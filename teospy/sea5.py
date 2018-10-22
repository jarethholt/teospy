"""Seawater potential and conservative temperature functions.

This module provides properties of seawater related to potential
temperature and conservative temperature. Most importantly, it provides
the thermobaric and cabbeling coefficients with respect to potential and
conservative temperature, which are not available in any other module.

Everywhere in this module, potential and conservative temperature are
with respect to an absolute pressure of 1 atm.

The conservative temperature is the potential enthalpy, scaled by a
constant heat capacity to have units of temperature. The potential
enthalpy is more nearly conserved under oceanic advection than the
potential temperature, and is the recommended basis for calculating
thermobaric and cabbeling effects than potential temperature.

:Examples:

>>> tconfromtpot(0.035,300.)
300.010069445
>>> tpotfromtcon(0.035,300.010069445)
300.
>>> expansion_tcon(0.035,300.,1e8)
3.91772847589e-4
>>> expansion_tpot(0.035,300.,1e8)
3.92515634064e-4
>>> contraction_tcon(0.035,300.,1e8)
0.649596383654
>>> contraction_tpot(0.035,300.,1e8)
0.663973579411
>>> cabb_tcon(0.035,300.,1e5)
8.61252567438267e-6
>>> cabb_tpot(0.035,300.,1e5)
8.33874537690444e-6
>>> thrmb_tcon(0.035,300.,1e5)
1.48109271668362e-12
>>> thrmb_tpot(0.035,300.,1e5)
1.45941010702991e-12

:Functions:

* :func:`tconfromtpot`: Calculate conservative temperature from
  potential temperature.
* :func:`tpotfromtcon`: Calculate potential temperature from
  conservative temperature.
* :func:`expansion_tcon`: Thermal expansion coefficient with respect to
  conservative temperature.
* :func:`expansion_tpot`: Thermal expansion coefficient with respect to
  potential temperature.
* :func:`expansion_t`: Thermal expansion coefficient with respect to
  in-situ temperature.
* :func:`contraction_tcon`: Haline contraction coefficient at constant
  conservative temperature.
* :func:`contraction_tpot`: Haline contraction coefficient at constant
  potential temperature.
* :func:`contraction_t`: Haline contraction coefficient at constant
  in-situ temperature.
* :func:`cabb_tcon`: Cabbeling coefficient with respect to conservative
  temperature.
* :func:`cabb_tpot`: Cabbeling coefficient with respect to potential
  temperature.
* :func:`thrmb_tcon`: Thermobaric coefficient with respect to
  conservative temperature.
* :func:`thrmb_tpot`: Thermobaric coefficient with respect to potential
  temperature.

"""

__all__ = ['tconfromtpot','tpotfromtcon',
    'expansion_tcon','expansion_tpot','expansion_t','contraction_tcon',
    'contraction_tpot','contraction_t',
    'cabb_tcon','cabb_tpot','thrmb_tcon','thrmb_tpot']

import constants0
import maths3
import flu3a
import sea3a
import sea3c

_CHKTOL = constants0.CHKTOL
_PATM = constants0.PATM
_TCELS = constants0.TCELS
_CSEA = constants0.CSEA
_SAL0 = constants0.SAL0
_SAL1 = constants0.SAL1
_chkflubnds = constants0.chkflubnds
_chksalbnds = constants0.chksalbnds
_C_THETA = (
    (-1.446013646344788e-2,-3.305308995852924e-3,1.062415929128982e-4,
        9.477566673794488e-1,2.166591947736613e-3,3.828842955039902e-3),
    (1.,6.506097115635800e-4,3.830289486850898e-3,1.247811760368034e-6)
)
_DTEMP = 1e-3
_DSALT = 1e-5
_DPRES = 1e3
_newton = maths3.newton
_dliq_default = flu3a._dliq_default
_eq_pot = sea3c.eq_pot


## Converting potential and conservative temperature
def tconfromtpot(salt,tpot,dlpot=None,chkvals=False,chktol=_CHKTOL,
    dlpot0=None,chkbnd=False,useext=False,mathargs=None):
    """Calculate conservative from potential temperature.
    
    Calculate the conservative temperature from the potential
    temperature.
    
    :arg float salt: Salinity in kg/kg.
    :arg float tpot: Potential temperature in K.
    :arg dlpot: Seawater liquid water density at the reference pressure
        in kg/m3. If unknown, pass None (default) and it will be
        calculated.
    :type dlpot: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg dlpot0: Initial guess for the potential liquid density in
        kg/m3. If None (default) then `flu3a._dliq_default` is used.
    :type dlpot0: float or str or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg bool useext: If False (default) then the salt contribution is
        calculated from _GSCOEFFS; if True, from _GSCOEFFS_EXT.
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Conservative temperature in K.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> tconfromtpot(0.035,300.)
    300.010069445
    """
    hpot = sea3a.enthalpy(salt,tpot,_PATM,dliq=dlpot,chkvals=chkvals,
        chktol=chktol,dliq0=dlpot0,chkbnd=chkbnd,useext=useext,
        mathargs=mathargs)
    tcon = _TCELS + hpot/_CSEA
    return tcon

def _approx_stc(salt,tcon):
    """Approximate Tp from STc.
    
    Approximate the potential temperature from the salinity and
    conservative temperature.
    
    :arg float s: Salinity in kg/kg.
    :arg float ct: Conservative temperature in K.
    :returns: Potential temperature in K.
    """
    s = salt * _SAL0/_SAL1
    ct = tcon - _TCELS
    CNUMS, CDENS = _C_THETA
    th0num = (CNUMS[0] + CNUMS[1]*s + CNUMS[2]*s**2
        + CNUMS[3]*ct + CNUMS[4]*s*ct + CNUMS[5]*ct**2)
    th0den = CDENS[0] + CDENS[1]*s + CDENS[2]*ct + CDENS[3]*ct**2
    th0 = th0num / th0den
    dthdct = (CNUMS[3] + CNUMS[4]*s + 2*CNUMS[5]*ct
        - th0*(CDENS[2] + 2*CDENS[3]*ct)) / th0den
    #ct0 = tconfromtpot(salt,_TCELS+th0) - _TCELS
    #th1 = th0 + (ct-ct0)*dthdct
    th1 = th0 + (ct-th0)*dthdct
    tpot = th1 + _TCELS
    return tpot

def tpotfromtcon(salt,tcon,tpot0=None,dlpot0=None,chkbnd=False,
    useext=False,mathargs=None):
    """Calculate potential from conservative temperature.
    
    Calculate the potential temperature from the conservative
    temperature.
    
    :arg float salt: Salinity in kg/kg.
    :arg float tcon: Conservative temperature in K.
    :arg tpot0: Initial guess for the potential temperature in K. If
        None (default) then `_approx_stc` is used.
    :type tpot0: float or None
    :arg dlpot0: Initial guess for the potential liquid density in
        kg/m3. If None (default) then `flu3a._dliq_default` is used.
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
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> tpotfromtcon(0.035,300.010069445)
    300.
    """
    if tpot0 is None:
        tpot0 = _approx_stc(salt,tcon)
    hpot = _CSEA * (tcon - _TCELS)
    __, tpot, dlpot = sea3c.eq_shp(salt,_PATM,enth=hpot,temp0=tpot0,
        dliq0=dlpot0,chkbnd=chkbnd,useext=useext,mathargs=mathargs)
    return tpot


## Expansion and contraction coefficients
def expansion_tcon(salt,temp,pres,dliq=None,tpot=None,dlpot=None,
    chkvals=False,chktol=_CHKTOL,dliq0=None,tpot0=None,dlpot0=None,
    chkbnd=False,useext=False,mathargs=None):
    """Calculate thermal expansion wrt conservative temperature.
    
    Calculate the thermal expansion coefficient of seawater with respect
    to conservative temperature. This is the rate of change in the
    relative specific volume of seawater with conservative temperature,
    keeping salinity and pressure constant.
    
    :arg float salt: Salinity in kg/kg.
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg dliq: Seawater liquid water density in kg/m3. If unknown, pass
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
    :arg dliq0: Initial guess for the liquid water density in kg/m3. If
        None (default) then `flu3a._dliq_default` is used.
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
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> expansion_tcon(0.035,300.,1e8)
    3.91772847589e-4
    """
    alpha_h = sea3c.expansion_h(salt,pres,temp=temp,ppot=_PATM,dliq=dliq,
        tpot=tpot,dlpot=dlpot,chkvals=chkvals,chktol=chktol,dliq0=dliq0,
        tpot0=tpot0,dlpot0=dlpot0,chkbnd=chkbnd,useext=useext,mathargs=mathargs)
    alpha = alpha_h * _CSEA
    return alpha

def expansion_tpot(salt,temp,pres,dliq=None,tpot=None,dlpot=None,
    chkvals=False,chktol=_CHKTOL,dliq0=None,tpot0=None,dlpot0=None,
    chkbnd=False,useext=False,mathargs=None):
    """Calculate thermal expansion wrt potential temperature.
    
    Calculate the thermal expansion coefficient of seawater with respect
    to potential temperature. This is the rate of change in the relative
    specific volume of seawater with potential temperature, keeping
    salinity and pressure constant.
    
    :arg float salt: Salinity in kg/kg.
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg dliq: Seawater liquid water density in kg/m3. If unknown, pass
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
    :arg dliq0: Initial guess for the liquid water density in kg/m3. If
        None (default) then `flu3a._dliq_default` is used.
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
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> expansion_tpot(0.035,300.,1e8)
    3.92515634064e-4
    """
    alpha = sea3c.expansion_theta(salt,pres,temp=temp,ppot=_PATM,dliq=dliq,
        tpot=tpot,dlpot=dlpot,chkvals=chkvals,chktol=chktol,dliq0=dliq0,
        tpot0=tpot0,dlpot0=dlpot0,chkbnd=chkbnd,useext=useext,mathargs=mathargs)
    return alpha

def expansion_t(salt,temp,pres,dliq=None,chkvals=False,chktol=_CHKTOL,
    dliq0=None,chkbnd=False,useext=False,mathargs=None,**kwargs):
    """Calculate thermal expansion wrt in-situ temperature.
    
    Calculate the thermal expansion coefficient of seawater with respect
    to in-situ temperature. This is the simplest form of the expansion
    coefficient.
    
    :arg float salt: Salinity in kg/kg.
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg dliq: Seawater liquid water density in kg/m3. If unknown, pass
        None (default) and it will be calculated.
    :type dliq: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg dliq0: Initial guess for the liquid water density in kg/m3. If
        None (default) then `flu3a._dliq_default` is used.
    :type dliq0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg bool useext: If False (default) then the salt contribution is
        calculated from _GSCOEFFS; if True, from _GSCOEFFS_EXT.
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :arg dict kwargs: Any additional keyword arguments passed in will be
        ignored, but allow all functions in this module to be called
        with the same arguments.
    :returns: Expansion coefficient in 1/K.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> expansion_t(0.035,300.,1e8)
    3.73608885178e-4
    """
    alpha = sea3a.expansion_t(salt,temp,pres,dliq=dliq,chkvals=chkvals,
        chktol=chktol,dliq0=dliq0,chkbnd=chkbnd,useext=useext,mathargs=mathargs)
    return alpha

def contraction_tcon(salt,temp,pres,dliq=None,tpot=None,dlpot=None,
    chkvals=False,chktol=_CHKTOL,dliq0=None,tpot0=None,dlpot0=None,
    chkbnd=False,useext=False,mathargs=None):
    """Calculate haline contraction wrt conservative temperature.
    
    Calculate the haline contraction coefficient of seawater with
    respect to conservative temperature. This is the rate of change in
    the relative density of seawater with salinity, keeping the pressure
    and conservative temperature constant.
    
    :arg float salt: Salinity in kg/kg.
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg dliq: Seawater liquid water density in kg/m3. If unknown, pass
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
    :arg dliq0: Initial guess for the liquid water density in kg/m3. If
        None (default) then `flu3a._dliq_default` is used.
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
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> contraction_tcon(0.035,300.,1e8)
    0.649596383654
    """
    beta = sea3c.contraction_h(salt,pres,temp=temp,ppot=_PATM,dliq=dliq,
        tpot=tpot,dlpot=dlpot,chkvals=chkvals,chktol=chktol,dliq0=dliq0,
        tpot0=tpot0,dlpot0=dlpot0,chkbnd=chkbnd,useext=useext,mathargs=mathargs)
    return beta

def contraction_tpot(salt,temp,pres,dliq=None,tpot=None,dlpot=None,
    chkvals=False,chktol=_CHKTOL,dliq0=None,tpot0=None,dlpot0=None,
    chkbnd=False,useext=False,mathargs=None):
    """Calculate haline contraction wrt potential temperature.
    
    Calculate the haline contraction coefficient of seawater with
    respect to potential temperature. This is the rate of change in the
    relative density of seawater with salinity, keeping the pressure and
    potential temperature constant.
    
    :arg float salt: Salinity in kg/kg.
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg dliq: Seawater liquid water density in kg/m3. If unknown, pass
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
    :arg dliq0: Initial guess for the liquid water density in kg/m3. If
        None (default) then `flu3a._dliq_default` is used.
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
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> contraction_tpot(0.035,300.,1e8)
    0.663973579411
    """
    beta = sea3c.contraction_theta(salt,pres,temp=temp,ppot=_PATM,dliq=dliq,
        tpot=tpot,dlpot=dlpot,chkvals=chkvals,chktol=chktol,dliq0=dliq0,
        tpot0=tpot0,dlpot0=dlpot0,chkbnd=chkbnd,useext=useext,mathargs=mathargs)
    return beta

def contraction_t(salt,temp,pres,dliq=None,chkvals=False,chktol=_CHKTOL,
    dliq0=None,chkbnd=False,useext=False,mathargs=None,**kwargs):
    """Calculate haline contraction wrt in-situ temperature.
    
    Calculate the haline contraction coefficient of seawater with
    respect to in-situ temperature. This is the simplest form of the
    contraction coefficient.
    
    :arg float salt: Salinity in kg/kg.
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg dliq: Seawater liquid water density in kg/m3. If unknown, pass
        None (default) and it will be calculated.
    :type dliq: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg dliq0: Initial guess for the liquid water density in kg/m3. If
        None (default) then `flu3a._dliq_default` is used.
    :type dliq0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg bool useext: If False (default) then the salt contribution is
        calculated from _GSCOEFFS; if True, from _GSCOEFFS_EXT.
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :arg dict kwargs: Any additional keyword arguments passed in will be
        ignored, but allow all functions in this module to be called
        with the same arguments.
    :returns: Contraction coefficient in 1/(kg/kg).
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> contraction_t(0.035,300.,1e8)
    0.666238827368
    """
    beta = sea3a.contraction_t(salt,temp,pres,dliq=dliq,chkvals=chkvals,
        chktol=chktol,dliq0=dliq0,chkbnd=chkbnd,useext=useext,mathargs=mathargs)
    return beta


## Cabbeling and thermobaric coefficients
def cabb_tcon(salt,temp,pres,dtemp=_DTEMP,dsalt=_DSALT,dliq=None,
    tpot=None,dlpot=None,chkvals=False,chktol=_CHKTOL,dliq0=None,
    tpot0=None,dlpot0=None,chkbnd=False,useext=False,mathargs=None):
    """Calculate cabbeling coefficient wrt conservative temperature.
    
    Calculate the cabbeling coefficient of seawater with respect to
    conservative temperature. This coefficient is used in calculating
    the dianeutral velocity of seawater from gradients of conservative
    temperature.
    
    :arg float salt: Salinity in kg/kg.
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg float dtemp: Amount to change the conservative temperature by
        when calculating finite differences (default _DTEMP).
    :arg float dsalt: Amount to change the salinity by when calculating
        finite differences (default _DSALT).
    :arg dliq: Seawater liquid water density in kg/m3. If unknown, pass
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
    :arg dliq0: Initial guess for the liquid water density in kg/m3. If
        None (default) then `flu3a._dliq_default` is used.
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
    :returns: Cabbeling coefficient in 1/K^2.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    :raises RuntimeWarning: If the given salinity is less than the
        salinity difference used in finite differences (salt < dsalt).
        Salinities of (0.,dsalt) are used instead.
    
    :Examples:
    
    >>> cabb_tcon(0.035,300.,1e5)
    8.61252567438267e-6
    """
    __, __, dliq, hpot, tpot, dlpot = _eq_pot(salt,pres,_PATM,temp=temp,
        dliq=dliq,tpot=tpot,dlpot=dlpot,chkvals=chkvals,chktol=chktol,
        dliq0=dliq0,tpot0=tpot0,dlpot0=dlpot0,chkbnd=chkbnd,useext=useext,
        mathargs=mathargs)
    tcon = hpot/_CSEA + _TCELS
    alpha0 = expansion_tcon(salt,temp,pres,dliq=dliq,tpot=tpot,dlpot=dlpot,
        useext=useext)
    beta0 = contraction_tcon(salt,temp,pres,dliq=dliq,tpot=tpot,dlpot=dlpot,
        useext=useext)
    ratio = alpha0/beta0
    
    # Calculate finite differences wrt temperature
    tcon_l = tcon - dtcon
    tcon_u = tcon + dtcon
    hpot_l = _CSEA*(tcon_l - _TCELS)
    hpot_u = _CSEA*(tcon_u - _TCELS)
    __, t_l, d_l, __, tp_l, dp_l = _eq_pot(salt,pres,_PATM,hpot=hpot_l,
        useext=useext,mathargs=mathargs)
    __, t_u, d_u, __, tp_u, dp_u = _eq_pot(salt,pres,_PATM,hpot=hpot_u,
        useext=useext,mathargs=mathargs)
    alpha_l = expansion_tcon(salt,t_l,pres,dliq=d_l,tpot=tp_l,dlpot=dp_l,
        useext=useext)
    alpha_u = expansion_tcon(salt,t_u,pres,dliq=d_u,tpot=tp_u,dlpot=dp_u,
        useext=useext)
    dalpha_tc = (alpha_u - alpha_l) / (tcon_u - tcon_l)
    
    # Calculate finite differences wrt salinity
    if salt < dsalt:
        warnmsg = ('Salinity {0} is less than change in salinity {1}. Using '
            'salinities of (0.,{1}) for finite difference').format(salt,dsalt)
        warnings.warn(warnmsg,RuntimeWarning)
        salt_l = 0.
        salt_u = dsalt
    else:
        salt_l = salt - dsalt
        salt_u = salt + dsalt
    __, t_l, d_l, __, tp_l, dp_l = _eq_pot(salt_l,pres,_PATM,hpot=hpot,
        useext=useext,mathargs=mathargs)
    __, t_u, d_u, __, tp_u, dp_u = _eq_pot(salt_u,pres,_PATM,hpot=hpot,
        useext=useext,mathargs=mathargs)
    alpha_l = expansion_tcon(salt_l,t_l,pres,dliq=d_l,tpot=tp_l,dlpot=dp_l,
        useext=useext)
    alpha_u = expansion_tcon(salt_u,t_u,pres,dliq=d_u,tpot=tp_u,dlpot=dp_u,
        useext=useext)
    dalpha_s = (alpha_u - alpha_l) / (salt_u - salt_l)
    beta_l = contraction_tcon(salt_l,t_l,pres,dliq=d_l,tpot=tp_l,dlpot=dp_l,
        useext=useext)
    beta_u = contraction_tcon(salt_u,t_u,pres,dliq=d_u,tpot=tp_u,dlpot=dp_u,
        useext=useext)
    dbeta_s = (beta_u - beta_l) / (salt_u - salt_l)
    cabb = dalpha_tc + 2*ratio*dalpha_s - ratio**2*dbeta_s
    return cabb

def cabb_tpot(salt,temp,pres,dtemp=_DTEMP,dsalt=_DSALT,dliq=None,
    tpot=None,dlpot=None,chkvals=False,chktol=_CHKTOL,dliq0=None,
    tpot0=None,dlpot0=None,chkbnd=False,useext=False,mathargs=None):
    """Calculate cabbeling coefficient wrt potential temperature.
    
    Calculate the cabbeling coefficient of seawater with respect to
    potential temperature. This coefficient is used in calculating the
    dianeutral velocity of seawater from gradients of potential
    temperature.
    
    :arg float salt: Salinity in kg/kg.
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg float dtemp: Amount to change the conservative temperature by
        when calculating finite differences (default _DTEMP).
    :arg float dsalt: Amount to change the salinity by when calculating
        finite differences (default _DSALT).
    :arg dliq: Seawater liquid water density in kg/m3. If unknown, pass
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
    :arg dliq0: Initial guess for the liquid water density in kg/m3. If
        None (default) then `flu3a._dliq_default` is used.
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
    :returns: Cabbeling coefficient in 1/K^2.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    :raises RuntimeWarning: If the given salinity is less than the
        salinity difference used in finite differences (salt < dsalt).
        Salinities of (0.,dsalt) are used instead.
    
    :Examples:
    
    >>> cabb_tpot(0.035,300.,1e5)
    8.33874537690444e-6
    """
    __, __, dliq, __, tpot, dlpot = _eq_pot(salt,pres,_PATM,temp=temp,
        dliq=dliq,tpot=tpot,dlpot=dlpot,chkvals=chkvals,chktol=chktol,
        dliq0=dliq0,tpot0=tpot0,dlpot0=dlpot0,chkbnd=chkbnd,useext=useext,
        mathargs=mathargs)
    alpha0 = expansion_tpot(salt,temp,pres,dliq=dliq,tpot=tpot,dlpot=dlpot,
        useext=useext)
    beta0 = contraction_tpot(salt,temp,pres,dliq=dliq,tpot=tpot,dlpot=dlpot,
        useext=useext)
    ratio = alpha0/beta0
    
    # Calculate finite differences wrt temperature
    tpot_l = tpot - dtpot
    tpot_u = tpot + dtpot
    __, t_l, d_l, __, __, dp_l = _eq_pot(salt,pres,_PATM,tpot=tpot_l,
        useext=useext,mathargs=mathargs)
    __, t_u, d_u, __, __, dp_u = _eq_pot(salt,pres,_PATM,tpot=tpot_u,
        useext=useext,mathargs=mathargs)
    alpha_l = expansion_tpot(salt,t_l,pres,dliq=d_l,tpot=tpot_l,dlpot=dp_l,
        useext=useext)
    alpha_u = expansion_tpot(salt,t_u,pres,dliq=d_u,tpot=tpot_u,dlpot=dp_u,
        useext=useext)
    dalpha_tp = (alpha_u - alpha_l) / (tpot_u - tpot_l)
    
    # Calculate finite differences wrt salinity
    if salt < dsalt:
        warnmsg = ('Salinity {0} is less than change in salinity {1}. Using '
            'salinities of (0.,{1}) for finite difference').format(salt,dsalt)
        warnings.warn(warnmsg,RuntimeWarning)
        salt_l = 0.
        salt_u = dsalt
    else:
        salt_l = salt - dsalt
        salt_u = salt + dsalt
    __, t_l, d_l, __, __, dp_l = _eq_pot(salt_l,pres,_PATM,tpot=tpot,
        useext=useext,mathargs=mathargs)
    __, t_u, d_u, __, __, dp_u = _eq_pot(salt_u,pres,_PATM,tpot=tpot,
        useext=useext,mathargs=mathargs)
    alpha_l = expansion_tpot(salt_l,t_l,pres,dliq=d_l,tpot=tpot,dlpot=dp_l,
        useext=useext)
    alpha_u = expansion_tpot(salt_u,t_u,pres,dliq=d_u,tpot=tpot,dlpot=dp_u,
        useext=useext)
    dalpha_s = (alpha_u - alpha_l) / (salt_u - salt_l)
    beta_l = contraction_tpot(salt_l,t_l,pres,dliq=d_l,tpot=tpot,dlpot=dp_l,
        useext=useext)
    beta_u = contraction_tpot(salt_u,t_u,pres,dliq=d_u,tpot=tpot,dlpot=dp_u,
        useext=useext)
    dbeta_s = (beta_u - beta_l) / (salt_u - salt_l)
    cabb = dalpha_tp + 2*ratio*dalpha_s - ratio**2*dbeta_s
    return cabb

def thrmb_tcon(salt,temp,pres,dpres=_DPRES,dliq=None,tpot=None,
    dlpot=None,chkvals=False,chktol=_CHKTOL,dliq0=None,tpot0=None,
    dlpot0=None,chkbnd=False,useext=False,mathargs=None):
    """Calculate thermobaric coefficient wrt conservative temperature.
    
    Calculate the thermobaric coefficient of seawater with respect to
    conservative temperature. This coefficient is used in calculating
    the dianeutral velocity of seawater from gradients of pressure and
    conservative temperature.
    
    :arg float salt: Salinity in kg/kg.
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg float dpres: Amount to change the pressure by when calculating
        finite differences (default _DPRES).
    :arg dliq: Seawater liquid water density in kg/m3. If unknown, pass
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
    :arg dliq0: Initial guess for the liquid water density in kg/m3. If
        None (default) then `flu3a._dliq_default` is used.
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
    :returns: Thermobaric coefficient in 1/K/Pa.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    :raises RuntimeWarning: If the given pressure is less than the
        pressure difference used in finite differences (pres <= dpres).
        Pressures of (.75*pres,dpres) are used instead.
    
    :Examples:
    
    >>> thrmb_tcon(0.035,300.,1e5)
    1.48109271668362e-12
    """
    __, __, dliq, hpot, tpot, dlpot = _eq_pot(salt,pres,_PATM,temp=temp,
        dliq=dliq,tpot=tpot,dlpot=dlpot,chkvals=chkvals,chktol=chktol,
        dliq0=dliq0,tpot0=tpot0,dlpot0=dlpot0,chkbnd=chkbnd,useext=useext,
        mathargs=mathargs)
    beta0 = contraction_tcon(salt,temp,pres,dliq=dliq,tpot=tpot,dlpot=dlpot,
        useext=useext)
    
    # Calculate finite differences wrt pressure
    if pres <= dpres:
        pres_l = .75*pres
        pres_u = dpres
        warnmsg = ('Pressure {0} is less than change in pressure {1}. '
            'Using pressures of ({2},{1}) for finite '
            'difference').format(pres,dpres,pres_l)
        warnings.warn(warnmsg,RuntimeWarning)
    else:
        pres_l = pres - dpres
        pres_u = pres + dpres
    __, t_l, d_l, __, __, __ = _eq_pot(salt,pres_l,_PATM,hpot=hpot,
        useext=useext,mathargs=mathargs)
    __, t_u, d_u, __, __, __ = _eq_pot(salt,pres_u,_PATM,hpot=hpot,
        useext=useext,mathargs=mathargs)
    alpha_l = expansion_tcon(salt,t_l,pres_l,dliq=d_l,tpot=tpot,dlpot=dlpot,
        useext=useext)
    alpha_u = expansion_tcon(salt,t_u,pres_u,dliq=d_u,tpot=tpot,dlpot=dlpot,
        useext=useext)
    beta_l = contraction_tcon(salt,t_l,pres_l,dliq=d_l,tpot=tpot,dlpot=dlpot,
        useext=useext)
    beta_u = contraction_tcon(salt,t_u,pres_u,dliq=d_u,tpot=tpot,dlpot=dlpot,
        useext=useext)
    thermb = beta0 * (alpha_u/beta_u - alpha_l/beta_l) / (pres_u - pres_l)
    return thermb

def thrmb_tpot(salt,temp,pres,dpres=_DPRES,dliq=None,tpot=None,
    dlpot=None,chkvals=False,chktol=_CHKTOL,dliq0=None,tpot0=None,
    dlpot0=None,chkbnd=False,useext=False,mathargs=None):
    """Calculate thermobaric coefficient wrt potential temperature.
    
    Calculate the thermobaric coefficient of seawater with respect to
    potential temperature. This coefficient is used in calculating the
    dianeutral velocity of seawater from gradients of pressure and
    potential temperature.
    
    :arg float salt: Salinity in kg/kg.
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg float dpres: Amount to change the pressure by when calculating
        finite differences (default _DPRES).
    :arg dliq: Seawater liquid water density in kg/m3. If unknown, pass
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
    :arg dliq0: Initial guess for the liquid water density in kg/m3. If
        None (default) then `flu3a._dliq_default` is used.
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
    :returns: Thermobaric coefficient in 1/K/Pa.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    :raises RuntimeWarning: If the given pressure is less than the
        pressure difference used in finite differences (pres <= dpres).
        Pressures of (.75*pres,dpres) are used instead.
    
    :Examples:
    
    >>> thrmb_tpot(0.035,300.,1e5)
    1.45941010702991e-12
    """
    __, __, dliq, __, tpot, dlpot = _eq_pot(salt,pres,_PATM,temp=temp,
        dliq=dliq,tpot=tpot,dlpot=dlpot,chkvals=chkvals,chktol=chktol,
        dliq0=dliq0,tpot0=tpot0,dlpot0=dlpot0,chkbnd=chkbnd,useext=useext,
        mathargs=mathargs)
    beta0 = contraction_tpot(salt,temp,pres,dliq=dliq,tpot=tpot,dlpot=dlpot,
        useext=useext)
    
    # Calculate finite differences wrt pressure
    if pres <= dpres:
        pres_l = .75*pres
        pres_u = dpres
        warnmsg = ('Pressure {0} is less than change in pressure {1}. '
            'Using pressures of ({2},{1}) for finite '
            'difference').format(pres,dpres,pres_l)
        warnings.warn(warnmsg,RuntimeWarning)
    else:
        pres_l = pres - dpres
        pres_u = pres + dpres
    __, t_l, d_l, __, __, __ = _eq_pot(salt,pres_l,_PATM,tpot=tpot,
        useext=useext,mathargs=mathargs)
    __, t_u, d_u, __, __, __ = _eq_pot(salt,pres_u,_PATM,tpot=tpot,
        useext=useext,mathargs=mathargs)
    alpha_l = expansion_tpot(salt,t_l,pres_l,dliq=d_l,tpot=tpot,dlpot=dlpot,
        useext=useext)
    alpha_u = expansion_tpot(salt,t_u,pres_u,dliq=d_u,tpot=tpot,dlpot=dlpot,
        useext=useext)
    beta_l = contraction_tpot(salt,t_l,pres_l,dliq=d_l,tpot=tpot,dlpot=dlpot,
        useext=useext)
    beta_u = contraction_tpot(salt,t_u,pres_u,dliq=d_u,tpot=tpot,dlpot=dlpot,
        useext=useext)
    thermb = beta0 * (alpha_u/beta_u - alpha_l/beta_l) / (pres_u - pres_l)
    return thermb

