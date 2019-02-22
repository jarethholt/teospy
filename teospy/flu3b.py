"""Fluid water properties from the Gibbs energy.

This module provides thermodynamic properties (e.g. heat capacity) of both liquid water and water vapour based on the Gibbs free energy, with temperature and pressure as the primary variables.

:Examples:

>>> liq_cp(300.,1e5)
4180.63952202
>>> liq_expansion(300.,1e5)
2.74803716256e-04
>>> liq_kappa_t(300.,1e5)
4.50515304336e-10
>>> vap_cp(300.,1e3)
1874.12149028
>>> vap_expansion(300.,1e3)
3.34352010567e-03
>>> vap_kappa_t(300.,1e3)
1.00048646242e-03

:Functions:

* :func:`liq_cp`: Liquid water isobaric heat capacity.
* :func:`liq_cv`: Liquid water isochoric heat capacity.
* :func:`liq_density`: Liquid water density.
* :func:`liq_enthalpy`: Liquid water enthalpy.
* :func:`liq_entropy`: Liquid water entropy.
* :func:`liq_expansion`: Liquid water thermal expansion coefficient.
* :func:`liq_gibbsenergy`: Liquid water Gibbs free energy.
* :func:`liq_helmholtzenergy`: Liquid water Helmholtz free energy.
* :func:`liq_internalenergy`: Liquid water internal energy.
* :func:`liq_kappa_s`: Liquid water isentropic compressibility.
* :func:`liq_kappa_t`: Liquid water isothermal compressibility.
* :func:`liq_lapserate`: Liquid water adiabatic lapse rate.
* :func:`liq_soundspeed`: Liquid water sound speed.
* :func:`vap_cp`: Water vapour isobaric heat capacity.
* :func:`vap_cv`: Water vapour isochoric heat capacity.
* :func:`vap_density`: Water vapour density.
* :func:`vap_enthalpy`: Water vapour enthalpy.
* :func:`vap_entropy`: Water vapour entropy.
* :func:`vap_expansion`: Water vapour thermal expansion coefficient.
* :func:`vap_gibbsenergy`: Water vapour Gibbs free energy.
* :func:`vap_helmholtzenergy`: Water vapour Helmholtz free energy.
* :func:`vap_internalenergy`: Water vapour internal energy.
* :func:`vap_kappa_s`: Water vapour isentropic compressibility.
* :func:`vap_kappa_t`: Water vapour isothermal compressibility.
* :func:`vap_lapserate`: Water vapour adiabatic lapse rate.
* :func:`vap_soundspeed`: Water vapour sound speed.

"""

__all__ = ['liq_cp','liq_cv','liq_density','liq_enthalpy','liq_entropy',
    'liq_expansion','liq_gibbsenergy','liq_helmholtzenergy',
    'liq_internalenergy','liq_kappa_s','liq_kappa_t','liq_lapserate',
    'liq_soundspeed','vap_cp','vap_cv','vap_density','vap_enthalpy',
    'vap_entropy','vap_expansion','vap_gibbsenergy','vap_helmholtzenergy',
    'vap_internalenergy','vap_kappa_s','vap_kappa_t','vap_lapserate',
    'vap_soundspeed']

from teospy import constants0
from teospy import flu2
from teospy import flu3a

_CHKTOL = constants0.CHKTOL
_eq_tp_liq = flu3a.eq_tp_liq
_eq_tp_vap = flu3a.eq_tp_vap
_liq_g = flu3a.liq_g
_vap_g = flu3a.vap_g


## Liquid water functions
def liq_cp(temp,pres,dliq=None,chkvals=False,chktol=_CHKTOL,dliq0=None,
    chkbnd=False,mathargs=None):
    """Calculate liquid water isobaric heat capacity.
    
    Calculate the isobaric (constant pressureheat capacity of liquid
    water at temperature and pressure.
    
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg dliq: Liquid water density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dliq: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg dliq0: Initial guess for the liquid water density in kg/m3. If
        None (default) then `_dliq_default` is used. A string specifier
        of the method (e.g. 'crit') can also be passed. See _LIQMETHODS
        in flu3a for valid specifiers.
    :type dliq0: float or str or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Heat capacity in J/kg/K.
    :raises RuntimeWarning: If a string is passed for `dliq0` that does
        not match an available method. The default is used instead.
    :raises UserWarning: If a string is passed for `dliq0` that
        specifies a function intended for water vapour.
    :raises RuntimeWarning: If the value of dliq is more consistent with
        water vapour in the subcritical region.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> liq_cp(300.,1e5)
    4180.63952202
    """
    dliq = _eq_tp_liq(temp,pres,dliq=dliq,chkvals=chkvals,chktol=chktol,
        dliq0=dliq0,chkbnd=chkbnd,mathargs=mathargs)
    g_tt = _liq_g(2,0,temp,pres,dliq=dliq)
    cp = -temp*g_tt
    return cp

def liq_cv(temp,pres,dliq=None,chkvals=False,chktol=_CHKTOL,dliq0=None,
    chkbnd=False,mathargs=None):
    """Calculate liquid water isochoric heat capacity.
    
    Calculate the isochoric (constant volume) heat capacity of liquid
    water at temperature and pressure.
    
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg dliq: Liquid water density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dliq: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg dliq0: Initial guess for the liquid water density in kg/m3. If
        None (default) then `_dliq_default` is used. A string specifier
        of the method (e.g. 'crit') can also be passed. See _LIQMETHODS
        in flu3a for valid specifiers.
    :type dliq0: float or str or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Heat capacity in J/kg/K.
    :raises RuntimeWarning: If a string is passed for `dliq0` that does
        not match an available method. The default is used instead.
    :raises UserWarning: If a string is passed for `dliq0` that
        specifies a function intended for water vapour.
    :raises RuntimeWarning: If the value of dliq is more consistent with
        water vapour in the subcritical region.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> liq_cv(300.,1e5)
    4130.17861503
    """
    dliq = _eq_tp_liq(temp,pres,dliq=dliq,chkvals=chkvals,chktol=chktol,
        dliq0=dliq0,chkbnd=chkbnd,mathargs=mathargs)
    g_tt = _liq_g(2,0,temp,pres,dliq=dliq)
    g_tp = _liq_g(1,1,temp,pres,dliq=dliq)
    g_pp = _liq_g(0,2,temp,pres,dliq=dliq)
    cv = temp * (-g_tt + g_tp**2/g_pp)
    return cv

def liq_density(temp,pres,dliq=None,chkvals=False,chktol=_CHKTOL,
    dliq0=None,chkbnd=False,mathargs=None):
    """Calculate liquid water density.
    
    Calculate the density of liquid water at temperature and pressure.
    
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg dliq: Liquid water density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dliq: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg dliq0: Initial guess for the liquid water density in kg/m3. If
        None (default) then `_dliq_default` is used. A string specifier
        of the method (e.g. 'crit') can also be passed. See _LIQMETHODS
        in flu3a for valid specifiers.
    :type dliq0: float or str or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Density in kg/m3.
    :raises RuntimeWarning: If a string is passed for `dliq0` that does
        not match an available method. The default is used instead.
    :raises UserWarning: If a string is passed for `dliq0` that
        specifies a function intended for water vapour.
    :raises RuntimeWarning: If the value of dliq is more consistent with
        water vapour in the subcritical region.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    """
    dliq = _eq_tp_liq(temp,pres,dliq=dliq,chkvals=chkvals,chktol=chktol,
        dliq0=dliq0,chkbnd=chkbnd,mathargs=mathargs)
    return dliq

def liq_enthalpy(temp,pres,dliq=None,chkvals=False,chktol=_CHKTOL,
    dliq0=None,chkbnd=False,mathargs=None):
    """Calculate liquid water enthalpy.
    
    Calculate the specific enthalpy of liquid water at temperature and
    pressure.
    
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg dliq: Liquid water density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dliq: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg dliq0: Initial guess for the liquid water density in kg/m3. If
        None (default) then `_dliq_default` is used. A string specifier
        of the method (e.g. 'crit') can also be passed. See _LIQMETHODS
        in flu3a for valid specifiers.
    :type dliq0: float or str or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Enthalpy in J/kg.
    :raises RuntimeWarning: If a string is passed for `dliq0` that does
        not match an available method. The default is used instead.
    :raises UserWarning: If a string is passed for `dliq0` that
        specifies a function intended for water vapour.
    :raises RuntimeWarning: If the value of dliq is more consistent with
        water vapour in the subcritical region.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> liq_enthalpy(300.,1e5)
    112653.679689
    """
    dliq = _eq_tp_liq(temp,pres,dliq=dliq,chkvals=chkvals,chktol=chktol,
        dliq0=dliq0,chkbnd=chkbnd,mathargs=mathargs)
    g = _liq_g(0,0,temp,pres,dliq=dliq)
    g_t = _liq_g(1,0,temp,pres,dliq=dliq)
    h = g - temp*g_t
    return h

def liq_entropy(temp,pres,dliq=None,chkvals=False,chktol=_CHKTOL,
    dliq0=None,chkbnd=False,mathargs=None):
    """Calculate liquid water entropy.
    
    Calculate the specific entropy of liquid water at temperature and
    pressure.
    
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg dliq: Liquid water density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dliq: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg dliq0: Initial guess for the liquid water density in kg/m3. If
        None (default) then `_dliq_default` is used. A string specifier
        of the method (e.g. 'crit') can also be passed. See _LIQMETHODS
        in flu3a for valid specifiers.
    :type dliq0: float or str or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Entropy in J/kg/K.
    :raises RuntimeWarning: If a string is passed for `dliq0` that does
        not match an available method. The default is used instead.
    :raises UserWarning: If a string is passed for `dliq0` that
        specifies a function intended for water vapour.
    :raises RuntimeWarning: If the value of dliq is more consistent with
        water vapour in the subcritical region.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> liq_entropy(300.,1e5)
    393.062433815
    """
    dliq = _eq_tp_liq(temp,pres,dliq=dliq,chkvals=chkvals,chktol=chktol,
        dliq0=dliq0,chkbnd=chkbnd,mathargs=mathargs)
    g_t = _liq_g(1,0,temp,pres,dliq=dliq)
    s = -g_t
    return s

def liq_expansion(temp,pres,dliq=None,chkvals=False,chktol=_CHKTOL,
    dliq0=None,chkbnd=False,mathargs=None):
    """Calculate liquid water thermal expansion coefficient.
    
    Calculate the thermal expansion coefficient of liquid water at
    temperature and pressure.
    
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg dliq: Liquid water density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dliq: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg dliq0: Initial guess for the liquid water density in kg/m3. If
        None (default) then `_dliq_default` is used. A string specifier
        of the method (e.g. 'crit') can also be passed. See _LIQMETHODS
        in flu3a for valid specifiers.
    :type dliq0: float or str or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Expansion coefficient in 1/K.
    :raises RuntimeWarning: If a string is passed for `dliq0` that does
        not match an available method. The default is used instead.
    :raises UserWarning: If a string is passed for `dliq0` that
        specifies a function intended for water vapour.
    :raises RuntimeWarning: If the value of dliq is more consistent with
        water vapour in the subcritical region.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> liq_expansion(300.,1e5)
    2.74803716256e-04
    """
    dliq = _eq_tp_liq(temp,pres,dliq=dliq,chkvals=chkvals,chktol=chktol,
        dliq0=dliq0,chkbnd=chkbnd,mathargs=mathargs)
    g_p = _liq_g(0,1,temp,pres,dliq=dliq)
    g_tp = _liq_g(1,1,temp,pres,dliq=dliq)
    alpha = g_tp / g_p
    return alpha

def liq_gibbsenergy(temp,pres,dliq=None,chkvals=False,chktol=_CHKTOL,
    dliq0=None,chkbnd=False,mathargs=None):
    """Calculate liquid water Gibbs free energy.
    
    Calculate the specific Gibbs free energy of liquid water at
    temperature and pressure.
    
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg dliq: Liquid water density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dliq: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg dliq0: Initial guess for the liquid water density in kg/m3. If
        None (default) then `_dliq_default` is used. A string specifier
        of the method (e.g. 'crit') can also be passed. See _LIQMETHODS
        in flu3a for valid specifiers.
    :type dliq0: float or str or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Gibbs energy in J/kg.
    :raises RuntimeWarning: If a string is passed for `dliq0` that does
        not match an available method. The default is used instead.
    :raises UserWarning: If a string is passed for `dliq0` that
        specifies a function intended for water vapour.
    :raises RuntimeWarning: If the value of dliq is more consistent with
        water vapour in the subcritical region.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> liq_gibbsenergy(300.,1e5)
    -5265.05045577
    """
    dliq = _eq_tp_liq(temp,pres,dliq=dliq,chkvals=chkvals,chktol=chktol,
        dliq0=dliq0,chkbnd=chkbnd,mathargs=mathargs)
    g = _liq_g(0,0,temp,pres,dliq=dliq)
    return g

def liq_helmholtzenergy(temp,pres,dliq=None,chkvals=False,
    chktol=_CHKTOL,dliq0=None,chkbnd=False,mathargs=None):
    """Calculate liquid water Helmholtz free energy.
    
    Calculate the specific Helmholtz free energy of liquid water at
    temperature and pressure.
    
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg dliq: Liquid water density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dliq: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg dliq0: Initial guess for the liquid water density in kg/m3. If
        None (default) then `_dliq_default` is used. A string specifier
        of the method (e.g. 'crit') can also be passed. See _LIQMETHODS
        in flu3a for valid specifiers.
    :type dliq0: float or str or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Helmholtz energy in J/kg.
    :raises RuntimeWarning: If a string is passed for `dliq0` that does
        not match an available method. The default is used instead.
    :raises UserWarning: If a string is passed for `dliq0` that
        specifies a function intended for water vapour.
    :raises RuntimeWarning: If the value of dliq is more consistent with
        water vapour in the subcritical region.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    """
    dliq = _eq_tp_liq(temp,pres,dliq=dliq,chkvals=chkvals,chktol=chktol,
        dliq0=dliq0,chkbnd=chkbnd,mathargs=mathargs)
    g = _liq_g(0,0,temp,pres,dliq=dliq)
    g_p = _liq_g(0,1,temp,pres,dliq=dliq)
    f = g - pres*g_p
    return f

def liq_internalenergy(temp,pres,dliq=None,chkvals=False,chktol=_CHKTOL,
    dliq0=None,chkbnd=False,mathargs=None):
    """Calculate liquid water internal energy.
    
    Calculate the specific internal energy of liquid water at
    temperature and pressure.
    
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg dliq: Liquid water density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dliq: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg dliq0: Initial guess for the liquid water density in kg/m3. If
        None (default) then `_dliq_default` is used. A string specifier
        of the method (e.g. 'crit') can also be passed. See _LIQMETHODS
        in flu3a for valid specifiers.
    :type dliq0: float or str or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Internal energy in J/kg.
    :raises RuntimeWarning: If a string is passed for `dliq0` that does
        not match an available method. The default is used instead.
    :raises UserWarning: If a string is passed for `dliq0` that
        specifies a function intended for water vapour.
    :raises RuntimeWarning: If the value of dliq is more consistent with
        water vapour in the subcritical region.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> liq_internalenergy(300.,1e5)
    112553.334133
    """
    dliq = _eq_tp_liq(temp,pres,dliq=dliq,chkvals=chkvals,chktol=chktol,
        dliq0=dliq0,chkbnd=chkbnd,mathargs=mathargs)
    g = _liq_g(0,0,temp,pres,dliq=dliq)
    g_t = _liq_g(1,0,temp,pres,dliq=dliq)
    g_p = _liq_g(0,1,temp,pres,dliq=dliq)
    u = g - temp*g_t - pres*g_p
    return u

def liq_kappa_s(temp,pres,dliq=None,chkvals=False,chktol=_CHKTOL,
    dliq0=None,chkbnd=False,mathargs=None):
    """Calculate liquid water isentropic compressibility.
    
    Calculate the isentropic compressibility of liquid water at
    temperature and pressure.
    
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg dliq: Liquid water density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dliq: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg dliq0: Initial guess for the liquid water density in kg/m3. If
        None (default) then `_dliq_default` is used. A string specifier
        of the method (e.g. 'crit') can also be passed. See _LIQMETHODS
        in flu3a for valid specifiers.
    :type dliq0: float or str or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Compressibility in 1/Pa.
    :raises RuntimeWarning: If a string is passed for `dliq0` that does
        not match an available method. The default is used instead.
    :raises UserWarning: If a string is passed for `dliq0` that
        specifies a function intended for water vapour.
    :raises RuntimeWarning: If the value of dliq is more consistent with
        water vapour in the subcritical region.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> liq_kappa_s(300.,1e5)
    4.45077521253e-10
    """
    dliq = _eq_tp_liq(temp,pres,dliq=dliq,chkvals=chkvals,chktol=chktol,
        dliq0=dliq0,chkbnd=chkbnd,mathargs=mathargs)
    g_p = _liq_g(0,1,temp,pres,dliq=dliq)
    g_tt = _liq_g(2,0,temp,pres,dliq=dliq)
    g_tp = _liq_g(1,1,temp,pres,dliq=dliq)
    g_pp = _liq_g(0,2,temp,pres,dliq=dliq)
    kappa = dliq*(-g_pp - g_tp**2/(-g_tt))
    return kappa

def liq_kappa_t(temp,pres,dliq=None,chkvals=False,chktol=_CHKTOL,
    dliq0=None,chkbnd=False,mathargs=None):
    """Calculate liquid water isothermal compressibility.
    
    Calculate the isothermal compressibility of liquid water at
    temperature and pressure.
    
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg dliq: Liquid water density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dliq: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg dliq0: Initial guess for the liquid water density in kg/m3. If
        None (default) then `_dliq_default` is used. A string specifier
        of the method (e.g. 'crit') can also be passed. See _LIQMETHODS
        in flu3a for valid specifiers.
    :type dliq0: float or str or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Compressibility in 1/Pa.
    :raises RuntimeWarning: If a string is passed for `dliq0` that does
        not match an available method. The default is used instead.
    :raises UserWarning: If a string is passed for `dliq0` that
        specifies a function intended for water vapour.
    :raises RuntimeWarning: If the value of dliq is more consistent with
        water vapour in the subcritical region.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> liq_kappa_t(300.,1e5)
    4.50515304336e-10
    """
    dliq = _eq_tp_liq(temp,pres,dliq=dliq,chkvals=chkvals,chktol=chktol,
        dliq0=dliq0,chkbnd=chkbnd,mathargs=mathargs)
    g_p = _liq_g(0,1,temp,pres,dliq=dliq)
    g_pp = _liq_g(0,2,temp,pres,dliq=dliq)
    kappa = -g_pp / g_p
    return kappa

def liq_lapserate(temp,pres,dliq=None,chkvals=False,chktol=_CHKTOL,
    dliq0=None,chkbnd=False,mathargs=None):
    """Calculate liquid water adiabatic lapse rate.
    
    Calculate the adiabatic lapse rate of liquid water at temperature
    and pressure.
    
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg dliq: Liquid water density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dliq: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg dliq0: Initial guess for the liquid water density in kg/m3. If
        None (default) then `_dliq_default` is used. A string specifier
        of the method (e.g. 'crit') can also be passed. See _LIQMETHODS
        in flu3a for valid specifiers.
    :type dliq0: float or str or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Lapse rate in K/Pa.
    :raises RuntimeWarning: If a string is passed for `dliq0` that does
        not match an available method. The default is used instead.
    :raises UserWarning: If a string is passed for `dliq0` that
        specifies a function intended for water vapour.
    :raises RuntimeWarning: If the value of dliq is more consistent with
        water vapour in the subcritical region.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> liq_lapserate(300.,1e5)
    1.97878804448e-08
    """
    dliq = _eq_tp_liq(temp,pres,dliq=dliq,chkvals=chkvals,chktol=chktol,
        dliq0=dliq0,chkbnd=chkbnd,mathargs=mathargs)
    g_tt = _liq_g(2,0,temp,pres,dliq=dliq)
    g_tp = _liq_g(1,1,temp,pres,dliq=dliq)
    gamma = g_tp/(-g_tt)
    return gamma

def liq_soundspeed(temp,pres,dliq=None,chkvals=False,chktol=_CHKTOL,
    dliq0=None,chkbnd=False,mathargs=None):
    """Calculate liquid water sound speed.
    
    Calculate the speed of sound in liquid water at temperature and
    pressure.
    
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg dliq: Liquid water density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dliq: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg dliq0: Initial guess for the liquid water density in kg/m3. If
        None (default) then `_dliq_default` is used. A string specifier
        of the method (e.g. 'crit') can also be passed. See _LIQMETHODS
        in flu3a for valid specifiers.
    :type dliq0: float or str or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Sound speed in m/s.
    :raises RuntimeWarning: If a string is passed for `dliq0` that does
        not match an available method. The default is used instead.
    :raises UserWarning: If a string is passed for `dliq0` that
        specifies a function intended for water vapour.
    :raises RuntimeWarning: If the value of dliq is more consistent with
        water vapour in the subcritical region.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> liq_soundspeed(300.,1e5)
    1501.52041506
    """
    dliq = _eq_tp_liq(temp,pres,dliq=dliq,chkvals=chkvals,chktol=chktol,
        dliq0=dliq0,chkbnd=chkbnd,mathargs=mathargs)
    g_p = _liq_g(0,1,temp,pres,dliq=dliq)
    g_tt = _liq_g(2,0,temp,pres,dliq=dliq)
    g_tp = _liq_g(1,1,temp,pres,dliq=dliq)
    g_pp = _liq_g(0,2,temp,pres,dliq=dliq)
    dvdp = g_pp + g_tp**2/(-g_tt)
    csqinv = -dliq**2 * dvdp
    c = csqinv**(-.5)
    return c


## Water vapour functions
def vap_cp(temp,pres,dvap=None,chkvals=False,chktol=_CHKTOL,dvap0=None,
    chkbnd=False,mathargs=None):
    """Calculate water vapour isobaric heat capacity.
    
    Calculate the isobaric heat capacity of water vapour at temperature
    and pressure.
    
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg dvap: Water vapour density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dvap: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg dvap0: Initial guess for the water vapour density in kg/m3. If
        None (default) then `_dvap_default` is used. A string specifier
        of the method (e.g. 'crit') can also be passed. See _VAPMETHODS
        in flu3a for valid specifiers.
    :type dvap0: float or str or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Heat capacity in J/kg/K.
    :raises RuntimeWarning: If a string is passed for `dvap0` that does
        not match an available method. The default is used instead.
    :raises UserWarning: If a string is passed for `dvap0` that
        specifies a function intended for water vapour.
    :raises RuntimeWarning: If the value of dvap is more consistent with
        water vapour in the subcritical region.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> vap_cp(300.,1e3)
    1874.12149028
    """
    dvap = _eq_tp_vap(temp,pres,dvap=dvap,chkvals=chkvals,chktol=chktol,
        dvap0=dvap0,chkbnd=chkbnd,mathargs=mathargs)
    g_tt = _vap_g(2,0,temp,pres,dvap=dvap)
    cp = -temp * g_tt
    return cp

def vap_cv(temp,pres,dvap=None,chkvals=False,chktol=_CHKTOL,dvap0=None,
    chkbnd=False,mathargs=None):
    """Calculate water vapour isochoric heat capacity.
    
    Calculate the isochoric heat capacity of water vapour at temperature
    and pressure.
    
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg dvap: Water vapour density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dvap: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg dvap0: Initial guess for the water vapour density in kg/m3. If
        None (default) then `_dvap_default` is used. A string specifier
        of the method (e.g. 'crit') can also be passed. See _VAPMETHODS
        in flu3a for valid specifiers.
    :type dvap0: float or str or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Heat capacity in J/kg/K.
    :raises RuntimeWarning: If a string is passed for `dvap0` that does
        not match an available method. The default is used instead.
    :raises UserWarning: If a string is passed for `dvap0` that
        specifies a function intended for water vapour.
    :raises RuntimeWarning: If the value of dvap is more consistent with
        water vapour in the subcritical region.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> vap_cv(300.,1e3)
    1410.22845789
    """
    dvap = _eq_tp_vap(temp,pres,dvap=dvap,chkvals=chkvals,chktol=chktol,
        dvap0=dvap0,chkbnd=chkbnd,mathargs=mathargs)
    g_tt = _vap_g(2,0,temp,pres,dvap=dvap)
    g_tp = _vap_g(1,1,temp,pres,dvap=dvap)
    g_pp = _vap_g(0,2,temp,pres,dvap=dvap)
    cv = temp*(-g_tt + g_tp**2/g_pp)
    return cv

def vap_density(temp,pres,dvap=None,chkvals=False,chktol=_CHKTOL,
    dvap0=None,chkbnd=False,mathargs=None):
    """Calculate water vapour density.
    
    Calculate the density of water vapour at temperature and pressure.
    
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg dvap: Water vapour density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dvap: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg dvap0: Initial guess for the water vapour density in kg/m3. If
        None (default) then `_dvap_default` is used. A string specifier
        of the method (e.g. 'crit') can also be passed. See _VAPMETHODS
        in flu3a for valid specifiers.
    :type dvap0: float or str or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Density in kg/m3.
    :raises RuntimeWarning: If a string is passed for `dvap0` that does
        not match an available method. The default is used instead.
    :raises UserWarning: If a string is passed for `dvap0` that
        specifies a function intended for water vapour.
    :raises RuntimeWarning: If the value of dvap is more consistent with
        water vapour in the subcritical region.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    """
    dvap = _eq_tp_vap(temp,pres,dvap=dvap,chkvals=chkvals,chktol=chktol,
        dvap0=dvap0,chkbnd=chkbnd,mathargs=mathargs)
    return dvap

def vap_enthalpy(temp,pres,dvap=None,chkvals=False,chktol=_CHKTOL,
    dvap0=None,chkbnd=False,mathargs=None):
    """Calculate water vapour enthalpy.
    
    Calculate the specific enthalpy of water vapour at temperature and
    pressure.
    
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg dvap: Water vapour density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dvap: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg dvap0: Initial guess for the water vapour density in kg/m3. If
        None (default) then `_dvap_default` is used. A string specifier
        of the method (e.g. 'crit') can also be passed. See _VAPMETHODS
        in flu3a for valid specifiers.
    :type dvap0: float or str or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Enthalpy in J/kg.
    :raises RuntimeWarning: If a string is passed for `dvap0` that does
        not match an available method. The default is used instead.
    :raises UserWarning: If a string is passed for `dvap0` that
        specifies a function intended for water vapour.
    :raises RuntimeWarning: If the value of dvap is more consistent with
        water vapour in the subcritical region.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> vap_enthalpy(300.,1e3)
    2551013.47892
    """
    dvap = _eq_tp_vap(temp,pres,dvap=dvap,chkvals=chkvals,chktol=chktol,
        dvap0=dvap0,chkbnd=chkbnd,mathargs=mathargs)
    g = _vap_g(0,0,temp,pres,dvap=dvap)
    g_t = _vap_g(1,0,temp,pres,dvap=dvap)
    h = g - temp*g_t
    return h

def vap_entropy(temp,pres,dvap=None,chkvals=False,chktol=_CHKTOL,
    dvap0=None,chkbnd=False,mathargs=None):
    """Calculate water vapour entropy.
    
    Calculate the specific entropy of water vapour at temperature and
    pressure.
    
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg dvap: Water vapour density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dvap: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg dvap0: Initial guess for the water vapour density in kg/m3. If
        None (default) then `_dvap_default` is used. A string specifier
        of the method (e.g. 'crit') can also be passed. See _VAPMETHODS
        in flu3a for valid specifiers.
    :type dvap0: float or str or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Entropy in J/kg/K.
    :raises RuntimeWarning: If a string is passed for `dvap0` that does
        not match an available method. The default is used instead.
    :raises UserWarning: If a string is passed for `dvap0` that
        specifies a function intended for water vapour.
    :raises RuntimeWarning: If the value of dvap is more consistent with
        water vapour in the subcritical region.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> vap_entropy(300.,1e3)
    9103.67940087
    """
    dvap = _eq_tp_vap(temp,pres,dvap=dvap,chkvals=chkvals,chktol=chktol,
        dvap0=dvap0,chkbnd=chkbnd,mathargs=mathargs)
    g_t = _vap_g(1,0,temp,pres,dvap=dvap)
    s = -g_t
    return s

def vap_expansion(temp,pres,dvap=None,chkvals=False,chktol=_CHKTOL,
    dvap0=None,chkbnd=False,mathargs=None):
    """Calculate water vapour thermal expansion coefficient.
    
    Calculate the thermal expansion coefficient of water vapour at
    temperature and pressure.
    
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg dvap: Water vapour density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dvap: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg dvap0: Initial guess for the water vapour density in kg/m3. If
        None (default) then `_dvap_default` is used. A string specifier
        of the method (e.g. 'crit') can also be passed. See _VAPMETHODS
        in flu3a for valid specifiers.
    :type dvap0: float or str or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Expansion coefficient in 1/K.
    :raises RuntimeWarning: If a string is passed for `dvap0` that does
        not match an available method. The default is used instead.
    :raises UserWarning: If a string is passed for `dvap0` that
        specifies a function intended for water vapour.
    :raises RuntimeWarning: If the value of dvap is more consistent with
        water vapour in the subcritical region.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> vap_expansion(300.,1e3)
    3.34352010567e-03
    """
    dvap = _eq_tp_vap(temp,pres,dvap=dvap,chkvals=chkvals,chktol=chktol,
        dvap0=dvap0,chkbnd=chkbnd,mathargs=mathargs)
    g_p = _vap_g(0,1,temp,pres,dvap=dvap)
    g_tp = _vap_g(1,1,temp,pres,dvap=dvap)
    alpha = g_tp / g_p
    return alpha

def vap_gibbsenergy(temp,pres,dvap=None,chkvals=False,chktol=_CHKTOL,
    dvap0=None,chkbnd=False,mathargs=None):
    """Calculate water vapour Gibbs free energy.
    
    Calculate the specific Gibbs free energy of water vapour at
    temperature and pressure.
    
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg dvap: Water vapour density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dvap: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg dvap0: Initial guess for the water vapour density in kg/m3. If
        None (default) then `_dvap_default` is used. A string specifier
        of the method (e.g. 'crit') can also be passed. See _VAPMETHODS
        in flu3a for valid specifiers.
    :type dvap0: float or str or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Gibbs energy in J/kg.
    :raises RuntimeWarning: If a string is passed for `dvap0` that does
        not match an available method. The default is used instead.
    :raises UserWarning: If a string is passed for `dvap0` that
        specifies a function intended for water vapour.
    :raises RuntimeWarning: If the value of dvap is more consistent with
        water vapour in the subcritical region.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> vap_gibbsenergy(300.,1e3)
    -180090.341338
    """
    dvap = _eq_tp_vap(temp,pres,dvap=dvap,chkvals=chkvals,chktol=chktol,
        dvap0=dvap0,chkbnd=chkbnd,mathargs=mathargs)
    g = _vap_g(0,0,temp,pres,dvap=dvap)
    return g

def vap_helmholtzenergy(temp,pres,dvap=None,chkvals=False,
    chktol=_CHKTOL,dvap0=None,chkbnd=False,mathargs=None):
    """Calculate water vapour Helmholtz free energy.
    
    Calculate the specific Helmholtz free energy of water vapour at
    temperature and pressure.
    
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg dvap: Water vapour density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dvap: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg dvap0: Initial guess for the water vapour density in kg/m3. If
        None (default) then `_dvap_default` is used. A string specifier
        of the method (e.g. 'crit') can also be passed. See _VAPMETHODS
        in flu3a for valid specifiers.
    :type dvap0: float or str or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Helmholtz energy in J/kg.
    :raises RuntimeWarning: If a string is passed for `dvap0` that does
        not match an available method. The default is used instead.
    :raises UserWarning: If a string is passed for `dvap0` that
        specifies a function intended for water vapour.
    :raises RuntimeWarning: If the value of dvap is more consistent with
        water vapour in the subcritical region.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    """
    dvap = _eq_tp_vap(temp,pres,dvap=dvap,chkvals=chkvals,chktol=chktol,
        dvap0=dvap0,chkbnd=chkbnd,mathargs=mathargs)
    g = _vap_g(0,0,temp,pres,dvap=dvap)
    g_p = _vap_g(0,1,temp,pres,dvap=dvap)
    f = g - pres*g_p
    return f

def vap_internalenergy(temp,pres,dvap=None,chkvals=False,chktol=_CHKTOL,
    dvap0=None,chkbnd=False,mathargs=None):
    """Calculate water vapour internal energy.
    
    Calculate the specific internal energy of water vapour at
    temperature and pressure.
    
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg dvap: Water vapour density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dvap: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg dvap0: Initial guess for the water vapour density in kg/m3. If
        None (default) then `_dvap_default` is used. A string specifier
        of the method (e.g. 'crit') can also be passed. See _VAPMETHODS
        in flu3a for valid specifiers.
    :type dvap0: float or str or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Internal energy in J/kg.
    :raises RuntimeWarning: If a string is passed for `dvap0` that does
        not match an available method. The default is used instead.
    :raises UserWarning: If a string is passed for `dvap0` that
        specifies a function intended for water vapour.
    :raises RuntimeWarning: If the value of dvap is more consistent with
        water vapour in the subcritical region.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> vap_internalenergy(300.,1e3)
    2412625.00085
    """
    dvap = _eq_tp_vap(temp,pres,dvap=dvap,chkvals=chkvals,chktol=chktol,
        dvap0=dvap0,chkbnd=chkbnd,mathargs=mathargs)
    g = _vap_g(0,0,temp,pres,dvap=dvap)
    g_t = _vap_g(1,0,temp,pres,dvap=dvap)
    g_p = _vap_g(0,1,temp,pres,dvap=dvap)
    u = g - temp*g_t - pres*g_p
    return u

def vap_kappa_s(temp,pres,dvap=None,chkvals=False,chktol=_CHKTOL,
    dvap0=None,chkbnd=False,mathargs=None):
    """Calculate water vapour isentropic compressibility.
    
    Calculate the isentropic compressibility of water vapour at
    temperature and pressure.
    
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg dvap: Water vapour density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dvap: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg dvap0: Initial guess for the water vapour density in kg/m3. If
        None (default) then `_dvap_default` is used. A string specifier
        of the method (e.g. 'crit') can also be passed. See _VAPMETHODS
        in flu3a for valid specifiers.
    :type dvap0: float or str or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Compressibility in 1/Pa.
    :raises RuntimeWarning: If a string is passed for `dvap0` that does
        not match an available method. The default is used instead.
    :raises UserWarning: If a string is passed for `dvap0` that
        specifies a function intended for water vapour.
    :raises RuntimeWarning: If the value of dvap is more consistent with
        water vapour in the subcritical region.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> vap_kappa_s(300.,1e3)
    7.52840457971e-04
    """
    dvap = _eq_tp_vap(temp,pres,dvap=dvap,chkvals=chkvals,chktol=chktol,
        dvap0=dvap0,chkbnd=chkbnd,mathargs=mathargs)
    g_tt = _vap_g(2,0,temp,pres,dvap=dvap)
    g_tp = _vap_g(1,1,temp,pres,dvap=dvap)
    g_pp = _vap_g(0,2,temp,pres,dvap=dvap)
    kappa = dvap*(-g_pp - g_tp**2/(-g_tt))
    return kappa

def vap_kappa_t(temp,pres,dvap=None,chkvals=False,chktol=_CHKTOL,
    dvap0=None,chkbnd=False,mathargs=None):
    """Calculate water vapour isothermal compressibility.
    
    Calculate the isothermal compressibility of water vapour at
    temperature and pressure.
    
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg dvap: Water vapour density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dvap: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg dvap0: Initial guess for the water vapour density in kg/m3. If
        None (default) then `_dvap_default` is used. A string specifier
        of the method (e.g. 'crit') can also be passed. See _VAPMETHODS
        in flu3a for valid specifiers.
    :type dvap0: float or str or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Compressibility in 1/Pa.
    :raises RuntimeWarning: If a string is passed for `dvap0` that does
        not match an available method. The default is used instead.
    :raises UserWarning: If a string is passed for `dvap0` that
        specifies a function intended for water vapour.
    :raises RuntimeWarning: If the value of dvap is more consistent with
        water vapour in the subcritical region.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> vap_kappa_t(300.,1e3)
    1.00048646242e-03
    """
    dvap = _eq_tp_vap(temp,pres,dvap=dvap,chkvals=chkvals,chktol=chktol,
        dvap0=dvap0,chkbnd=chkbnd,mathargs=mathargs)
    g_p = _vap_g(0,1,temp,pres,dvap=dvap)
    g_pp = _vap_g(0,2,temp,pres,dvap=dvap)
    kappa = -g_pp / g_p
    return kappa

def vap_lapserate(temp,pres,dvap=None,chkvals=False,chktol=_CHKTOL,
    dvap0=None,chkbnd=False,mathargs=None):
    """Calculate water vapour adiabatic lapse rate.
    
    Calculate the adiabatic lapse rate of water vapour at temperature
    and pressure.
    
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg dvap: Water vapour density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dvap: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg dvap0: Initial guess for the water vapour density in kg/m3. If
        None (default) then `_dvap_default` is used. A string specifier
        of the method (e.g. 'crit') can also be passed. See _VAPMETHODS
        in flu3a for valid specifiers.
    :type dvap0: float or str or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Lapse rate in K/Pa.
    :raises RuntimeWarning: If a string is passed for `dvap0` that does
        not match an available method. The default is used instead.
    :raises UserWarning: If a string is passed for `dvap0` that
        specifies a function intended for water vapour.
    :raises RuntimeWarning: If the value of dvap is more consistent with
        water vapour in the subcritical region.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> vap_lapserate(300.,1e3)
    7.40674488635e-02
    """
    dvap = _eq_tp_vap(temp,pres,dvap=dvap,chkvals=chkvals,chktol=chktol,
        dvap0=dvap0,chkbnd=chkbnd,mathargs=mathargs)
    g_tt = _vap_g(2,0,temp,pres,dvap=dvap)
    g_tp = _vap_g(1,1,temp,pres,dvap=dvap)
    gamma = g_tp/(-g_tt)
    return gamma

def vap_soundspeed(temp,pres,dvap=None,chkvals=False,chktol=_CHKTOL,
    dvap0=None,chkbnd=False,mathargs=None):
    """Calculate water vapour sound speed.
    
    Calculate the speed of sound in water vapour at temperature and
    pressure.
    
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg dvap: Water vapour density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dvap: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg dvap0: Initial guess for the water vapour density in kg/m3. If
        None (default) then `_dvap_default` is used. A string specifier
        of the method (e.g. 'crit') can also be passed. See _VAPMETHODS
        in flu3a for valid specifiers.
    :type dvap0: float or str or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Sound speed in m/s.
    :raises RuntimeWarning: If a string is passed for `dvap0` that does
        not match an available method. The default is used instead.
    :raises UserWarning: If a string is passed for `dvap0` that
        specifies a function intended for water vapour.
    :raises RuntimeWarning: If the value of dvap is more consistent with
        water vapour in the subcritical region.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> vap_soundspeed(300.,1e3)
    428.744430495
    """
    dvap = _eq_tp_vap(temp,pres,dvap=dvap,chkvals=chkvals,chktol=chktol,
        dvap0=dvap0,chkbnd=chkbnd,mathargs=mathargs)
    g_p = _vap_g(0,1,temp,pres,dvap=dvap)
    g_tt = _vap_g(2,0,temp,pres,dvap=dvap)
    g_tp = _vap_g(1,1,temp,pres,dvap=dvap)
    g_pp = _vap_g(0,2,temp,pres,dvap=dvap)
    dvdp = g_pp + g_tp**2/(-g_tt)
    csqinv = -dvap**2 * dvdp
    c = csqinv**(-.5)
    return c

