"""Seawater Gibbs free energy and related properties.

This module provides the Gibbs free energy of seawater (liquid water and
salt) and its derivatives with respect to salinity, temperature, and
pressure. It also provides properties (e.g. heat capacity) derived from
the Gibbs energy. Finally, this module can be called as a function::

    python sea3a.py

which will compare results from this module to reference values in IAPWS
2008, tables 8a, 8b, and 8c.

:Examples:

>>> sea_g(0,0,0,0.035,300.,1e5)
-5114.99619856
>>> sea_g(0,0,1,0.035,300.,1e5)
9.77858615182e-04
>>> sea_g(0,2,0,0.035,300.,1e5)
-13.3358324655
>>> density(0.035,300.,1e5)
1022.64272613
>>> lapserate(0.035,300.,1e5)
2.28411342567e-08
>>> temp_maxdensity(0.01,1e5)
274.950121498

:Functions:

* :func:`sea_g`: Seawater Gibbs free energy with derivatives.
* :func:`liqpot`: Chemical potential of liquid water in seawater.
* :func:`salpot`: Chemical potential of salt in seawater.
* :func:`contraction_t`: Haline contraction coefficient of seawater at
  constant in-situ temperature (isothermal).
* :func:`cp`: Seawater isobaric heat capacity.
* :func:`density`: Seawater density.
* :func:`enthalpy`: Seawater enthalpy.
* :func:`helmholtzenergy`: Seawater Helmholtz free energy.
* :func:`entropy`: Seawater entropy.
* :func:`expansion_t`: Seawater thermal expansion coefficient with
  respect to in-situ temperature.
* :func:`gibbsenergy`: Seawater Gibbs free energy.
* :func:`internalenergy`: Seawater internal energy.
* :func:`kappa_s`: Seawater isentropic compressibility.
* :func:`kappa_t`: Seawater isothermal compressibility.
* :func:`lapserate`: Seawater adiabatic lapse rate.
* :func:`osmcoeff`: Seawater osmotic coefficient.
* :func:`soundspeed`: Seawater sound speed.
* :func:`temp_maxdensity`: Temperature of maximum seawater density.
* :func:`chkiapws08table8a`: Check module against IAPWS 2008, table 8a.
* :func:`chkiapws08table8b`: Check module against IAPWS 2008, table 8b.
* :func:`chkiapws08table8c`: Check module against IAPWS 2008, table 8c.

"""

__all__ = ['sea_g','liqpot','salpot','contraction_t','cp','density','enthalpy',
    'helmholtzenergy','entropy','expansion_t','gibbsenergy','internalenergy',
    'kappa_s','kappa_t','lapserate','osmcoeff','soundspeed','temp_maxdensity',
    'chkiapws08table8a','chkiapws08table8b','chkiapws08table8c']

import warnings
import constants0
import flu2
import sal2
import maths3
import flu3a

_CHKTOL = constants0.CHKTOL
_chkflubnds = constants0.chkflubnds
_chksalbnds = constants0.chksalbnds
_sal_g = sal2.sal_g
_newton = maths3.newton
_eq_tp_liq = flu3a.eq_tp_liq
_liq_g = flu3a.liq_g
_MDT = 277.
_MDD = 1e3

## Gibbs energy
def sea_g(drvs,drvt,drvp,salt,temp,pres,dliq=None,chkvals=False,
    chktol=_CHKTOL,dliq0=None,chkbnd=False,useext=False,mathargs=None):
    """Calculate seawater Gibbs free energy with derivatives.
    
    Calculate the specific Gibbs free energy of seawater or its
    derivatives with respect to salinity, temperature, and pressure.
    
    :arg int drvs: Number of salinity derivatives.
    :arg int drvt: Number of temperature derivatives.
    :arg int drvp: Number of pressure derivatives.
    :arg float salt: Salinity in kg/kg.
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
    :arg bool useext: If False (default) then the salt contribution is
        calculated from _GSCOEFFS; if True, from _GSCOEFFS_EXT.
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Gibbs energy in units of
        (J/kg) / (kg/kg)^drvs / K^drvt / Pa^drvp.
    :raises RuntimeWarning: If a string is passed for `dliq0` that does
        not match an available method. The default is used instead.
    :raises UserWarning: If a string is passed for `dliq0` that
        specifies a function intended for water vapour.
    :raises RuntimeWarning: If the value of dliq is more consistent with
        water vapour in the subcritical region.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> sea_g(0,0,0,0.035,300.,1e5)
    -5114.99619856
    >>> sea_g(0,0,1,0.035,300.,1e5)
    9.77858615182e-04
    >>> sea_g(0,2,0,0.035,300.,1e5)
    -13.3358324655
    """
    dliq = _eq_tp_liq(temp,pres,dliq=dliq,chkvals=chkvals,chktol=chktol,
        dliq0=dliq0,mathargs=mathargs)
    _chkflubnds(temp,dliq,chkbnd=chkbnd)
    _chksalbnds(salt,temp,pres,chkbnd=chkbnd)
    g = _sal_g(drvs,drvt,drvp,salt,temp,pres,useext=useext)
    if drvs == 0:
        g += _liq_g(drvt,drvp,temp,pres,dliq=dliq)
    return g


## Thermodynamic properties
def liqpot(salt,temp,pres,dliq=None,chkvals=False,chktol=_CHKTOL,
    dliq0=None,chkbnd=False,useext=False,mathargs=None):
    """Calculate liquid water chemical potential.
    
    Calculate the chemical potential of liquid water in seawater at
    salinity, temperature, and pressure.
    
    :arg float salt: Salinity in kg/kg.
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
    :arg bool useext: If False (default) then the salt contribution is
        calculated from _GSCOEFFS; if True, from _GSCOEFFS_EXT.
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Chemical potential in J/kg.
    :raises RuntimeWarning: If a string is passed for `dliq0` that does
        not match an available method. The default is used instead.
    :raises UserWarning: If a string is passed for `dliq0` that
        specifies a function intended for water vapour.
    :raises RuntimeWarning: If the value of dliq is more consistent with
        water vapour in the subcritical region.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> liqpot(0.035,300.,1e5)
    -7865.77834936
    """
    dliq = _eq_tp_liq(temp,pres,dliq=dliq,chkvals=chkvals,chktol=chktol,
        dliq0=dliq0,mathargs=mathargs)
    _chkflubnds(temp,dliq,chkbnd=chkbnd)
    _chksalbnds(salt,temp,pres,chkbnd=chkbnd)
    g = sea_g(0,0,0,salt,temp,pres,dliq=dliq,useext=useext)
    g_s = sea_g(1,0,0,salt,temp,pres,dliq=dliq,useext=useext)
    lpot = g - salt*g_s
    return lpot

def salpot(salt,temp,pres,dliq=None,chkvals=False,chktol=_CHKTOL,
    dliq0=None,chkbnd=False,useext=False,mathargs=None):
    """Calculate salt chemical potential.
    
    Calculate the chemical potential of salt in seawater at salinity,
    temperature, and pressure.
    
    :arg float salt: Salinity in kg/kg.
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
    :arg bool useext: If False (default) then the salt contribution is
        calculated from _GSCOEFFS; if True, from _GSCOEFFS_EXT.
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Chemical potential in J/kg.
    :raises RuntimeWarning: If a string is passed for `dliq0` that does
        not match an available method. The default is used instead.
    :raises UserWarning: If a string is passed for `dliq0` that
        specifies a function intended for water vapour.
    :raises RuntimeWarning: If the value of dliq is more consistent with
        water vapour in the subcritical region.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> salpot(0.035,300.,1e5)
    78593.7757371
    """
    _chksalbnds(salt,temp,pres,chkbnd=chkbnd)
    gsal = sal2.salpot(salt,temp,pres,useext=useext)
    return gsal

def contraction_t(salt,temp,pres,dliq=None,chkvals=False,chktol=_CHKTOL,
    dliq0=None,chkbnd=False,useext=False,mathargs=None):
    """Calculate seawater haline contraction coefficient.
    
    Calculate the haline contraction coefficient of seawater at constant
    in-situ temperature (isothermal) at salinity, temperature, and
    pressure.
    
    :arg float salt: Salinity in kg/kg.
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
    :arg bool useext: If False (default) then the salt contribution is
        calculated from _GSCOEFFS; if True, from _GSCOEFFS_EXT.
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Contraction coefficient in kg/kg.
    :raises RuntimeWarning: If a string is passed for `dliq0` that does
        not match an available method. The default is used instead.
    :raises UserWarning: If a string is passed for `dliq0` that
        specifies a function intended for water vapour.
    :raises RuntimeWarning: If the value of dliq is more consistent with
        water vapour in the subcritical region.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> contraction_t(0.035,300.,1e5)
    0.732910044599
    """
    dliq = _eq_tp_liq(temp,pres,dliq=dliq,chkvals=chkvals,chktol=chktol,
        dliq0=dliq0,mathargs=mathargs)
    _chkflubnds(temp,dliq,chkbnd=chkbnd)
    _chksalbnds(salt,temp,pres,chkbnd=chkbnd)
    g_p = sea_g(0,0,1,salt,temp,pres,dliq=dliq,useext=useext)
    g_sp = sea_g(1,0,1,salt,temp,pres,dliq=dliq,useext=useext)
    beta = -g_sp / g_p
    return beta

def cp(salt,temp,pres,dliq=None,chkvals=False,chktol=_CHKTOL,dliq0=None,
    chkbnd=False,useext=False,mathargs=None):
    """Calculate seawater isobaric heat capacity.
    
    Calculate the isobaric (constant pressure) heat capacity of seawater
    at salinity, temperature, and pressure.
    
    :arg float salt: Salinity in kg/kg.
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
    :arg bool useext: If False (default) then the salt contribution is
        calculated from _GSCOEFFS; if True, from _GSCOEFFS_EXT.
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
    
    >>> cp(0.035,300.,1e5)
    4000.74973964
    """
    dliq = _eq_tp_liq(temp,pres,dliq=dliq,chkvals=chkvals,chktol=chktol,
        dliq0=dliq0,mathargs=mathargs)
    _chkflubnds(temp,dliq,chkbnd=chkbnd)
    _chksalbnds(salt,temp,pres,chkbnd=chkbnd)
    g_tt = sea_g(0,2,0,salt,temp,pres,dliq=dliq,useext=useext)
    cp = -temp * g_tt
    return cp

def density(salt,temp,pres,dliq=None,chkvals=False,chktol=_CHKTOL,
    dliq0=None,chkbnd=False,useext=False,mathargs=None):
    """Calculate seawater density.
    
    Calculate the density of seawater at salinity, temperature, and
    pressure.
    
    :arg float salt: Salinity in kg/kg.
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
    :arg bool useext: If False (default) then the salt contribution is
        calculated from _GSCOEFFS; if True, from _GSCOEFFS_EXT.
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
    
    :Examples:
    
    >>> density(0.035,300.,1e5)
    1022.64272613
    """
    dliq = _eq_tp_liq(temp,pres,dliq=dliq,chkvals=chkvals,chktol=chktol,
        dliq0=dliq0,mathargs=mathargs)
    _chkflubnds(temp,dliq,chkbnd=chkbnd)
    _chksalbnds(salt,temp,pres,chkbnd=chkbnd)
    g_p = sea_g(0,0,1,salt,temp,pres,dliq=dliq,useext=useext)
    dsea = g_p**(-1)
    return dsea

def enthalpy(salt,temp,pres,dliq=None,chkvals=False,chktol=_CHKTOL,
    dliq0=None,chkbnd=False,useext=False,mathargs=None):
    """Calculate seawater enthalpy.
    
    Calculate the specific enthalpy of seawater at salinity,
    temperature, and pressure.
    
    :arg float salt: Salinity in kg/kg.
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
    :arg bool useext: If False (default) then the salt contribution is
        calculated from _GSCOEFFS; if True, from _GSCOEFFS_EXT.
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
    
    >>> enthalpy(0.035,300.,1e5)
    107220.675963
    """
    dliq = _eq_tp_liq(temp,pres,dliq=dliq,chkvals=chkvals,chktol=chktol,
        dliq0=dliq0,mathargs=mathargs)
    _chkflubnds(temp,dliq,chkbnd=chkbnd)
    _chksalbnds(salt,temp,pres,chkbnd=chkbnd)
    g = sea_g(0,0,0,salt,temp,pres,dliq=dliq,useext=useext)
    g_t = sea_g(0,1,0,salt,temp,pres,dliq=dliq,useext=useext)
    h = g - temp*g_t
    return h

def entropy(salt,temp,pres,dliq=None,chkvals=False,chktol=_CHKTOL,
    dliq0=None,chkbnd=False,useext=False,mathargs=None):
    """Calculate seawater entropy.
    
    Calculate the specific entropy of seawater at salinity, temperature,
    and pressure.
    
    :arg float salt: Salinity in kg/kg.
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
    :arg bool useext: If False (default) then the salt contribution is
        calculated from _GSCOEFFS; if True, from _GSCOEFFS_EXT.
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
    
    >>> entropy(0.035,300.,1e5)
    374.452240540
    """
    dliq = _eq_tp_liq(temp,pres,dliq=dliq,chkvals=chkvals,chktol=chktol,
        dliq0=dliq0,mathargs=mathargs)
    _chkflubnds(temp,dliq,chkbnd=chkbnd)
    _chksalbnds(salt,temp,pres,chkbnd=chkbnd)
    g_t = sea_g(0,1,0,salt,temp,pres,dliq=dliq,useext=useext)
    s = -g_t
    return s

def expansion_t(salt,temp,pres,dliq=None,chkvals=False,chktol=_CHKTOL,
    dliq0=None,chkbnd=False,useext=False,mathargs=None):
    """Calculate seawater thermal expansion coefficient.
    
    Calculate the thermal expansion coefficient of seawater with respect
    to in-situ temperature at salinity, temperature, and pressure.
    
    :arg float salt: Salinity in kg/kg.
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
    :arg bool useext: If False (default) then the salt contribution is
        calculated from _GSCOEFFS; if True, from _GSCOEFFS_EXT.
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
    
    >>> expansion_t(0.035,300.,1e5)
    3.11502639583e-04
    """
    dliq = _eq_tp_liq(temp,pres,dliq=dliq,chkvals=chkvals,chktol=chktol,
        dliq0=dliq0,mathargs=mathargs)
    _chkflubnds(temp,dliq,chkbnd=chkbnd)
    _chksalbnds(salt,temp,pres,chkbnd=chkbnd)
    g_p = sea_g(0,0,1,salt,temp,pres,dliq=dliq,useext=useext)
    g_tp = sea_g(0,1,1,salt,temp,pres,dliq=dliq,useext=useext)
    alpha = g_tp / g_p
    return alpha

def gibbsenergy(salt,temp,pres,dliq=None,chkvals=False,chktol=_CHKTOL,
    dliq0=None,chkbnd=False,useext=False,mathargs=None):
    """Calculate seawater Gibbs free energy.
    
    Calculate the specific Gibbs free energy of seawater at salinity,
    temperature, and pressure.
    
    :arg float salt: Salinity in kg/kg.
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
    :arg bool useext: If False (default) then the salt contribution is
        calculated from _GSCOEFFS; if True, from _GSCOEFFS_EXT.
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
    
    >>> gibbsenergy(0.035,300.,1e5)
    -5114.99619856
    """
    dliq = _eq_tp_liq(temp,pres,dliq=dliq,chkvals=chkvals,chktol=chktol,
        dliq0=dliq0,mathargs=mathargs)
    _chkflubnds(temp,dliq,chkbnd=chkbnd)
    _chksalbnds(salt,temp,pres,chkbnd=chkbnd)
    g = sea_g(0,0,0,salt,temp,pres,dliq=dliq,useext=useext)
    return g

def helmholtzenergy(salt,temp,pres,dliq=None,chkvals=False,
    chktol=_CHKTOL,dliq0=None,chkbnd=False,useext=False,mathargs=None):
    """Calculate seawater Helmholtz free energy.
    
    Calculate the specific Helmholtz free energy of seawater at
    salinity, temperature, and pressure.
    
    :arg float salt: Salinity in kg/kg.
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
    :arg bool useext: If False (default) then the salt contribution is
        calculated from _GSCOEFFS; if True, from _GSCOEFFS_EXT.
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
        dliq0=dliq0,mathargs=mathargs)
    _chkflubnds(temp,dliq,chkbnd=chkbnd)
    _chksalbnds(salt,temp,pres,chkbnd=chkbnd)
    g = sea_g(0,0,0,salt,temp,pres,dliq=dliq,useext=useext)
    g_p = sea_g(0,0,1,salt,temp,pres,dliq=dliq,useext=useext)
    f = g - pres*g_p
    return f

def internalenergy(salt,temp,pres,dliq=None,chkvals=False,
    chktol=_CHKTOL,dliq0=None,chkbnd=False,useext=False,mathargs=None):
    """Calculate seawater internal energy.
    
    Calculate the specific internal energy of seawater at salinity,
    temperature, and pressure.
    
    :arg float salt: Salinity in kg/kg.
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
    :arg bool useext: If False (default) then the salt contribution is
        calculated from _GSCOEFFS; if True, from _GSCOEFFS_EXT.
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
    
    >>> internalenergy(0.035,300.,1e5)
    107122.890102
    """
    dliq = _eq_tp_liq(temp,pres,dliq=dliq,chkvals=chkvals,chktol=chktol,
        dliq0=dliq0,mathargs=mathargs)
    _chkflubnds(temp,dliq,chkbnd=chkbnd)
    _chksalbnds(salt,temp,pres,chkbnd=chkbnd)
    g = sea_g(0,0,0,salt,temp,pres,dliq=dliq,useext=useext)
    g_t = sea_g(0,1,0,salt,temp,pres,dliq=dliq,useext=useext)
    g_p = sea_g(0,0,1,salt,temp,pres,dliq=dliq,useext=useext)
    u = g - temp*g_t - pres*g_p
    return u

def kappa_s(salt,temp,pres,dliq=None,chkvals=False,chktol=_CHKTOL,
    dliq0=None,chkbnd=False,useext=False,mathargs=None):
    """Calculate seawater isentropic compressibility.
    
    Calculate the isentropic compressibility of seawater at salinity,
    temperature, and pressure.
    
    :arg float salt: Salinity in kg/kg.
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
    :arg bool useext: If False (default) then the salt contribution is
        calculated from _GSCOEFFS; if True, from _GSCOEFFS_EXT.
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
    
    >>> kappa_s(0.035,300.,1e5)
    4.13135667732e-10
    """
    dliq = _eq_tp_liq(temp,pres,dliq=dliq,chkvals=chkvals,chktol=chktol,
        dliq0=dliq0,mathargs=mathargs)
    _chkflubnds(temp,dliq,chkbnd=chkbnd)
    _chksalbnds(salt,temp,pres,chkbnd=chkbnd)
    g_p = sea_g(0,0,1,salt,temp,pres,dliq=dliq,useext=useext)
    g_tt = sea_g(0,2,0,salt,temp,pres,dliq=dliq,useext=useext)
    g_tp = sea_g(0,1,1,salt,temp,pres,dliq=dliq,useext=useext)
    g_pp = sea_g(0,0,2,salt,temp,pres,dliq=dliq,useext=useext)
    kappa = (g_tp**2/g_tt - g_pp) / g_p
    return kappa

def kappa_t(salt,temp,pres,dliq=None,chkvals=False,chktol=_CHKTOL,
    dliq0=None,chkbnd=False,useext=False,mathargs=None):
    """Calculate seawater isothermal compressibility.
    
    Calculate the isothermal (constant in-situ temperature)
    compressibility of seawater at salinity, temperature, and pressure.
    
    :arg float salt: Salinity in kg/kg.
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
    :arg bool useext: If False (default) then the salt contribution is
        calculated from _GSCOEFFS; if True, from _GSCOEFFS_EXT.
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
    
    >>> kappa_t(0.035,300.,1e5)
    4.20250741344e-10
    """
    dliq = _eq_tp_liq(temp,pres,dliq=dliq,chkvals=chkvals,chktol=chktol,
        dliq0=dliq0,mathargs=mathargs)
    _chkflubnds(temp,dliq,chkbnd=chkbnd)
    _chksalbnds(salt,temp,pres,chkbnd=chkbnd)
    g_p = sea_g(0,0,1,salt,temp,pres,dliq=dliq,useext=useext)
    g_pp = sea_g(0,0,2,salt,temp,pres,dliq=dliq,useext=useext)
    kappat = -g_pp / g_p
    return kappat

def lapserate(salt,temp,pres,dliq=None,chkvals=False,chktol=_CHKTOL,
    dliq0=None,chkbnd=False,useext=False,mathargs=None):
    """Calculate seawater adiabatic lapse rate.
    
    Calculate the adiabatic lapse rate of seawater at salinity,
    temperature, and pressure.
    
    :arg float salt: Salinity in kg/kg.
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
    :arg bool useext: If False (default) then the salt contribution is
        calculated from _GSCOEFFS; if True, from _GSCOEFFS_EXT.
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
    
    >>> lapserate(0.035,300.,1e5)
    2.28411342567e-08
    """
    dliq = _eq_tp_liq(temp,pres,dliq=dliq,chkvals=chkvals,chktol=chktol,
        dliq0=dliq0,mathargs=mathargs)
    _chkflubnds(temp,dliq,chkbnd=chkbnd)
    _chksalbnds(salt,temp,pres,chkbnd=chkbnd)
    g_tt = sea_g(0,2,0,salt,temp,pres,dliq=dliq,useext=useext)
    g_tp = sea_g(0,1,1,salt,temp,pres,dliq=dliq,useext=useext)
    gamma = -g_tp / g_tt
    return gamma

def osmcoeff(salt,temp,pres,dliq=None,chkvals=False,chktol=_CHKTOL,
    dliq0=None,chkbnd=False,useext=False,mathargs=None):
    """Calculate seawater osmotic coefficient.
    
    Calculate the osmotic coefficient of salt in seawater at salinity,
    temperature, and pressure.
    
    :arg float salt: Salinity in kg/kg.
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
    :arg bool useext: If False (default) then the salt contribution is
        calculated from _GSCOEFFS; if True, from _GSCOEFFS_EXT.
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Osmotic coefficient, unitless.
    :raises RuntimeWarning: If a string is passed for `dliq0` that does
        not match an available method. The default is used instead.
    :raises UserWarning: If a string is passed for `dliq0` that
        specifies a function intended for water vapour.
    :raises RuntimeWarning: If the value of dliq is more consistent with
        water vapour in the subcritical region.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> osmcoeff(0.035,300.,1e5)
    0.902777495349
    """
    _chksalbnds(salt,temp,pres,chkbnd=chkbnd)
    phi = sal2.osmcoeff(salt,temp,pres,useext=useext)
    return phi

def soundspeed(salt,temp,pres,dliq=None,chkvals=False,chktol=_CHKTOL,
    dliq0=None,chkbnd=False,useext=False,mathargs=None):
    """Calculate seawater sound speed.
    
    Calculate the speed of sound in seawater at salinity, temperature,
    and pressure.
    
    :arg float salt: Salinity in kg/kg.
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
    :arg bool useext: If False (default) then the salt contribution is
        calculated from _GSCOEFFS; if True, from _GSCOEFFS_EXT.
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
    
    >>> soundspeed(0.035,300.,1e5)
    1538.47940766
    """
    dliq = _eq_tp_liq(temp,pres,dliq=dliq,chkvals=chkvals,chktol=chktol,
        dliq0=dliq0,mathargs=mathargs)
    _chkflubnds(temp,dliq,chkbnd=chkbnd)
    _chksalbnds(salt,temp,pres,chkbnd=chkbnd)
    g_p = sea_g(0,0,1,salt,temp,pres,dliq=dliq,useext=useext)
    g_tt = sea_g(0,2,0,salt,temp,pres,dliq=dliq,useext=useext)
    g_tp = sea_g(0,1,1,salt,temp,pres,dliq=dliq,useext=useext)
    g_pp = sea_g(0,0,2,salt,temp,pres,dliq=dliq,useext=useext)
    csqinv = (g_tp**2/g_tt - g_pp) / g_p**2
    c = csqinv**(-.5)
    return c


## Temperature of maximum density
def _volfun(td,salt,pres,useext=False):
    """Calculate seawater specific volume at TD/SP.
    
    Calculate the specific volume of seawater from temperature and
    liquid water density for fixed salinity and pressure. This function
    is formatted for use in temp_maxdensity, specifically by
    scipy.optimize.fmin_slsqp.
    
    :arg td: Temperature in K and liquid water density in kg/m3.
    :type td: list[float,float]
    :arg float salt: Salinity in kg/kg.
    :arg float pres: Pressure in Pa.
    :arg bool useext: If False (default) then the salt contribution is
        calculated from _GSCOEFFS; if True, from _GSCOEFFS_EXT.
    :returns: Specific volume in m3/kg.
    :rtype: float
    """
    temp, dliq = td[:]
    gl_p = dliq**(-1)
    gs_p = _sal_g(0,0,1,salt,temp,pres,useext=useext)
    vol = gl_p + gs_p
    return vol

def _volder(td,salt,pres,chkbnd=False,useext=False):
    """Calculate seawater specific volume Jacobian at TD/SP.
    
    Calculate the derivatives of the specific volume of seawater with
    respect to temperature and liquid water density for fixed salinity
    and pressure. This function is formatted for use in temp_maxdensity,
    specifically by scipy.optimize.fmin_slsqp.
    
    :arg td: Temperature in K and liquid water density in kg/m3.
    :type td: list[float,float]
    :arg float salt: Salinity in kg/kg.
    :arg float pres: Pressure in Pa.
    :arg bool useext: If False (default) then the salt contribution is
        calculated from _GSCOEFFS; if True, from _GSCOEFFS_EXT.
    :returns: Volume derivatives in m3/kg/K and (m3/kg)/(kg/m3).
    :rtype: list[float,float]
    """
    temp, dliq = td[:]
    gl_pd = -dliq**(-2)
    gs_tp = _sal_g(0,1,1,salt,temp,pres,useext=useext)
    dvol = [gs_tp, gl_pd]
    return dvol

def _pdiffun(td,pres):
    """Calculate seawater disequilibrium at TD.
    
    Calculate the difference (pres-pliq) between the given pressure and
    the pressure of liquid water at the given temperature and density.
    This function is formatted for use in temp_maxdensity, specifically
    by scipy.optimize.fmin_slsqp.
    
    :arg td: Temperature in K and liquid water density in kg/m3.
    :type td: list[float,float]
    :arg float salt: Salinity in kg/kg.
    :arg float pres: Pressure in Pa.
    :arg bool useext: If False (default) then the salt contribution is
        calculated from _GSCOEFFS; if True, from _GSCOEFFS_EXT.
    :returns: Pressure difference in Pa.
    :rtype: float
    """
    temp, dliq = td[:]
    pliq = flu2.eq_pressure(0,0,temp,dliq)
    pdif = pliq - pres
    return pdif

def _pdifder(td,pres):
    """Calculate seawater disequilibrium derivative at TD.
    
    Calculate the derivatives of the difference (pres-pliq) with respect
    to temperature and liquid water density, where pliq is the pressure
    in liquid water. This function is formatted for use in
    temp_maxdensity, specifically by scipy.optimize.fmin_slsqp.
    
    :arg td: Temperature in K and liquid water density in kg/m3.
    :type td: list[float,float]
    :arg float salt: Salinity in kg/kg.
    :arg float pres: Pressure in Pa.
    :arg bool useext: If False (default) then the salt contribution is
        calculated from _GSCOEFFS; if True, from _GSCOEFFS_EXT.
    :returns: Pressude derivatives in Pa/K and Pa/(kg/m3).
    :rtype: list[float,float]
    """
    temp, dliq = td[:]
    pliq_t = flu2.eq_pressure(1,0,temp,dliq)
    pliq_d = flu2.eq_pressure(0,1,temp,dliq)
    dpdif = [pliq_t, pliq_d]
    return dpdif

def temp_maxdensity(salt,pres,temp0=_MDT,dliq0=_MDD,chkbnd=False,
    useext=False,mathargs=None):
    """Calculate the temperature of maximum seawater density.
    
    Calculate the temperature at which seawater at the given salinity
    and pressure reaches maximum density. This function uses
    :func:`~scipy.optimize.fmin_slsqp` for optimization.
    
    :arg float salt: Salinity in kg/kg.
    :arg float pres: Pressure in Pa.
    :arg float temp0: Initial guess for the temperature in K
        (default _MDT).
    :arg float dliq0: Initial guess for the liquid water density in
        kg/m3 (default _MDD).
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg bool useext: If False (default) then the salt contribution is
        calculated from _GSCOEFFS; if True, from _GSCOEFFS_EXT.
    :arg mathargs: Keyword arguments to the root-finder
        :func:`~scipy.optimize.fmin_slsqp`. If None (default) then no
        arguments are passed and default parameters will be used. If
        `full_output` is passed and True, then temp_maxdensity will also
        return the liquid water density and all the summary information
        from `fmin_slsqp`.
    :returns: Temperature of maximum density in K. If `full_output` is
        passed in `mathargs` and is True, then the liquid water density
        will also be returned along with the summary information from
        `fmin_slsqp`.
    :rtype: float or tuple(float,float,tuple).
    :raises ImportError: If scipy.optimize is not accessible.
    
    :Examples:
    
    >>> temp_maxdensity(0.01,1e5)
    274.950121498
    """
    try:
        import scipy.optimize
    except ImportError:
        errmsg = ('Scipy is currently required for sea_temp_maxdens')
        raise ImportError(errmsg)
    
    if mathargs is None:
        mathargs = {'disp': 0}
        full_output = False
    else:
        full_output = mathargs.get('full_output',False)
        if 'disp' not in mathargs.keys():
            mathargs['disp'] = 0
    
    # Cast the problem as constrained minimization
    f = lambda x: _volfun(x,salt,pres,useext=useext)
    fprime = lambda x: _volder(x,salt,pres,useext=useext)
    f_eqcons = lambda x: _pdiffun(x,pres)
    fprime_eqcons = lambda x: _pdifder(x,pres)
    x0 = [temp0, dliq0]
    res = scipy.optimize.fmin_slsqp(f,x0,f_eqcons=f_eqcons,fprime=fprime,
        fprime_eqcons=fprime_eqcons,**mathargs)
    
    if full_output:
        t_maxd, dl_maxd = res[0][:]
        summary = res[1:]
        return (t_maxd, dl_maxd, summary)
    t_maxd, dl_maxd = res[:]
    return t_maxd


## Functions to check results
def chkiapws08table8a(printresult=True,chktol=_CHKTOL):
    """Check accuracy against IAPWS 2008 table 8a.
    
    Evaluate the functions in this module and compare to reference
    values from IAPWS (2008), table 8a. These tests are for the
    properties of salt, liquid water, and seawater at standard salinity,
    0 Celsius, and 1 atm.
    
    :arg bool printresult: If True (default) and any results are outside
        of the given tolerance, then the function name, reference value,
        result value, and relative error are printed.
    :arg float chktol: Tolerance to use when choosing to print results
        (default _LEMMONTOL). The default tolerance is lowered due to
        the low number of significant figures for the reference values.
    :returns: :class:`~tester.Tester` instances containing the
        functions, arguments, reference values, results, and relative
        errors from the tests. The first three instances are for the
        Gibbs energy functions of liquid water, salt, and seawater; the
        other three instances are for other thermodynamic functions for
        liquid water, salt, and seawater.
    """
    from tester import Tester
    import flu3b
    args1 = (0.03516504, 273.15, 101325.)
    DERS2 = [(0,0),(1,0),(0,1),(2,0),(1,1),(0,2)]
    DERS3 = [(0,0,0),(1,0,0),(0,1,0),(0,0,1),(1,0,1),(0,2,0),(0,1,1),(0,0,2)]
    
    funs = _liq_g
    fargs = [(der+args1[1:]) for der in DERS2]
    refs = [1.01342742e2,1.47643376e-1,1.00015693912169e-3,-1.54473542320e1,
        -6.777003179e-8,-5.08928894643e-13]
    fnames = 'liq_g'
    argfmt = '({0:1d},{1:1d},{2:6.2f},{3:6g})'
    header = 'Liquid Gibbs derivatives'
    eqfun = _eq_tp_liq
    eqargs = args1[1:]
    eqkeys = ['dliq']
    testliqder = Tester(funs,fargs,refs,fnames,argfmt,header=header,
        eqfun=eqfun,eqargs=eqargs,eqkeys=eqkeys)
    
    funs = _sal_g
    fargs = [(der+args1) for der in DERS3]
    refs = [-1.0134274172939e2,6.39974067312299e4,-1.47643376346e-1,
        -2.74957224268433e-5,-7.59615411515309e-4,8.5286115117592e-1,
        1.1928678741396e-7,5.8153517233288e-14]
    fnames = 'sal_g'
    argfmt = '({0:1d},{1:1d},{2:1d},{3:10.8f},{4:6.2f},{5:6g})'
    header = 'Salt Gibbs derivatives'
    testsalder = Tester(funs,fargs,refs,fnames,argfmt,header=header)
    
    funs = sea_g
    fargs = [(der+args1) for der in DERS3]
    refs = [4.e-9,6.39974067312299e4,-6.e-11,9.7266121669485e-4,
        -7.59615411515309e-4,-1.459449308080e1,5.1516755627e-8,
        -4.507753774102e-13]
    fnames = 'sea_g'
    argfmt = '({0:1d},{1:1d},{2:1d},{3:10.8f},{4:6.2f},{5:6g})'
    header = 'Sea Gibbs derivatives'
    eqfun = _eq_tp_liq
    eqargs = args1[1:]
    eqkeys = ['dliq']
    testseader = Tester(funs,fargs,refs,fnames,argfmt,header=header,
        eqfun=eqfun,eqargs=eqargs,eqkeys=eqkeys)
    
    funs = [flu3b.liq_enthalpy,flu3b.liq_helmholtzenergy,
        flu3b.liq_internalenergy,flu3b.liq_entropy,flu3b.liq_density,
        flu3b.liq_cp,flu3b.liq_soundspeed,flu3b.liq_gibbsenergy]
    fargs = args1[1:]
    refs = [6.10139535e1,1.83989364064e-3,-4.03269484e1,-1.47643376e-1,
        9.9984308550433e2,4.21944480846e3,1.40238253109e3,1.01342742e2]
    refs_alt = [None,1.8399e-3,None,None,None,None,None,None]
    fnames = ['enthalpy','helmholtzenergy','internalenergy','entropy','density',
        'cp','soundspeed','gibbsenergy']
    argfmt = '({0:6.2f},{1:6g})'
    header = 'Liquid functions'
    eqfun = _eq_tp_liq
    eqargs = args1[1:]
    eqkeys = ['dliq']
    testliqfun = Tester(funs,fargs,refs,fnames,argfmt,header=header,
        eqfun=eqfun,eqargs=eqargs,eqkeys=eqkeys)
    
    funs = [sal2._auxenthalpy,sal2._auxhelmholtzenergy,sal2._auxinternalenergy,
        sal2._auxentropy,sal2._auxcp,sal2.liqpot]
    fargs = args1
    refs = [-6.10139534804e1,-9.8556737654491e1,-5.82279494055e1,
        1.47643376346e-1,-2.3295902344370e2,-2.35181410932936e3]
    fnames = ['enthalpy','helmholtzenergy','internalenergy','entropy','cp',
        'liqpot']
    argfmt = '({0:10.8f},{1:6.2f},{2:6g})'
    header = 'Salt functions'
    testsalfun = Tester(funs,fargs,refs,fnames,argfmt,header=header)
    
    funs = [enthalpy,helmholtzenergy,internalenergy,entropy,density,cp,
        soundspeed,liqpot]
    fargs = args1
    refs = [2.e-8,-9.855489778e1,-9.85548978e1,6.e-11,1.02810719995401e3,
        3.98648578502e3,1.44900246362e3,-2.250471368e3]
    fnames = ['enthalpy','helmholtzenergy','internalenergy','entropy','density',
        'cp','soundspeed','liqpot']
    argfmt = '({0:10.8f},{1:6.2f},{2:6g})'
    header = 'Sea functions'
    eqfun = _eq_tp_liq
    eqargs = args1[1:]
    eqkeys = ['dliq']
    testseafun = Tester(funs,fargs,refs,fnames,argfmt,header=header,
        eqfun=eqfun,eqargs=eqargs,eqkeys=eqkeys)
    
    tests = (testliqder,testsalder,testseader,testliqfun,testsalfun,
        testseafun)
    for test in tests:
        test.run()
        if printresult:
            test.printresults(chktol=chktol)
    return tests

def chkiapws08table8b(printresult=True,chktol=_CHKTOL):
    """Check accuracy against IAPWS 2008 table 8b.
    
    Evaluate the functions in this module and compare to reference
    values from IAPWS (2008), table 8a. These tests are for the
    properties of salt, liquid water, and seawater at a salinity of 0.1
    kg/kg, temperature of 353 K, and pressure of 1 atm.
    
    :arg bool printresult: If True (default) and any results are outside
        of the given tolerance, then the function name, reference value,
        result value, and relative error are printed.
    :arg float chktol: Tolerance to use when choosing to print results
        (default _LEMMONTOL). The default tolerance is lowered due to
        the low number of significant figures for the reference values.
    :returns: :class:`~tester.Tester` instances containing the
        functions, arguments, reference values, results, and relative
        errors from the tests. The first three instances are for the
        Gibbs energy functions of liquid water, salt, and seawater; the
        other three instances are for other thermodynamic functions for
        liquid water, salt, and seawater.
    """
    from tester import Tester
    import flu3b
    args1 = (0.1, 353., 101325.)
    DERS2 = [(0,0),(1,0),(0,1),(2,0),(1,1),(0,2)]
    DERS3 = [(0,0,0),(1,0,0),(0,1,0),(0,0,1),(1,0,1),(0,2,0),(0,1,1),(0,0,2)]
    
    funs = _liq_g
    fargs = [(der+args1[1:]) for der in DERS2]
    refs = [-4.46114968996e4,-1.0737599318875e3,1.02892955635611e-3,
        -1.1888500004755e1,6.59051552339e-7,-4.746728193611e-13]
    fnames = 'liq_g'
    argfmt = '({0:1d},{1:1d},{2:6.2f},{3:6g})'
    header = 'Liquid Gibbs derivatives'
    eqfun = _eq_tp_liq
    eqargs = args1[1:]
    eqkeys = ['dliq']
    testliqder = Tester(funs,fargs,refs,fnames,argfmt,header=header,
        eqfun=eqfun,eqargs=eqargs,eqkeys=eqkeys)
    
    funs = _sal_g
    fargs = [(der+args1) for der in DERS3]
    refs = [1.5087174003705e4,2.51957275851413e5,1.56230907404291e2,
        -5.7922728577126e-5,-3.0595780244234e-4,1.27922649315507,
        8.0306159575153e-7,2.1308615424374e-13]
    fnames = 'sal_g'
    argfmt = '({0:1d},{1:1d},{2:1d},{3:10.8f},{4:6.2f},{5:6g})'
    header = 'Salt Gibbs derivatives'
    testsalder = Tester(funs,fargs,refs,fnames,argfmt,header=header)
    
    funs = sea_g
    fargs = [(der+args1) for der in DERS3]
    refs = [-2.95243228959e4,2.51957275851413e5,-9.175290244832e2,
        9.7100682777898e-4,-3.0595780244234e-4,-1.0609273511600e1,
        1.462113148091e-6,-2.61586665117e-13]
    fnames = 'sea_g'
    argfmt = '({0:1d},{1:1d},{2:1d},{3:10.8f},{4:6.2f},{5:6g})'
    header = 'Sea Gibbs derivatives'
    eqfun = _eq_tp_liq
    eqargs = args1[1:]
    eqkeys = ['dliq']
    testseader = Tester(funs,fargs,refs,fnames,argfmt,header=header,
        eqfun=eqfun,eqargs=eqargs,eqkeys=eqkeys)
    
    funs = [flu3b.liq_enthalpy,flu3b.liq_helmholtzenergy,
        flu3b.liq_internalenergy,flu3b.liq_entropy,flu3b.liq_density,
        flu3b.liq_cp,flu3b.liq_soundspeed,flu3b.liq_gibbsenergy]
    fargs = args1[1:]
    refs = [3.344257590567e5,-4.47157531869e4,3.343215027694e5,
        1.0737599318875e3,9.7188383191308e2,4.1966405016784e3,
        1.5544629665347e3,-4.46114968996e4]
    fnames = ['enthalpy','helmholtzenergy','internalenergy','entropy','density',
        'cp','soundspeed','gibbsenergy']
    argfmt = '({0:6.2f},{1:6g})'
    header = 'Liquid functions'
    eqfun = _eq_tp_liq
    eqargs = args1[1:]
    eqkeys = ['dliq']
    testliqfun = Tester(funs,fargs,refs,fnames,argfmt,header=header,
        eqfun=eqfun,eqargs=eqargs,eqkeys=eqkeys)
    
    funs = [sal2._auxenthalpy,sal2._auxhelmholtzenergy,sal2._auxinternalenergy,
        sal2._auxentropy,sal2._auxcp,sal2.liqpot]
    fargs = args1
    refs = [-4.006233631001e4,1.5093043024178e4,-4.0056467289536e4,
        -1.56230907404291e2,-4.51566952083741e2,-1.01085535814360e4]
    fnames = ['enthalpy','helmholtzenergy','internalenergy','entropy','cp',
        'liqpot']
    argfmt = '({0:10.8f},{1:6.2f},{2:6g})'
    header = 'Salt functions'
    testsalfun = Tester(funs,fargs,refs,fnames,argfmt,header=header)
    
    funs = [enthalpy,helmholtzenergy,internalenergy,entropy,density,cp,
        soundspeed,liqpot]
    fargs = args1
    refs = [2.94363422747e5,-2.96227101627e4,2.942650354799e5,
        9.175290244832e2,1.02985887574790e3,3.745073549595e3,
        3.9612783529e3,-5.47200504810e4]
    fnames = ['enthalpy','helmholtzenergy','internalenergy','entropy','density',
        'cp','soundspeed','liqpot']
    argfmt = '({0:10.8f},{1:6.2f},{2:6g})'
    header = 'Sea functions'
    eqfun = _eq_tp_liq
    eqargs = args1[1:]
    eqkeys = ['dliq']
    testseafun = Tester(funs,fargs,refs,fnames,argfmt,header=header,
        eqfun=eqfun,eqargs=eqargs,eqkeys=eqkeys)
    
    tests = (testliqder,testsalder,testseader,testliqfun,testsalfun,
        testseafun)
    for test in tests:
        test.run()
        if printresult:
            test.printresults(chktol=chktol)
    return tests

def chkiapws08table8c(printresult=True,chktol=_CHKTOL):
    """Check accuracy against IAPWS 2008 table 8c.
    
    Evaluate the functions in this module and compare to reference
    values from IAPWS (2008), table 8c. These tests are for the
    properties of salt, liquid water, and seawater at standard salinity,
    0 Celsius, and a pressure of 1e8 Pa.
    
    :arg bool printresult: If True (default) and any results are outside
        of the given tolerance, then the function name, reference value,
        result value, and relative error are printed.
    :arg float chktol: Tolerance to use when choosing to print results
        (default _LEMMONTOL). The default tolerance is lowered due to
        the low number of significant figures for the reference values.
    :returns: :class:`~tester.Tester` instances containing the
        functions, arguments, reference values, results, and relative
        errors from the tests. The first three instances are for the
        Gibbs energy functions of liquid water, salt, and seawater; the
        other three instances are for other thermodynamic functions for
        liquid water, salt, and seawater.
    """
    from tester import Tester
    import flu3b
    args1 = (0.03516504, 273.15, 1e8)
    DERS2 = [(0,0),(1,0),(0,1),(2,0),(1,1),(0,2)]
    DERS3 = [(0,0,0),(1,0,0),(0,1,0),(0,0,1),(1,0,1),(0,2,0),(0,1,1),(0,0,2)]
    
    funs = _liq_g
    fargs = [(der+args1[1:]) for der in DERS2]
    refs = [9.773038621954e4,8.5146650206,9.5668332915351e-4,
        -1.429698733876e1,1.99079570803e-7,-3.715308894234e-13]
    fnames = 'liq_g'
    argfmt = '({0:1d},{1:1d},{2:6.2f},{3:6g})'
    header = 'Liquid Gibbs derivatives'
    eqfun = _eq_tp_liq
    eqargs = args1[1:]
    eqkeys = ['dliq']
    testliqder = Tester(funs,fargs,refs,fnames,argfmt,header=header,
        eqfun=eqfun,eqargs=eqargs,eqkeys=eqkeys)
    
    funs = _sal_g
    fargs = [(der+args1) for der in DERS3]
    refs = [-2.60093050730637e3,-5.4586158064880e3,7.5404568488117,
        -2.2912384179113e-5,-6.4075761854575e-4,4.88076973942251e-1,
        4.6628441224121e-8,3.57345735845327e-14]
    fnames = 'sal_g'
    argfmt = '({0:1d},{1:1d},{2:1d},{3:10.8f},{4:6.2f},{5:6g})'
    header = 'Salt Gibbs derivatives'
    testsalder = Tester(funs,fargs,refs,fnames,argfmt,header=header)
    
    funs = sea_g
    fargs = [(der+args1) for der in DERS3]
    refs = [9.512945571223e4,-5.4586158064880e3,1.60551218694e1,
        9.3377094497440e-4,-6.40757618545748e-4,-1.380891036482e1,
        2.45708012027e-7,-3.35796315839e-13]
    fnames = 'sea_g'
    argfmt = '({0:1d},{1:1d},{2:1d},{3:10.8f},{4:6.2f},{5:6g})'
    header = 'Sea Gibbs derivatives'
    eqfun = _eq_tp_liq
    eqargs = args1[1:]
    eqkeys = ['dliq']
    testseader = Tester(funs,fargs,refs,fnames,argfmt,header=header,
        eqfun=eqfun,eqargs=eqargs,eqkeys=eqkeys)
    
    funs = [flu3b.liq_enthalpy,flu3b.liq_helmholtzenergy,
        flu3b.liq_internalenergy,flu3b.liq_entropy,flu3b.liq_density,
        flu3b.liq_cp,flu3b.liq_soundspeed,flu3b.liq_gibbsenergy]
    fargs = args1[1:]
    refs = [9.54046054692e4,2.06205330419e3,-2.637274462e2,-8.5146650206,
        1.04527796139692e3,3.905222091582e3,1.575422398486e3,9.77303862195e4]
    fnames = ['enthalpy','helmholtzenergy','internalenergy','entropy','density',
        'cp','soundspeed','gibbsenergy']
    argfmt = '({0:6.2f},{1:6g})'
    header = 'Liquid functions'
    eqfun = _eq_tp_liq
    eqargs = args1[1:]
    eqkeys = ['dliq']
    testliqfun = Tester(funs,fargs,refs,fnames,argfmt,header=header,
        eqfun=eqfun,eqargs=eqargs,eqkeys=eqkeys)
    
    funs = [sal2._auxenthalpy,sal2._auxhelmholtzenergy,sal2._auxinternalenergy,
        sal2._auxentropy,sal2._auxcp,sal2.liqpot]
    fargs = args1
    refs = [-4.6606062955593e3,-3.0969208939506e2,-2.3693678776480e3,
        -7.5404568488117,-1.33318225432326e2,-2.4089780641266e3]
    fnames = ['enthalpy','helmholtzenergy','internalenergy','entropy','cp',
        'liqpot']
    argfmt = '({0:10.8f},{1:6.2f},{2:6g})'
    header = 'Salt functions'
    testsalfun = Tester(funs,fargs,refs,fnames,argfmt,header=header)
    
    funs = [enthalpy,helmholtzenergy,internalenergy,entropy,density,cp,
        soundspeed,liqpot]
    fargs = args1
    refs = [9.07439991736e4,1.75236121479e3,-2.6330953238e3,-1.6055121869e1,
        1.0709264465574e3,3.77190386615e3,1.621989976499e3,9.532140815541e4]
    fnames = ['enthalpy','helmholtzenergy','internalenergy','entropy','density',
        'cp','soundspeed','liqpot']
    argfmt = '({0:10.8f},{1:6.2f},{2:6g})'
    header = 'Sea functions'
    eqfun = _eq_tp_liq
    eqargs = args1[1:]
    eqkeys = ['dliq']
    testseafun = Tester(funs,fargs,refs,fnames,argfmt,header=header,
        eqfun=eqfun,eqargs=eqargs,eqkeys=eqkeys)
    
    tests = (testliqder,testsalder,testseader,testliqfun,testsalfun,
        testseafun)
    for test in tests:
        test.run()
        if printresult:
            test.printresults(chktol=chktol)
    return tests


## Main function: Check tables
if __name__ == '__main__':
    testsa = chkiapws08table8a();
    testsb = chkiapws08table8b();
    testsc = chkiapws08table8c();

