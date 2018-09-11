"""Humid air properties from the Gibbs energy.

This module provides thermodynamic properties (e.g. heat capacity) of
humid air based on the Gibbs free energy, with dry air mass fraction,
temperature, and pressure as primary variables. This module can also be
called as a function,

    python air_3b.py

which will compare results from this module to reference values in Lemmon et al. (2000).

:Examples:

>>> contraction(0.9,300.,1e5)
0.530280645260
>>> cp(0.9,300.,1e5)
1246.34991644
>>> density(0.9,300.,1e5)
1.09708772444

:Functions:

* compressibility: Deviation of humid air from ideal gas behavior.
* compressibility_lemmon: Deviation of humid air from ideal gas
    behavior; uses the Lemmon et al. gas constant.
* contraction: Humid air contraction coefficient.
* cp: Humid air isobaric heat capacity.
* cv: Humid air isochoric heat capacity.
* density: Humid air density.
* enthalpy: Humid air enthalpy.
* entropy: Humid air entropy.
* expansion: Humid air thermal expansion coefficient.
* gibbsenergy: Humid air Gibbs free energy.
* internalenergy: Humid air internal energy.
* kappas: Humid air isentropic compressibility.
* kappat: Humid air isothermal compressibility.
* lapserate: Humid air adiabatic lapse rate.
* soundspeed: Humid air sound speed.
* vappot: Water vapour chemical potential.
* chklemmon2000: Check module against Lemmon et al. 2000.

"""

__all__ = ['compressibility','compressibility_lemmon','contraction','cp','cv',
    'density','enthalpy','entropy','expansion','gibbsenergy','internalenergy',
    'kappas','kappat','lapserate','soundspeed','vappot',
    'chklemmon2000']

import constants0
import convert0
import air1
import air2
import air3a

_CHKTOL = constants0.CHKTOL
_eq_atp = air3a.eq_atp
_air_g = air3a.air_g
_LEMMONTOL = 3e-4  # Special tolerance level for the Lemmon tests


## Thermodynamic quantities
def compressibility(airf,temp,pres,dhum=None,chkvals=False,
    chktol=_CHKTOL,dhum0=None,chkbnd=False,mathargs=None):
    """Calculate humid air compressibility factor.
    
    Calculate the compressibility factor of humid air, a measure of the
    deviation from ideal gas behavior, from dry air mass fraction,
    temperature, and pressure. Uses the standard value for the universal
    gas constant.
    
    :arg float airf: Dry air mass fraction in kg/kg.
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg dhum: Humid air density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dhum: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg dhum0: Initial guess for the humid air density in kg/m3. If
        None (default) then `_approx_atp` is used.
    :type dhum0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Compressibility factor, unitless.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> compressibility(0.9,300.,1e5)
    0.997830231485
    """
    dhum = _eq_atp(airf,temp,pres,dhum=dhum,chkvals=chkvals,chktol=chktol,
        dhum0=dhum0,chkbnd=chkbnd,mathargs=mathargs)
    
    rav = constants0.GAS_CONSTANT_MOLAR_SI  # Standard constant
    rav /= convert0.air_molarmass(airf)
    compress = pres / (rav*dhum*temp)
    return compress

def compressibility_lemmon(airf,temp,pres,dhum=None,chkvals=False,
    chktol=_CHKTOL,dhum0=None,chkbnd=False,mathargs=None):
    """Calculate humid air compressibility with Lemmon gas constant.
    
    Calculate the compressibility factor of humid air, a measure of the
    deviation from ideal gas behavior, from dry air mass fraction,
    temperature, and pressure. Uses the Lemmon et al. (2000) value for
    the universal gas constant.
    
    :arg float airf: Dry air mass fraction in kg/kg.
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg dhum: Humid air density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dhum: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg dhum0: Initial guess for the humid air density in kg/m3. If
        None (default) then `_approx_atp` is used.
    :type dhum0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Compressibility factor, unitless.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> compressibility_lemmon(0.9,300.,1e5)
    0.997825670991
    """
    dhum = _eq_atp(airf,temp,pres,dhum=dhum,chkvals=chkvals,chktol=chktol,
        dhum0=dhum0,chkbnd=chkbnd,mathargs=mathargs)
    
    rav = constants0.GAS_CONSTANT_MOLAR_L2000  # Lemmon constant
    rav /= convert0.air_molarmass(airf)
    compress = pres / (rav*dhum*temp)
    return compress

def contraction(airf,temp,pres,dhum=None,chkvals=False,chktol=_CHKTOL,
    dhum0=None,chkbnd=False,mathargs=None):
    """Calculate humid air contraction coefficient.
    
    Calculate the contraction coefficient of humid air, the relative
    increase in density that would occur by replacing water vapour with
    dry air, from the dry air mass fraction, temperature, and pressure.
    
    :arg float airf: Dry air mass fraction in kg/kg.
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg dhum: Humid air density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dhum: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg dhum0: Initial guess for the humid air density in kg/m3. If
        None (default) then `_approx_atp` is used.
    :type dhum0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Contraction coefficient in 1/(kg/kg) (unitless).
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> contraction(0.9,300.,1e5)
    0.530280645260
    """
    dhum = _eq_atp(airf,temp,pres,dhum=dhum,chkvals=chkvals,chktol=chktol,
        dhum0=dhum0,chkbnd=chkbnd,mathargs=mathargs)
    
    g_p = _air_g(0,0,1,airf,temp,pres,dhum=dhum)
    g_ap = _air_g(1,0,1,airf,temp,pres,dhum=dhum)
    contract = -g_ap / g_p
    return contract

def cp(airf,temp,pres,dhum=None,chkvals=False,chktol=_CHKTOL,dhum0=None,
    chkbnd=False,mathargs=None):
    """Calculate humid air isobaric heat capacity.
    
    Calculate the isobaric heat capacity of humid air from dry air mass
    fraction, temperature, and pressure.
    
    :arg float airf: Dry air mass fraction in kg/kg.
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg dhum: Humid air density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dhum: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg dhum0: Initial guess for the humid air density in kg/m3. If
        None (default) then `_approx_atp` is used.
    :type dhum0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Heat capacity in J/kg/K.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> cp(0.9,300.,1e5)
    1246.34991644
    """
    dhum = _eq_atp(airf,temp,pres,dhum=dhum,chkvals=chkvals,chktol=chktol,
        dhum0=dhum0,chkbnd=chkbnd,mathargs=mathargs)
    
    g_tt = _air_g(0,2,0,airf,temp,pres,dhum=dhum)
    cp = -temp * g_tt
    return cp

def cv(airf,temp,pres,dhum=None,chkvals=False,chktol=_CHKTOL,dhum0=None,
    chkbnd=False,mathargs=None):
    """Calculate humid air isochoric heat capacity.
    
    Calculate the isochoric (constant volume) heat capacity of humid air
    from dry air mass fraction, temperature, and pressure.
    
    :arg float airf: Dry air mass fraction in kg/kg.
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg dhum: Humid air density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dhum: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg dhum0: Initial guess for the humid air density in kg/m3. If
        None (default) then `_approx_atp` is used.
    :type dhum0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Heat capacity in J/kg/K.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> cv(0.9,300.,1e5)
    920.600781012
    """
    dhum = _eq_atp(airf,temp,pres,dhum=dhum,chkvals=chkvals,chktol=chktol,
        dhum0=dhum0,chkbnd=chkbnd,mathargs=mathargs)
    
    '''
    # Use the air2 function
    cv = air2.cv(airf,temp,dhum)
    '''
    
    # Use Gibbs energy function
    g_tt = _air_g(0,2,0,airf,temp,pres,dhum=dhum)
    g_tp = _air_g(0,1,1,airf,temp,pres,dhum=dhum)
    g_pp = _air_g(0,0,2,airf,temp,pres,dhum=dhum)
    cv = -temp*g_tt + temp*g_tp**2/g_pp
    return cv

def density(airf,temp,pres,dhum=None,chkvals=False,chktol=_CHKTOL,
    dhum0=None,chkbnd=False,mathargs=None):
    """Calculate humid air density.
    
    Calculate the density of humid air from dry air mass fraction,
    temperature, and pressure.
    
    :arg float airf: Dry air mass fraction in kg/kg.
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg dhum: Humid air density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dhum: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg dhum0: Initial guess for the humid air density in kg/m3. If
        None (default) then `_approx_atp` is used.
    :type dhum0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Density in kg/m3.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> density(0.9,300.,1e5)
    1.09708772444
    """
    dhum = _eq_atp(airf,temp,pres,dhum=dhum,chkvals=chkvals,chktol=chktol,
        dhum0=dhum0,chkbnd=chkbnd,mathargs=mathargs)
    return dhum

def enthalpy(airf,temp,pres,dhum=None,chkvals=False,chktol=_CHKTOL,
    dhum0=None,chkbnd=False,mathargs=None):
    """Calculate humid air enthalpy.
    
    Calculate the specific enthalpy of humid air from dry air mass
    fraction, temperature, and pressure.
    
    :arg float airf: Dry air mass fraction in kg/kg.
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg dhum: Humid air density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dhum: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg dhum0: Initial guess for the humid air density in kg/m3. If
        None (default) then `_approx_atp` is used.
    :type dhum0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Enthalpy in J/kg.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> enthalpy(0.9,300.,1e5)
    277928.954795
    """
    dhum = _eq_atp(airf,temp,pres,dhum=dhum,chkvals=chkvals,chktol=chktol,
        dhum0=dhum0,chkbnd=chkbnd,mathargs=mathargs)
    
    '''
    # Use the air2 function
    h = air2.enthalpy(airf,temp,dhum)
    '''
    
    # Use Gibbs energy function
    g = _air_g(0,0,0,airf,temp,pres,dhum=dhum)
    g_t = _air_g(0,1,0,airf,temp,pres,dhum=dhum)
    h = g - temp*g_t
    return h

def entropy(airf,temp,pres,dhum=None,chkvals=False,chktol=_CHKTOL,
    dhum0=None,chkbnd=False,mathargs=None):
    """Calculate humid air entropy.
    
    Calculate the specific entropy of humid air from dry air mass
    fraction, temperature, and pressure.
    
    :arg float airf: Dry air mass fraction in kg/kg.
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg dhum: Humid air density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dhum: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg dhum0: Initial guess for the humid air density in kg/m3. If
        None (default) then `_approx_atp` is used.
    :type dhum0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Entropy in J/kg/K.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> entropy(0.9,300.,1e5)
    911.170080461
    """
    dhum = _eq_atp(airf,temp,pres,dhum=dhum,chkvals=chkvals,chktol=chktol,
        dhum0=dhum0,chkbnd=chkbnd,mathargs=mathargs)
    
    g_t = _air_g(0,1,0,airf,temp,pres,dhum=dhum)
    s = -g_t
    return s

def expansion(airf,temp,pres,dhum=None,chkvals=False,chktol=_CHKTOL,
    dhum0=None,chkbnd=False,mathargs=None):
    """Calculate humid air thermal expansion coefficient.
    
    Calculate the thermal expansion coefficient of humid air from the
    dry air mass fraction, temperature, and pressure.
    
    :arg float airf: Dry air mass fraction in kg/kg.
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg dhum: Humid air density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dhum: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg dhum0: Initial guess for the humid air density in kg/m3. If
        None (default) then `_approx_atp` is used.
    :type dhum0: float or None
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
    
    >>> expansion(0.9,300.,1e5)
    3.45704654420e-03
    """
    dhum = _eq_atp(airf,temp,pres,dhum=dhum,chkvals=chkvals,chktol=chktol,
        dhum0=dhum0,chkbnd=chkbnd,mathargs=mathargs)
    
    g_p = _air_g(0,0,1,airf,temp,pres,dhum=dhum)
    g_tp = _air_g(0,1,1,airf,temp,pres,dhum=dhum)
    alpha = g_tp / g_p
    return alpha

def gibbsenergy(airf,temp,pres,dhum=None,chkvals=False,chktol=_CHKTOL,
    dhum0=None,chkbnd=False,mathargs=None):
    """Calculate humid air Gibbs energy.
    
    Calculate the specific Gibbs free energy of humid air from the dry
    air mass fraction, temperature, and pressure.
    
    :arg float airf: Dry air mass fraction in kg/kg.
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg dhum: Humid air density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dhum: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg dhum0: Initial guess for the humid air density in kg/m3. If
        None (default) then `_approx_atp` is used.
    :type dhum0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Gibbs energy in J/kg.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> gibbsenergy(0.9,300.,1e5)
    4577.93065689
    """
    dhum = _eq_atp(airf,temp,pres,dhum=dhum,chkvals=chkvals,chktol=chktol,
        dhum0=dhum0,chkbnd=chkbnd,mathargs=mathargs)
    
    g = _air_g(0,0,0,airf,temp,pres,dhum=dhum)
    return g

def internalenergy(airf,temp,pres,dhum=None,chkvals=False,
    chktol=_CHKTOL,dhum0=None,chkbnd=False,mathargs=None):
    """Calculate humid air internal energy.
    
    Calculate the specific internal energy of humid air from the dry air
    mass fraction, temperature, and pressure.
    
    :arg float airf: Dry air mass fraction in kg/kg.
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg dhum: Humid air density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dhum: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg dhum0: Initial guess for the humid air density in kg/m3. If
        None (default) then `_approx_atp` is used.
    :type dhum0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Internal energy in J/kg.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> internalenergy(0.9,300.,1e5)
    186778.541048
    """
    dhum = _eq_atp(airf,temp,pres,dhum=dhum,chkvals=chkvals,chktol=chktol,
        dhum0=dhum0,chkbnd=chkbnd,mathargs=mathargs)
    
    '''
    # Use air2 function
    u = air2.air_f_internal_energy(airf,temp,dhum,chkbnd=chkbnd)
    '''
    
    # Use Gibbs energy function
    g = _air_g(0,0,0,airf,temp,pres,dhum=dhum)
    g_t = _air_g(0,1,0,airf,temp,pres,dhum=dhum)
    g_p = _air_g(0,0,1,airf,temp,pres,dhum=dhum)
    u = g - temp*g_t - pres*g_p
    return u

def kappas(airf,temp,pres,dhum=None,chkvals=False,chktol=_CHKTOL,
    dhum0=None,chkbnd=False,mathargs=None):
    """Calculate humid air isentropic compressibility.
    
    Calculate the isentropic compressibility of humid air from the dry
    air mass fraction, temperature, and pressure.
    
    :arg float airf: Dry air mass fraction in kg/kg.
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg dhum: Humid air density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dhum: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg dhum0: Initial guess for the humid air density in kg/m3. If
        None (default) then `_approx_atp` is used.
    :type dhum0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Compressibility in 1/Pa.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> kappas(0.9,300.,1e5)
    7.41034505449e-06
    """
    dhum = _eq_atp(airf,temp,pres,dhum=dhum,chkvals=chkvals,chktol=chktol,
        dhum0=dhum0,chkbnd=chkbnd,mathargs=mathargs)
    
    '''
    # Use air2 function
    kappa = air2.air_f_kappa_s(airf,temp,dhum,chkbnd=chkbnd)
    '''
    
    # Use Gibbs energy function
    g_p = _air_g(0,0,1,airf,temp,pres,dhum=dhum)
    g_tt = _air_g(0,2,0,airf,temp,pres,dhum=dhum)
    g_tp = _air_g(0,1,1,airf,temp,pres,dhum=dhum)
    g_pp = _air_g(0,0,2,airf,temp,pres,dhum=dhum)
    kappa = -(g_pp - g_tp**2/g_tt) / g_p
    return kappa

def kappat(airf,temp,pres,dhum=None,chkvals=False,chktol=_CHKTOL,
    dhum0=None,chkbnd=False,mathargs=None):
    """Calculate humid air isothermal compressibility.
    
    Calculate the isothermal compressibility of humid air from the dry
    air mass fraction, temperature, and pressure.
    
    :arg float airf: Dry air mass fraction in kg/kg.
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg dhum: Humid air density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dhum: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg dhum0: Initial guess for the humid air density in kg/m3. If
        None (default) then `_approx_atp` is used.
    :type dhum0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Compressibility in 1/Pa.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> kappat(0.9,300.,1e5)
    1.00324517749e-05
    """
    dhum = _eq_atp(airf,temp,pres,dhum=dhum,chkvals=chkvals,chktol=chktol,
        dhum0=dhum0,chkbnd=chkbnd,mathargs=mathargs)
    
    g_p = _air_g(0,0,1,airf,temp,pres,dhum=dhum)
    g_pp = _air_g(0,0,2,airf,temp,pres,dhum=dhum)
    kappa = -g_pp / g_p
    return kappa

def lapserate(airf,temp,pres,dhum=None,chkvals=False,chktol=_CHKTOL,
    dhum0=None,chkbnd=False,mathargs=None):
    """Calculate humid air adiabatic lapse rate.
    
    Calculate the dry adiabatic lapse rate of humid air from the dry air
    mass fraction, temperature, and pressure. 'Dry' here means that
    water vapour does not condense into liquid or ice; see
    :func:`liqair4b.lapserate` and :func:`iceair4b.lapserate` for those
    lapse rates.
    
    :arg float airf: Dry air mass fraction in kg/kg.
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg dhum: Humid air density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dhum: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg dhum0: Initial guess for the humid air density in kg/m3. If
        None (default) then `_approx_atp` is used.
    :type dhum0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Lapse rate in K/Pa.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> lapserate(0.9,300.,1e5)
    7.58481752251e-04
    """
    dhum = _eq_atp(airf,temp,pres,dhum=dhum,chkvals=chkvals,chktol=chktol,
        dhum0=dhum0,chkbnd=chkbnd,mathargs=mathargs)
    
    '''
    # Use air2 function
    gamma = air2.air_f_lapserate(airf,temp,dhum,chkbnd=False)
    '''
    
    # Use Gibbs energy function
    g_p = _air_g(0,0,1,airf,temp,pres,dhum=dhum)
    g_tt = _air_g(0,2,0,airf,temp,pres,dhum=dhum)
    g_tp = _air_g(0,1,1,airf,temp,pres,dhum=dhum)
    g_pp = _air_g(0,0,2,airf,temp,pres,dhum=dhum)
    gamma = -g_tp / g_tt
    return gamma

def soundspeed(airf,temp,pres,dhum=None,chkvals=False,chktol=_CHKTOL,
    dhum0=None,chkbnd=False,mathargs=None):
    """Calculate humid air sound speed.
    
    Calculate the speed of sound in humid air from the dry air mass
    fraction, temperature, and pressure.
    
    :arg float airf: Dry air mass fraction in kg/kg.
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg dhum: Humid air density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dhum: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg dhum0: Initial guess for the humid air density in kg/m3. If
        None (default) then `_approx_atp` is used.
    :type dhum0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Sound speed in m/s.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> soundspeed(0.9,300.,1e5)
    350.719656182
    """
    dhum = _eq_atp(airf,temp,pres,dhum=dhum,chkvals=chkvals,chktol=chktol,
        dhum0=dhum0,chkbnd=chkbnd,mathargs=mathargs)
    
    '''
    # Use air2 function
    c = air2.air_f_soundspeed(airf,temp,dhum,chkbnd=False)
    '''
    
    # Use Gibbs energy function
    g_p = _air_g(0,0,1,airf,temp,pres,dhum=dhum)
    g_tt = _air_g(0,2,0,airf,temp,pres,dhum=dhum)
    g_tp = _air_g(0,1,1,airf,temp,pres,dhum=dhum)
    g_pp = _air_g(0,0,2,airf,temp,pres,dhum=dhum)
    drhodp = -dhum**2 * (g_pp - g_tp**2/g_tt)
    c = drhodp**(-.5)
    return c

def vappot(airf,temp,pres,dhum=None,chkvals=False,chktol=_CHKTOL,
    dhum0=None,chkbnd=False,mathargs=None):
    """Calculate water vapour chemical potential.
    
    Calculate the chemical potential of water vapour in humid air from
    dry air mass fraction, temperature, and pressure.
    
    :arg float airf: Dry air mass fraction in kg/kg.
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg dhum: Humid air density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dhum: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg dhum0: Initial guess for the humid air density in kg/m3. If
        None (default) then `_approx_atp` is used.
    :type dhum0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Chemical potential in J/kg.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> vappot(0.9,300.,1e5)
    193705.688576
    """
    dhum = _eq_atp(airf,temp,pres,dhum=dhum,chkvals=chkvals,chktol=chktol,
        dhum0=dhum0,chkbnd=chkbnd,mathargs=mathargs)
    
    g = _air_g(0,0,0,airf,temp,pres,dhum=dhum)
    g_a = _air_g(1,0,0,airf,temp,pres,dhum=dhum)
    gvap = g - airf*g_a
    return gvap


## Functions to check results
def chklemmon2000(mode,printresult=True,chktol=_LEMMONTOL):
    """Check accuracy against Lemmon et al. (2000).
    
    Evaluate the functions in this module and compare to reference
    values from Lemmon et al. (2000). The reference values are for
    thermodynamic properties calculated from either the Helmholtz or
    Gibbs free energies, as specified by the mode.
    
    :arg int mode: Either 1 to use the Helmholtz energy functions from
        :mod:`air2` or 2 for the Gibbs energy functions (this module).
    :arg bool printresult: If True (default) and any results are outside
        of the given tolerance, then the function name, reference value,
        result value, and relative error are printed.
    :arg float chktol: Tolerance to use when choosing to print results
        (default _LEMMONTOL). The default tolerance is lowered due to
        the low number of significant figures for the reference values.
    :returns: :class:`~tester.Tester` instance containing the
        functions, arguments, reference values, results, and relative
        errors from the tests.
    :rtype: Tester
    :raises ValueError: If `mode` is not 1 or 2.
    """
    if mode not in (1,2):
        errmsg = ('Mode {0} must be either 1 (for Helmholtz functions) or 2 '
            '(for Gibbs functions)').format(mode)
        raise ValueError(errmsg)
    from tester import Tester
    MDRY = constants0.MDRY
    
    if mode==1:
        fargs_scaled = [
            (1., 140., 0.087718), (1., 270., 0.045164), (1., 400., 0.030461),
            (1.,2000., 0.006092), (1., 140.,19.810   ), (1., 270., 4.6064  ),
            (1., 400., 2.9202  ), (1.,2000., 0.59094 ), (1., 270.,47.327   ),
            (1., 400.,45.208   ), (1.,2000.,32.893   )
        ]
        fargscales = [1.,1.,1e3*MDRY]
        fargs = [tuple(a*scl for (a,scl) in zip(args,fargscales))
            for args in fargs_scaled]
        argfmt = '({1:4g},{2:10.4e})'
        funs = [air2.internalenergy,air2.enthalpy,air2.entropy,air2.cv,air2.cp,
            air2.soundspeed]
    else:
        fargs = [
            (1., 140.,101325.), (1., 270.,101325.), (1., 400.,101325.),
            (1.,2000.,101325.), (1., 140.,    1e7), (1., 270.,    1e7),
            (1., 400.,    1e7), (1.,2000.,    1e7), (1., 270.,    2e9),
            (1., 400.,    2e9), (1.,2000.,    2e9)
        ]
        argfmt = '({1:4g},{2:6g})'
        funs = [internalenergy,enthalpy,entropy,cv,cp,soundspeed]
    fnames = ['internalenergy','enthalpy','entropy','cv','cp','soundspeed']
    
    refs_scaled = [
        [2873.2,5578.9,8294.3,48610.0,-329.46,4911.3,7923.0,48600.0,4354.8,
            8076.2,53433.0],
        [4028.3,7822.4,11621.0,65242.0,175.34,7082.2,11347.0,65522.0,46614.0,
            52316.0,114240.0],
        [176.60,195.78,207.26,259.62,118.30,155.25,168.19,221.44,96.232,
            113.52,176.45],
        [20.81,20.76,21.04,27.90,25.82,21.61,21.38,27.93,37.64,34.45,31.86],
        [29.38,29.13,29.38,36.21,72.88,35.28,31.50,36.25,45.74,42.27,38.21],
        [236.4,329.6,400.5,863.5,418.1,349.7,425.6,878.6,2899.8,2822.9,2472.1]
    ]
    refscales = [MDRY**(-1)]*5 + [1.]
    
    # Need to add offsets from Lemmon et al. (2000)
    TRED = air1._TRED_DRY
    RDRY_L2000 = constants0.GAS_CONSTANT_AIR_L2000
    N4_CURR, N5_CURR = air1._C_DRYF0[0]
    N4_L2000 = -13.841928076
    N5_L2000 = 17.275266575
    S0 = RDRY_L2000 * (N4_CURR - N4_L2000)
    F0 = RDRY_L2000*TRED * (N5_L2000 - N5_CURR)
    refoffsets = [F0,F0,S0,0.,0.,0.]
    refs = [[(r*scl-off) for r in ref]
        for (ref,scl,off) in zip(refs_scaled,refscales,refoffsets)]
    
    # TODO: Incorporate equilibrium + kwargs into Tester
    test = Tester(funs,fargs,refs,fnames,argfmt)
    test.run()
    if printresult:
        if mode==1:
            msg = 'Dry air Helmholtz functions'
        else:
            msg = 'Dry air Gibbs functions'
        print(msg)
        test.printresults(chktol=chktol)
    return test


## Main function: Check tables
if __name__ == '__main__':
    testf = chklemmon2000(1)
    testg = chklemmon2000(2)

