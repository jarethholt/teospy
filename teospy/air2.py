"""Humid air Helmholtz function and thermodynamic properties.

This module implements the Helmholtz free energy for humid air (dry air
and water vapour) as a function of dry air mass fraction, temperature,
and humid air density. It also provides functions for thermodynamic
properties calculated directly from the free energy. This module can
also be called as a function,

    python air_2.py

which will compare results from this module to the reference values in
IAPWS 2010, tables 13, 14, and 15.

Functions with the prefix `eq_` are used in higher-level modules. They
provide the given quantities or derivatives with respect to the primary
variables.

:Examples:

>>> air_f(0,0,0,0.9,300.,1.)
-95019.5943231
>>> lapserate(0.9,300.,1.)
8.50371537341e-04
>>> pressure(0.9,300.,1.)
91175.3848662

:Functions:

* air_f: Helmholtz free energy of humid air with derivatives.
* cp: Isobaric heat capacity of humid air.
* cv: Isochoric heat capacity of humid air.
* enthalpy: Specific enthalpy of humid air.
* entropy: Specific entropy of humid air.
* expansion: Thermal expansion coefficient of humid air.
* gibbsenergy: Specific Gibbs free energy of humid air.
* internalenergy: Specific internal energy of humid air.
* kappas: Isentropic compressibility of humid air.
* kappat: Isothermal compressibility of humid air.
* lapserate: Adiabatic lapse rate of humid air.
* pressure: Pressure of humid air.
* soundspeed: Speed of sound in humid air.
* vappot: Water vapour chemical potential in humid air.
* eq_entropy: Humid air entropy with derivatives.
* eq_pressure: Humid air pressure with derivatives.
* eq_vappot: Water vapour chemical potential with derivatives.
* chkiapws10table: Check module accuracy against an IAPWS 2010 reference
    table.
* chkiapws10table13: Check module accuracy against IAPWS10 table 13.
* chkiapws10table14: Check module accuracy against IAPWS10 table 14.
* chkiapws10table15: Check module accuracy against IAPWS10 table 15.

"""

__all__ = ['air_f','cp','cv','enthalpy','entropy','expansion','gibbsenergy',
    'internalenergy','kappas','kappat','lapserate','pressure','soundspeed',
    'vappot',
    'eq_entropy','eq_pressure','eq_vappot',
    'chkiapws10table','chkiapws10table13','chkiapws10table14',
    'chkiapws10table15']

import constants0
import air1
import flu1

# Single constants
_MDRY = constants0.MDRY
_MWAT = constants0.MWAT
_RUNIV = constants0.RUNIV
_CHKTOL = constants0.CHKTOL

# Auxiliary functions
_chkhumbnds = constants0.chkhumbnds
_dry_f = air1.dry_f
_air_baw = air1.air_baw
_air_caaw = air1.air_caaw
_air_caww = air1.air_caww
_flu_f = flu1.flu_f


## Helmholtz free energy functions
def _air_f_mix(drva,drvt,drvd,airf,temp,dhum,chkbnd=False):
    """Calculate free energy of mixing with derivatives.
    
    Calculate the excess free energy of mixing for humid air (dry air
    and water vapour) and its derivatives with respect to air mass
    fraction, temperature, and density. This is the free energy due to
    non-ideal air-water mixing, which is determined by the virial
    coefficients. Derivatives up to second order are available.
    
    :arg int drva: Number of dry fraction derivatives.
    :arg int drvt: Number of temperature derivatives.
    :arg int drvd: Number of density derivatives.
    :arg float airf: Dry air mass fraction in kg/kg.
    :arg float temp: Temperature in K.
    :arg float dhum: Humid air density in kg/m3.
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :returns: Free energy of mixing in units of
        (J/kg) / (kg/kg)^drva / K^drvt / (kg/m3)^drvd.
    :raises ValueError: If temp or dhum are nonpositive or if airf is
        not between 0 and 1.
    :raises RuntimeWarning: If temp or dhum are outside the recommended
        bounds and chkbnd is True.
    :raises ValueError: If any of (drva,drvt,drvd) are negative or if
        (drva+drvt+drvd)>2.
    
    :Examples:
    
    >>> air_f_mix(0,0,0,0.9,300.,1.)
    -25.8379179182
    >>> air_f_mix(1,0,0,0.9,300.,1.)
    233.827370434
    >>> air_f_mix(0,1,0,0.9,300.,1.)
    0.164195952060
    >>> air_f_mix(0,0,1,0.9,300.,1.)
    -26.2357498619
    >>> air_f_mix(2,0,0,0.9,300.,1.)
    500.273928155
    >>> air_f_mix(1,1,0,0.9,300.,1.)
    -1.53932744055
    >>> air_f_mix(1,0,1,0.9,300.,1.)
    241.520643317
    >>> air_f_mix(0,2,0,0.9,300.,1.)
    -0.687329742959e-3
    >>> air_f_mix(0,1,1,0.9,300.,1.)
     0.172192606103
    >>> air_f_mix(0,0,2,0.9,300.,1.)
    -0.795663887493
    """
    _chkhumbnds(airf,temp,dhum,chkbnd=chkbnd)
    if not (all((drv>=0)*(drv<=2) for drv in (drva,drvt,drvd))
            and (drva+drvt+drvd)<=2):
        errmsg = 'Derivatives {0} not recognized'.format((drva,drvt,drvd))
        raise ValueError(errmsg)
    
    # Calculate virial terms with temperature derivatives
    baw = _air_baw(drvt,temp) * temp
    caaw = _air_caaw(drvt,temp) * temp
    caww = _air_caww(drvt,temp) * temp
    if drvt == 0:
        pass
    elif drvt == 1:
        baw += _air_baw(0,temp)
        caaw += _air_caaw(0,temp)
        caww += _air_caww(0,temp)
    elif drvt == 2:
        baw += 2 * _air_baw(1,temp)
        caaw += 2 * _air_caaw(1,temp)
        caww += 2 * _air_caww(1,temp)
    baw *= 2 * _RUNIV / (_MDRY * _MWAT)
    caaw *= 1.5 * _RUNIV / (_MDRY * _MWAT) / _MDRY
    caww *= 1.5 * _RUNIV / (_MDRY * _MWAT) / _MWAT
    
    # Modify terms with density derivatives
    if drvd == 0:
        baw *= dhum
        caaw *= dhum**2
        caww *= dhum**2
    elif drvd == 1:
        caaw *= 2 * dhum
        caww *= 2 * dhum
    elif drvd == 2:
        baw *= 0.
        caaw *= 2
        caww *= 2
    
    # Modify terms with dry fraction derivatives
    if drva == 0:
        baw *= airf * (1 - airf)
        caaw *= airf**2 * (1 - airf)
        caww *= airf * (1 - airf)**2
    elif drva == 1:
        baw *= 1 - 2*airf
        caaw *= airf * (2 - 3*airf)
        caww *= 1 - 4*airf + 3*airf**2
    elif drva == 2:
        baw *= -2
        caaw *= 2 - 6*airf
        caww *= -4 + 6*airf
    fmix = baw + caaw + caww
    return fmix

def air_f(drva,drvt,drvd,airf,temp,dhum,chkbnd=False):
    """Calculate humid air Helmholtz free energy with derivatives.
    
    Calculate the Helmholtz free energy of humid air or its derivatives
    with respect to dry air mass fraction, temperature, or density.
    Derivatives up to second order are available.
    
    :arg int drva: Number of dry fraction derivatives.
    :arg int drvt: Number of temperature derivatives.
    :arg int drvd: Number of density derivatives.
    :arg float airf: Dry air mass fraction in kg/kg.
    :arg float temp: Temperature in K.
    :arg float dhum: Humid air density in kg/m3.
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :returns: Helmholtz free energy in units of
        (J/kg) / (kg/kg)^drva / K^drvt / (kg/m3)^drvd.
    :raises ValueError: If temp or dhum are nonpositive or if airf is
        not between 0 and 1.
    :raises RuntimeWarning: If temp or dhum are outside the recommended
        bounds and chkbnd is True.
    :raises ValueError: If any of (drva,drvt,drvd) are negative or if
        (drva+drvt+drvd)>2.
    
    :Examples:
    
    >>> air_f(0,0,0,0.9,300.,1.)
    -95019.5943231
    >>> air_f(1,0,0,0.9,300.,1.)
    -205645.554995
    >>> air_f(0,1,0,0.9,300.,1.)
    -940.175394023
    >>> air_f(0,0,1,0.9,300.,1.)
    91175.3848662
    >>> air_f(2,0,0,0.9,300.,1.)
    1447768.46379
    >>> air_f(1,1,0,0.9,300.,1.)
    7443.09771950
    >>> air_f(1,0,1,0.9,300.,1.)
    -48847.9096826
    >>> air_f(0,2,0,0.9,300.,1.)
    -2.96482218054
    >>> air_f(0,1,1,0.9,300.,1.)
    312.063110700
    >>> air_f(0,0,2,0.9,300.,1.)
    -91421.4440689
    """
    _chkhumbnds(airf,temp,dhum,chkbnd=chkbnd)
    if not (all((drv>=0)*(drv<=2) for drv in (drva,drvt,drvd))
            and (drva+drvt+drvd)<=2):
        errmsg = 'Derivatives {0} not recognized'.format((drva,drvt,drvd))
        raise ValueError(errmsg)
    f = 0.
    ddry = airf * dhum
    dvap = (1-airf) * dhum
    
    # Is water vapour present?
    fvap = 0.
    if dvap > 0:
        if drva == 0:
            fvap += (_flu_f(drvt,drvd,temp,dvap)
                * (1-airf)**(drvd+1))
        elif drva == 1:
            if drvd == 0:
                fvap -= _flu_f(drvt,0,temp,dvap)
                fvap -= dvap * _flu_f(drvt,1,temp,dvap)
            elif drvd == 1:
                fvap -= 2*(1-airf) * _flu_f(drvt,1,temp,dvap)
                fvap -= dvap*(1-airf) * _flu_f(drvt,2,temp,dvap)
        elif drva == 2:
            fvap += 2*dhum * _flu_f(drvt,1,temp,dvap)
            fvap += dvap*dhum * _flu_f(drvt,2,temp,dvap)
    
    # Is dry air present?
    fdry = 0.
    if ddry > 0:
        if drva == 0:
            fdry += (_dry_f(drvt,drvd,temp,ddry)
                * airf**(drvd+1))
        elif drva == 1:
            if drvd == 0:
                fdry += _dry_f(drvt,0,temp,ddry)
                fdry += ddry * _dry_f(drvt,1,temp,ddry)
            elif drvd == 1:
                fdry += 2*airf * _dry_f(drvt,1,temp,ddry)
                fdry += ddry*airf * _dry_f(drvt,2,temp,ddry)
        elif drva == 2:
            fdry += 2*dhum * _dry_f(drvt,1,temp,ddry)
            fdry += ddry*dhum * _dry_f(drvt,2,temp,ddry)
    
    # Are both water vapour and dry air present?
    fmix = 0.
    if (dvap > 0) and (ddry > 0):
        fmix += _air_f_mix(drva,drvt,drvd,airf,temp,dhum)
    f = fdry + fvap + fmix
    return f


## Thermodynamic properties
def cp(airf,temp,dhum,chkbnd=False):
    """Calculate humid air isobaric heat capacity.
    
    Calculate the isobaric (constant pressure) heat capacity of humid
    air from dry air mass fraction, temperature, and humid air density.
    
    :arg float airf: Dry air mass fraction in kg/kg.
    :arg float temp: Temperature in K.
    :arg float dhum: Humid air density in kg/m3.
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :returns: Heat capacity in J/kg/K.
    :raises ValueError: If temp or dhum are nonpositive or if airf is
        not between 0 and 1.
    :raises RuntimeWarning: If temp or dhum are outside the recommended
        bounds and chkbnd is True.
    
    :Examples:
    
    >>> cp(0.9,300.,1.)
    1210.74031058
    """
    _chkhumbnds(airf,temp,dhum,chkbnd=chkbnd)
    f_d = air_f(0,0,1,airf,temp,dhum)
    f_tt = air_f(0,2,0,airf,temp,dhum)
    f_td = air_f(0,1,1,airf,temp,dhum)
    f_dd = air_f(0,0,2,airf,temp,dhum)
    cp = f_td**2 - f_dd*f_tt - 2*f_d*f_tt/dhum
    cp /= f_dd + 2*f_d/dhum
    cp *= temp
    return cp

def cv(airf,temp,dhum,chkbnd=False):
    """Calculate humid air isochoric heat capacity.
    
    Calculate the isochoric (constant volume) heat capacity of humid air
    from dry air mass fraction, temperature, and humid air density.
    
    :arg float airf: Dry air mass fraction in kg/kg.
    :arg float temp: Temperature in K.
    :arg float dhum: Humid air density in kg/m3.
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :returns: Heat capacity in J/kg/K.
    :raises ValueError: If temp or dhum are nonpositive or if airf is
        not between 0 and 1.
    :raises RuntimeWarning: If temp or dhum are outside the recommended
        bounds and chkbnd is True.
    
    :Examples:
    
    >>> air_f_cv(0.9,300.,1.)
    889.446654163
    """
    _chkhumbnds(airf,temp,dhum,chkbnd=chkbnd)
    f_tt = air_f(0,2,0,airf,temp,dhum)
    cv = -temp * f_tt
    return cv

def enthalpy(airf,temp,dhum,chkbnd=False):
    """Calculate humid air enthalpy.
    
    Calculate the specific enthalpy of humid air from dry air mass
    fraction, temperature, and humid air density.
    
    :arg float airf: Dry air mass fraction in kg/kg.
    :arg float temp: Temperature in K.
    :arg float dhum: Humid air density in kg/m3.
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :returns: Enthalpy in J/kg.
    :raises ValueError: If temp or dhum are nonpositive or if airf is
        not between 0 and 1.
    :raises RuntimeWarning: If temp or dhum are outside the recommended
        bounds and chkbnd is True.
    
    :Examples:
    
    >>> enthalpy(0.9,300.,1.)
    278208.408750
    """
    _chkhumbnds(airf,temp,dhum,chkbnd=chkbnd)
    f = air_f(0,0,0,airf,temp,dhum)
    f_t = air_f(0,1,0,airf,temp,dhum)
    f_d = air_f(0,0,1,airf,temp,dhum)
    h = f - temp*f_t + dhum*f_d
    return h

def entropy(airf,temp,dhum,chkbnd=False):
    """Calculate humid air entropy.
    
    Calculate the specific entropy of humid air from dry air mass
    fraction, temperature, and humid air density.
    
    :arg float airf: Dry air mass fraction in kg/kg.
    :arg float temp: Temperature in K.
    :arg float dhum: Humid air density in kg/m3.
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :returns: Entropy in J/kg/K.
    :raises ValueError: If temp or dhum are nonpositive or if airf is
        not between 0 and 1.
    :raises RuntimeWarning: If temp or dhum are outside the recommended
        bounds and chkbnd is True.
    
    :Examples:
    
    >>> entropy(0.9,300.,1.)
    940.175394023
    """
    _chkhumbnds(airf,temp,dhum,chkbnd=chkbnd)
    f_t = air_f(0,1,0,airf,temp,dhum)
    s = -f_t
    return s

def expansion(airf,temp,dhum,chkbnd=False):
    """Calculate humid air thermal expansion coefficient.
    
    Calculate the thermal expansion coefficient of humid air from dry
    air mass fraction, temperature, and humid air density.
    
    :arg float airf: Dry air mass fraction in kg/kg.
    :arg float temp: Temperature in K.
    :arg float dhum: Humid air density in kg/m3.
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :returns: Expansion coefficient in 1/K.
    :raises ValueError: If temp or dhum are nonpositive or if airf is
        not between 0 and 1.
    :raises RuntimeWarning: If temp or dhum are outside the recommended
        bounds and chkbnd is True.
    
    :Examples:
    
    >>> expansion(0.9,300.,1.)
    3.43193033077e-03
    """
    _chkhumbnds(airf,temp,dhum,chkbnd=chkbnd)
    f_d = air_f(0,0,1,airf,temp,dhum)
    f_td = air_f(0,1,1,airf,temp,dhum)
    f_dd = air_f(0,0,2,airf,temp,dhum)
    denom = 2*f_d + dhum*f_dd
    alpha = f_td / denom
    return alpha

def gibbsenergy(airf,temp,dhum,chkbnd=False):
    """Calculate humid air Gibbs free energy.
    
    Calculate the specific Gibbs free energy of humid air from dry air
    mass fraction, temperature, and humid air density.
    
    :arg float airf: Dry air mass fraction in kg/kg.
    :arg float temp: Temperature in K.
    :arg float dhum: Humid air density in kg/m3.
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :returns: Gibbs free energy in J/kg.
    :raises ValueError: If temp or dhum are nonpositive or if airf is
        not between 0 and 1.
    :raises RuntimeWarning: If temp or dhum are outside the recommended
        bounds and chkbnd is True.
    
    :Examples:
    
    >>> gibbsenergy(0.9,300.,1.)
    -3844.20945693
    """
    _chkhumbnds(airf,temp,dhum,chkbnd=chkbnd)
    f = air_f(0,0,0,airf,temp,dhum)
    f_d = air_f(0,0,1,airf,temp,dhum)
    g = f + dhum*f_d
    return g

def internalenergy(airf,temp,dhum,chkbnd=False):
    """Calculate humid air internal energy.
    
    Calculate the specific internal energy of humid air from dry air
    mass fraction, temperature, and humid air density.
    
    :arg float airf: Dry air mass fraction in kg/kg.
    :arg float temp: Temperature in K.
    :arg float dhum: Humid air density in kg/m3.
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :returns: Internal energy in J/kg.
    :raises ValueError: If temp or dhum are nonpositive or if airf is
        not between 0 and 1.
    :raises RuntimeWarning: If temp or dhum are outside the recommended
        bounds and chkbnd is True.
    
    :Examples:
    
    >>> internalenergy(0.9,300.,1.)
    187033.023884
    """
    _chkhumbnds(airf,temp,dhum,chkbnd=chkbnd)
    f = air_f(0,0,0,airf,temp,dhum)
    f_t = air_f(0,1,0,airf,temp,dhum)
    u = f - temp*f_t
    return u

def kappas(airf,temp,dhum,chkbnd=False):
    """Calculate humid air isentropic compressibility.
    
    Calculate the isentropic compressibility of humid air from the dry
    air mass fraction, temperature, and humid air density.
    
    :arg float airf: Dry air mass fraction in kg/kg.
    :arg float temp: Temperature in K.
    :arg float dhum: Humid air density in kg/m3.
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :returns: Compressibility in 1/Pa.
    :raises ValueError: If temp or dhum are nonpositive or if airf is
        not between 0 and 1.
    :raises RuntimeWarning: If temp or dhum are outside the recommended
        bounds and chkbnd is True.
    
    :Examples:
    
    >>> kappas(0.9,300.,1.)
    8.07913626816e-06
    """
    _chkhumbnds(airf,temp,dhum,chkbnd=chkbnd)
    f_d = air_f(0,0,1,airf,temp,dhum)
    f_tt = air_f(0,2,0,airf,temp,dhum)
    f_td = air_f(0,1,1,airf,temp,dhum)
    f_dd = air_f(0,0,2,airf,temp,dhum)
    denom = dhum**2 * (f_tt*(2*f_d + dhum*f_dd) - dhum*f_td**2)
    kappa = f_tt / denom
    return kappa

def kappat(airf,temp,dhum,chkbnd=False):
    """Calculate humid air isothermal compressibility.
    
    Calculate the isothermal compressibility of humid air from dry air
    mass fraction, temperature, and humid air density.
    
    :arg float airf: Dry air mass fraction in kg/kg.
    :arg float temp: Temperature in K.
    :arg float dhum: Humid air density in kg/m3.
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :returns: Compressibility in 1/Pa.
    :raises ValueError: If temp or dhum are nonpositive or if airf is
        not between 0 and 1.
    :raises RuntimeWarning: If temp or dhum are outside the recommended
        bounds and chkbnd is True.
    
    :Examples:
    
    >>> kappat(0.9,300.,1.)
    1.09975521396e-05
    """
    _chkhumbnds(airf,temp,dhum,chkbnd=chkbnd)
    f_d = air_f(0,0,1,airf,temp,dhum)
    f_dd = air_f(0,0,2,airf,temp,dhum)
    denom = dhum**2 * (2*f_d + dhum*f_dd)
    kappa = denom**(-1)
    return kappa

def lapserate(airf,temp,dhum,chkbnd=False):
    """Calculate humid air adiabatic lapse rate.
    
    Calculate the 'dry' adiabatic lapse rate of humid air from the dry
    air mass fraction, temperature, and humid air density. 'Dry' here
    means that water does not condense or evaporate as pressure changes.
    
    :arg float airf: Dry air mass fraction in kg/kg.
    :arg float temp: Temperature in K.
    :arg float dhum: Humid air density in kg/m3.
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :returns: Lapse rate in K/Pa.
    :raises ValueError: If temp or dhum are nonpositive or if airf is
        not between 0 and 1.
    :raises RuntimeWarning: If temp or dhum are outside the recommended
        bounds and chkbnd is True.
    
    :Examples:
    
    >>> lapserate(0.9,300.,1.)
    8.50371537341e-04
    """
    _chkhumbnds(airf,temp,dhum,chkbnd=chkbnd)
    f_d = air_f(0,0,1,airf,temp,dhum)
    f_tt = air_f(0,2,0,airf,temp,dhum)
    f_td = air_f(0,1,1,airf,temp,dhum)
    f_dd = air_f(0,0,2,airf,temp,dhum)
    denom = dhum * (f_tt*(2*f_d + dhum*f_dd) - dhum*f_td**2)
    gamma = -f_td / denom
    return gamma

def pressure(airf,temp,dhum,chkbnd=False):
    """Calculate humid air pressure.
    
    Calculate the pressure in humid air from dry air mass fraction,
    temperature, and humid air density.
    
    :arg float airf: Dry air mass fraction in kg/kg.
    :arg float temp: Temperature in K.
    :arg float dhum: Humid air density in kg/m3.
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :returns: Pressure in Pa.
    :raises ValueError: If temp or dhum are nonpositive or if airf is
        not between 0 and 1.
    :raises RuntimeWarning: If temp or dhum are outside the recommended
        bounds and chkbnd is True.
    
    :Examples:
    
    >>> pressure(0.9,300.,1.)
    91175.3848662
    """
    _chkhumbnds(airf,temp,dhum,chkbnd=chkbnd)
    f_d = air_f(0,0,1,airf,temp,dhum,chkbnd=chkbnd)
    p = dhum**2 * f_d
    return p

def soundspeed(airf,temp,dhum,chkbnd=False):
    """Calculate humid air sound speed.
    
    Calculate the speed of sound in humid air from dry air mass
    fraction, temperature, and humid air density.
    
    :arg float airf: Dry air mass fraction in kg/kg.
    :arg float temp: Temperature in K.
    :arg float dhum: Humid air density in kg/m3.
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :returns: Sound speed in m/s.
    :raises ValueError: If temp or dhum are nonpositive or if airf is
        not between 0 and 1.
    :raises RuntimeWarning: If temp or dhum are outside the recommended
        bounds and chkbnd is True.
    
    :Examples:
    
    >>> soundspeed(0.9,300.,1.)
    351.817577078
    """
    _chkhumbnds(airf,temp,dhum,chkbnd=chkbnd)
    f_d = air_f(0,0,1,airf,temp,dhum)
    f_tt = air_f(0,2,0,airf,temp,dhum)
    f_td = air_f(0,1,1,airf,temp,dhum)
    f_dd = air_f(0,0,2,airf,temp,dhum)
    c_sq = dhum**2 * (f_dd - f_td**2 / f_tt) + 2*dhum*f_d
    c = c_sq**(.5)
    return c

def vappot(airf,temp,dhum,chkbnd=False):
    """Calculate humid air vapour potential.
    
    Calculate the chemical potential of water vapour in humid air from
    dry air mass fraction, temperature, and humid air density.
    
    :arg float airf: Dry air mass fraction in kg/kg.
    :arg float temp: Temperature in K.
    :arg float dhum: Humid air density in kg/m3.
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :returns: Vapour potential in J/kg.
    :raises ValueError: If temp or dhum are nonpositive or if airf is
        not between 0 and 1.
    :raises RuntimeWarning: If temp or dhum are outside the recommended
        bounds and chkbnd is True.
    """
    _chkhumbnds(airf,temp,dhum,chkbnd=chkbnd)
    f = air_f(0,0,0,airf,temp,dhum)
    f_a = air_f(1,0,0,airf,temp,dhum)
    f_d = air_f(0,0,1,airf,temp,dhum)
    gvap = f + dhum*f_d - airf*f_a
    return gvap


## Equilibrium calculation auxiliary functions
def eq_entropy(drva,drvt,drvd,airf,temp,dhum,chkbnd=False):
    """Calculate humid air entropy with derivatives.
    
    Calculate the specific entropy of humid air or its derivatives with
    respect to dry air mass fraction, temperature, and humid air
    density.
    
    :arg int drva: Number of dry fraction derivatives.
    :arg int drvt: Number of temperature derivatives.
    :arg int drvd: Number of density derivatives.
    :arg float airf: Dry air mass fraction in kg/kg.
    :arg float temp: Temperature in K.
    :arg float dhum: Humid air density in kg/m3.
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :returns: Entropy in units of
        (J/kg/K) / (kg/kg)^drva / K^drvt / (kg/m3)^drvd.
    :raises ValueError: If temp or dhum are nonpositive or if airf is
        not between 0 and 1.
    :raises RuntimeWarning: If temp or dhum are outside the recommended
        bounds and chkbnd is True.
    :raises ValueError: If any of (drva,drvt,drvd) are <0 or if
        drva+drvt+drvd>1.
    """
    _chkhumbnds(airf,temp,dhum,chkbnd=chkbnd)
    entr = -air_f(drva,drvt+1,drvd,airf,temp,dhum)
    return entr

def eq_pressure(drva,drvt,drvd,airf,temp,dhum,chkbnd=False):
    """Calculate humid air pressure with derivatives.
    
    Calculate the pressure of humid air or its derivatives with respect
    to dry air mass fraction, temperature, and humid air density.
    
    :arg int drva: Number of dry fraction derivatives.
    :arg int drvt: Number of temperature derivatives.
    :arg int drvd: Number of density derivatives.
    :arg float airf: Dry air mass fraction in kg/kg.
    :arg float temp: Temperature in K.
    :arg float dhum: Humid air density in kg/m3.
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :returns: Pressure in units of
        Pa / (kg/kg)^drva / K^drvt / (kg/m3)^drvd.
    :raises ValueError: If temp or dhum are nonpositive or if airf is
        not between 0 and 1.
    :raises RuntimeWarning: If temp or dhum are outside the recommended
        bounds and chkbnd is True.
    :raises ValueError: If any of (drva,drvt,drvd) are <0 or if
        drva+drvt+drvd>1.
    """
    _chkhumbnds(airf,temp,dhum,chkbnd=chkbnd)
    f_d = air_f(0,0,1,airf,temp,dhum)
    
    # Run through derivative cases
    if (drva,drvt,drvd) == (0,0,0):
        pres = dhum**2 * f_d
    elif (drva,drvt,drvd) == (1,0,0):
        f_ad = air_f(1,0,1,airf,temp,dhum)
        pres = dhum**2 * f_ad
    elif (drva,drvt,drvd) == (0,1,0):
        f_td = air_f(0,1,1,airf,temp,dhum)
        pres = dhum**2 * f_td
    elif (drva,drvt,drvd) == (0,0,1):
        f_dd = air_f(0,0,2,airf,temp,dhum)
        pres = 2*dhum*f_d + dhum**2 * f_dd
    else:
        errmsg = 'Derivatives {0} not recognized'.format((drva,drvt,drvd))
        raise ValueError(errmsg)
    return pres

def eq_vappot(drva,drvt,drvd,airf,temp,dhum,chkbnd=False):
    """Calculate vapour potential with derivatives.
    
    Calculate the chemical potential of water vapour in humid air or its
    derivatives with respect to dry air mass fraction, temperature, and
    humid air density. Only first-order derivatives are available.
    
    :arg int drva: Number of dry fraction derivatives.
    :arg int drvt: Number of temperature derivatives.
    :arg int drvd: Number of density derivatives.
    :arg float airf: Dry air mass fraction in kg/kg.
    :arg float temp: Temperature in K.
    :arg float dhum: Humid air density in kg/m3.
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :returns: Vapour potential in units of
        (J/kg) / (kg/kg)^drva / K^drvt / (kg/m3)^drvd.
    :raises ValueError: If temp or dhum are nonpositive or if airf is
        not between 0 and 1.
    :raises RuntimeWarning: If temp or dhum are outside the recommended
        bounds and chkbnd is True.
    :raises ValueError: If any of (drva,drvt,drvd) are <0 or if
        drva+drvt+drvd>1.
    """
    _chkhumbnds(airf,temp,dhum,chkbnd=chkbnd)
    f = air_f(0,0,0,airf,temp,dhum)
    f_a = air_f(1,0,0,airf,temp,dhum)
    f_d = air_f(0,0,1,airf,temp,dhum)
    
    # Run through derivative cases
    if (drva,drvt,drvd) == (0,0,0):
        gvap = f + dhum*f_d - airf*f_a
    elif (drva,drvt,drvd) == (1,0,0):
        f_aa = air_f(2,0,0,airf,temp,dhum)
        f_ad = air_f(1,0,1,airf,temp,dhum)
        gvap = dhum*f_ad - airf*f_aa
    elif (drva,drvt,drvd) == (0,1,0):
        f_t = air_f(0,1,0,airf,temp,dhum)
        f_at = air_f(1,1,0,airf,temp,dhum)
        f_td = air_f(0,1,1,airf,temp,dhum)
        gvap = f_t + dhum*f_td - airf*f_at
    elif (drva,drvt,drvd) == (0,0,1):
        f_ad = air_f(1,0,1,airf,temp,dhum)
        f_dd = air_f(0,0,2,airf,temp,dhum)
        gvap = 2*f_d + dhum*f_dd - airf*f_ad
    else:
        errmsg = 'Derivatives {0} not recognized'.format((drva,drvt,drvd))
        raise ValueError(errmsg)
    return gvap


## Functions to check results
def chkiapws10table(number,printresult=True,chktol=_CHKTOL):
    """Check accuracy against IAPWS 2010.
    
    Evaluate the functions in this module and compare to reference
    values in IAPWS 2010. This function is simply a wrapper for the
    three chkiapws10table* functions.
    
    :arg int number: Table number to check; 13, 14, or 15.
    :arg bool printresult: If True (default) and any results are outside
        of the given tolerance, then the function name, reference value,
        result value, and relative error are printed.
    :arg float chktol: Tolerance to use when choosing to print results
        (default _CHKTOL).
    :returns: :class:`~tester.Tester` instances containing the
        functions, arguments, reference values, results, and relative
        errors from the tests. See the individual functions to find what
        each instance corresponds to.
    :raises ValueError: If `number` is not 13, 14, or 15.
    """
    if number == 13:
        fun = chkiapws10table13
    elif number == 14:
        fun = chkiapws10table14
    elif number == 15:
        fun = chkiapws10table15
    else:
        errmsg = 'Table number {0} not available'.format(number)
        raise ValueError(errmsg)
    testers = fun(printresult=printresult,chktol=chktol)
    return testers

def chkiapws10table13(printresult=True,chktol=_CHKTOL):
    """Check accuracy against IAPWS 2010 table 13.
    
    Evaluate the functions in this module and compare to reference
    values from IAPWS 2010, table 13. The functions include the free
    energy air_f and related thermodynamic quantities (e.g. heat
    capacity, entropy).
    
    :arg bool printresult: If True (default) and any results are outside
        of the given tolerance, then the function name, reference value,
        result value, and relative error are printed.
    :arg float chktol: Tolerance to use when choosing to print results
        (default _CHKTOL).
    :returns: :class:`~tester.Tester` instances containing the
        functions, arguments, reference values, results, and relative
        errors from the tests. The first instance deals with derivatives
        of air_f; the second instance with other thermodynamic
        functions.
    :rtype: (Tester,Tester)
    """
    from tester import Tester
    ders = ((0,0,0),(1,0,0),(0,1,0),(0,0,1),(2,0,0),(1,1,0),(1,0,1),
        (0,2,0),(0,1,1),(0,0,2))
    fargs0 = (0.892247719,200.,0.163479657e-4)
    fargs1 = (0.977605798,300.,0.114614216e+1)
    fargs2 = (0.825565291,400.,0.793354063e+1)
    propfargs = [fargs0,fargs1,fargs2]
    
    # Tester instance for derivative of air_f
    derfuns = air_f
    derfargs = [(der+fargs) for fargs in propfargs for der in ders]
    # Derivatives change before arguments do here
    derrefs = [-0.682093392e+6,-0.572680404e+6,-0.405317966e+4,0.374173101e+10,
        0.920967684e+6,0.915653743e+4,-0.213442099e+10,-0.394011921e+1,
        0.187087034e+8,-0.228880603e+15,-0.927718178e+5,-0.263453864e+3,
        -0.296711481e+3,0.761242496e+5,0.624886233e+7,0.822733446e+4,
        -0.450004399e+5,-0.244742952e+1,0.254456302e+3,-0.664465525e+5,
        0.240345570e+5,0.311096733e+6,-0.106891931e+4,0.158878781e+5,
        0.113786423e+7,0.702631471e+4,-0.727972651e+4,-0.222449294e+1,
        0.414350772e+2,-0.201886184e+4]
    derfnames = 'air_f'
    derargfmt = '({0:1g},{1:1g},{2:1g},{3:8.6f},{4:3g},{5:10.4e})'
    header = 'Humid air Helmholtz energy'
    dertest = Tester(derfuns,derfargs,derrefs,derfnames,derargfmt,
        header=header)
    
    # Tester instance for other air properties
    propfuns = [pressure,enthalpy,gibbsenergy,entropy,cp,soundspeed,vappot]
    proprefs = [
        [0.999999998,0.1e+6,0.1e+7],
        [0.189712231e+6,0.834908383e+5,0.577649408e+6],
        [-0.620923701e+6,-0.552260595e+4,0.150081684e+6],
        [0.405317966e+4,0.296711481e+3,0.106891931e+4],
        [0.109387397e+4,0.102681324e+4,0.123552454e+4],
        [0.291394959e+3,0.349234196e+3,0.416656820e+3],
        [-0.109950917e+6,-0.526505193e+4,-0.106748981e+6]
    ]
    propfnames = ['pressure','enthalpy','gibbsenergy','entropy','cp',
        'soundspeed','vappot']
    propargfmt = '({0:8.6f},{1:3g},{2:10.4e})'
    header = 'Humid air thermodynamic properties'
    proptest = Tester(propfuns,propfargs,proprefs,propfnames,propargfmt,
        header=header)
    
    # Run Tester instances and print results
    dertest.run()
    proptest.run()
    if printresult:
        dertest.printresults(chktol=chktol)
        proptest.printresults(chktol=chktol)
    return dertest, proptest

def chkiapws10table14(printresult=True,chktol=_CHKTOL):
    """Check accuracy against IAPWS 2010 table 14.
    
    Evaluate the functions in this module and compare to reference
    values from IAPWS 2010, table 14. These values separate the
    contributions from dry air and water vapour.
    
    :arg bool printresult: If True (default) and any results are outside
        of the given tolerance, then the function name, reference value,
        result value, and relative error are printed.
    :arg float chktol: Tolerance to use when choosing to print results
        (default _CHKTOL).
    :returns: :class:`~tester.Tester` instances containing the
        functions, arguments, reference values, results, and relative
        errors from the tests. The first instance deals with dry air and
        the second instance with water vapour.
    :rtype: (Tester,Tester)
    """
    from tester import Tester
    ders = ((0,0),(1,0),(0,1),(2,0),(1,1),(0,2))
    temps = (200.,300.,400.)
    dhums = (0.163479657e-4,0.114614216e+1,0.793354063e+1)
    airfs = (0.892247719,0.977605798,0.825565291)
    argfmt = '({0:1g},{1:1g},{2:3g},{3:10.4e})'
    # Derivatives change before arguments do here
    dryfargs = [(der+(temp,dhum*airf))
        for (airf,temp,dhum) in zip(airfs,temps,dhums) for der in ders]
    vapfargs = [(der+(temp,dhum*(1-airf)))
        for (airf,temp,dhum) in zip(airfs,temps,dhums) for der in ders]
    
    # Dry air tester
    dryfuns = _dry_f
    dryrefs = [-0.740041144e+6,-0.304774177e+4,0.393583654e+10,-0.357677878e+1,
        0.196791837e+8,-0.269828549e+15,-0.916103453e+5,-0.108476220e+3,
        0.768326795e+5,-0.239319940e+1,0.256683306e+3,-0.685917373e+5,
        0.895561286e+5,0.193271394e+3,0.175560114e+5,-0.181809877e+1,
        0.442769673e+2,-0.267635928e+4]
    dryfnames = 'dry_f'
    header = 'Dry air Helmholtz energy'
    drytest = Tester(dryfuns,dryfargs,dryrefs,dryfnames,argfmt,header=header)
    
    # Water vapour tester
    vapfuns = _flu_f
    vaprefs = [-0.202254351e+6,-0.123787544e+5,0.523995674e+11,-0.694877601e+1,
        0.262001885e+9,-0.297466671e+17,-0.143157426e+6,-0.851598213e+4,
        0.538480619e+7,-0.480817011e+1,0.181489502e+5,-0.210184992e+9,
        -0.285137534e+6,-0.705288048e+4,0.129645039e+6,-0.411710659e+1,
        0.361784086e+3,-0.965539462e+5]
    vapfnames = 'flu_f'
    header = 'Water vapour Helmholtz energy'
    vaptest = Tester(vapfuns,vapfargs,vaprefs,vapfnames,argfmt,header=header)
    
    # Run Tester instances and print results
    drytest.run()
    vaptest.run()
    if printresult:
        drytest.printresults(chktol=chktol)
        vaptest.printresults(chktol=chktol)
    return drytest, vaptest

def chkiapws10table15(printresult=True,chktol=_CHKTOL):
    """Check accuracy against IAPWS 2010 table 15.
    
    Evaluate the functions in this module and compare to reference
    values from IAPWS 2010, table 15. These values are for the non-ideal
    mixing of dry air and water vapour, i.e. the virial coefficients.
    
    :arg bool printresult: If True (default) and any results are outside
        of the given tolerance, then the function name, reference value,
        result value, and relative error are printed.
    :arg float chktol: Tolerance to use when choosing to print results
        (default _CHKTOL).
    :returns: :class:`~tester.Tester` instances containing the
        functions, arguments, reference values, results, and relative
        errors from the tests. The first instance involves each virial
        coefficient separately; the second involves the entire excess
        energy of mixing.
    :rtype: (Tester,Tester)
    """
    from tester import Tester
    
    # Test virial coefficients
    virfuns = [_air_baw, _air_caaw, _air_caww]
    virders = (0,1,2)
    temps = (200.,300.,400.)
    virfargs = [(der,temp) for temp in temps for der in virders]
    virrefs = [
        [-0.784874277752e-4,0.848076624222e-6,-0.122622146106e-7,
            -0.295672747428e-4,0.280097360438e-6,-0.242599241306e-8,
            -0.100804610474e-4,0.135021228495e-6,-0.839901728946e-9],
        [0.105493575000e-8,-0.152535000000e-11,-0.113436375000e-12,
            0.801977740741e-9,-0.196103456790e-11,0.170055637860e-13,
            0.672018171875e-9,-0.812416406250e-12,0.683147460938e-14],
        [-0.349872634207e-5,0.188025052349e-6,-0.124996855887e-7,
            -0.115552783680e-6,0.261363277754e-8,-0.751334581804e-10,
            -0.200806020909e-7,0.274535402840e-9,-0.491763909891e-11]
    ]
    virfnames = ['air_baw','air_caaw','air_caww']
    virargfmt = '({0:1g},{1:3g})'
    header = 'Virial coefficients'
    virtest = Tester(virfuns,virfargs,virrefs,virfnames,virargfmt,
        header=header)
    
    # Test energy of mixing
    mixfuns = _air_f_mix
    mixders = ((0,0,0),(1,0,0),(0,1,0),(0,0,1),(2,0,0),(1,1,0),(1,0,1),
        (0,2,0),(0,1,1),(0,0,2))
    mixargs = [
        (0.892224944,200.,0.163445112e-04),
        (0.977600624,300.,0.114587678e+1),
        (0.825531379,400.,0.793198757e+1)
    ]
    mixfargs = [(der+arg) for arg in mixargs for der in mixders]
    mixrefs = [-0.786211837111e-3,0.641377589024e-2,0.456427011454e-5,
        -0.481026562160e+2,0.163518396765e-1,-0.372355250910e-4,
        0.392414346187e+3,-0.378866038500e-7,0.279261664337,
        -0.192118914304e+2,-0.711673565972e+1,0.311768500757e+3,
        0.441245367098e-1,-0.623171267739e+1,0.534139178420e+3,
        -0.195026097643e+1,0.274152648677e+3,-0.148782305257e-3,
        0.390100461800e-1,-0.366162709444e-1,-0.161985033872e+3,
        0.830802876130e+3,0.178961265299e+1,-0.223365431713e+2,
        0.135815609516e+4,-0.916586082354e+1,0.125823777783e+3,
        -0.536718535916e-2,0.249618216264,-0.482803925450]
    mixfnames = 'air_f_mix'
    mixargfmt = '({0:1g},{1:1g},{2:1g},{3:8.6f},{4:3g},{5:10.4e})'
    header = 'Free energy of mixing'
    mixtest = Tester(mixfuns,mixfargs,mixrefs,mixfnames,mixargfmt,
        header=header)
    
    # Run Tester instances and print results
    virtest.run()
    mixtest.run()
    if printresult:
        virtest.printresults(chktol=chktol)
        mixtest.printresults(chktol=chktol)
    return virtest, mixtest


## Main function: Check tables
if __name__ == '__main__':
    dertest, proptest = chkiapws10table13()
    drytest, vaptest = chkiapws10table14()
    virtest, mixtest = chkiapws10table15()

