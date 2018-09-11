"""Seawater salt Gibbs free energy and thermodynamic properties.

This module provides the contribution of salt to the Gibbs free energy of seawater and its derivatives with respect to salinity, temperature, and pressure. It also provides thermodynamic properties that can be derived directly from this energy, such as the dilution and osmotic coefficients.

:Examples:

>>> sal_g(0,0,0,0.035,300.,1e6)
127.033640309
>>> sal_g(0,0,1,0.035,300.,1e6)
-2.55600080319e-05
>>> sal_g(0,2,0,0.035,300.,1e6)
0.597842170749
>>> activityw(0.035,300.,1e6)
0.981388410188
>>> saltenthalpy(0.035,300.,1e6)
-156107.959196
>>> mixvolume(0.01,0.035,0.6,300.,1e6)
-5.94174956892e-08

:Functions:

* sal_g: Gibbs free energy of salt in seawater.
* actcoeff: Activity coefficient of salt in seawater.
* actpotential: Activity potential of salt in seawater.
* dilution: Dilution coefficient of salt in seawater.
* osmcoeff: Osmotic coefficient of salt in seawater.
* activityw: Activity of water in seawater.
* liqpot: Chemical potential of water in seawater.
* salpot: Chemical potential of salt in seawater.
* chemcoeff: Chemical coefficient of salt in seawater.
* saltenthalpy: Specific enthalpy of salt in seawater.
* saltentropy: Specific entropy of salt in seawater.
* saltvolume: Specific volume of salt in seawater.
* mixenthalpy: Specific enthalpy of mixing of seawater parcels.
* mixentropy: Specific entropy of mixing of seawater parcels.
* mixvolume: Specific volume of mixing of seawater parcels.
* eq_liqpot: Chemical potential of water in seawater with derivatives.
* eq_entropy: Specific entropy of seawater with derivatives.
* eq_enthalpy: Specific enthalpy of seawater with derivatives.

"""

__all__ = ['sal_g',
    'actcoeff','activityw','actpotential','chemcoeff','dilution','liqpot',
    'osmcoeff','salpot','saltenthalpy','saltentropy','saltvolume',
    'mixenthalpy','mixentropy','mixvolume',
    'eq_liqpot','eq_enthalpy','eq_entropy']

import numpy
import constants0
import sal1

_chksalbnds = constants0.chksalbnds
_MWAT = constants0.MWAT
_MSAL = constants0.MSAL
_NSALTERMS = sal1.NSALTERMS
_sal_g_term = sal1.sal_g_term


## Gibbs energy function
def sal_g(drvs,drvt,drvp,salt,temp,pres,chkbnd=False,useext=False):
    """Calculate salt Gibbs free energy.
    
    Calculate the contribution of salt to the specific Gibbs free energy
    of seawater, or its derivatives with respect to salinity,
    temperature, and pressure. Derivatives of any order can be
    specified, but only derivatives up to second order have been error-
    checked.
    
    :arg int drvs: Number of salinity derivatives.
    :arg int drvt: Number of temperature derivatives.
    :arg int drvp: Number of pressure derivatives.
    :arg float salt: Absolute salinity in kg/kg.
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg bool useext: If False (default) then the polynomial is
        calculated from _GSCOEFFS; if True, from _GSCOEFFS_EXT.
    :returns: Gibbs free energy in units of
        (J/kg) / (kg/kg)^drvs / K^drvt / Pa^drvp.
    :raises ValueError: If temp or pres are nonpositive, or if salt is
        not between 0 and 1.
    :raises RuntimeWarning: If salt, temp, or pres are outside the
        recommended bounds and chkbnd is True.
    :raises ValueError: If any of drvs, drvt, or drvp are negative.
    :raises ValueError: If salt is 0, unless (drvs is 0) or
        ((drvs is 1) and (drvt>1) or (drvp>0)).
    
    :Examples:
    
    >>> sal_g(0,0,0,0.035,300.,1e6)
    127.033640309
    >>> sal_g(0,0,1,0.035,300.,1e6)
    -2.55600080319e-05
    >>> sal_g(0,2,0,0.035,300.,1e6)
    0.597842170749
    """
    _chksalbnds(salt,temp,pres,chkbnd=chkbnd)
    if min(drvs,drvt,drvp) < 0:
        errmsg = 'Derivatives {0} not recognized'.format((drvs,drvt,drvp))
        raise ValueError(errmsg)
    
    # Handle the case where salinity approaches 0
    if salt == 0:
        if drvs == 0:
            g = 0.
            return g
        elif (drvs == 1) and (drvp > 0 or drvt > 1):
            g2 = _sal_g_term(2,drvt,drvp,temp,pres,useext=useext)
            g = g2
            return g
        else:
            # All other terms have a singularity
            errmsg = ('Salinity of 0 is not compatible with derivatives '
                '{0}').format((drvs,drvt,drvp))
            raise ValueError(errmsg)
    
    # Calculate the Gibbs energy using salinity expansion
    gis = [_sal_g_term(iTerm,drvt,drvp,temp,pres,useext=useext)
        for iTerm in range(1,_NSALTERMS+1)]
    g = 0.0
    if gis[0] != 0:
        if drvs == 0:
            g += gis[0] * salt * numpy.log(salt)
        elif drvs == 1:
            g += gis[0] * (1 + numpy.log(salt))
        else:
            g += gis[0] * (-1.)**(drvs) / salt**(drvs-1)
    iStart = max(2*drvs-1,1)
    for iTerm in range(iStart,_NSALTERMS,2):
        if gis[iTerm] != 0:
            term = salt**((iTerm+1)/2 - drvs)
            for l in range(drvs):
                term *= ((iTerm+1)/2 - l)
            g += gis[iTerm] * term
    rs = salt**.5
    for iTerm in range(2,_NSALTERMS,2):
        if gis[iTerm] != 0:
            term = rs**((iTerm+1) - 2*drvs)
            for l in range(drvs):
                term *= (.5*(iTerm+1) - l)
            g += gis[iTerm] * term
    return g


## Thermodynamic properties
def actcoeff(salt,temp,pres,chkbnd=False,useext=False):
    """Calculate salt activity coefficient.
    
    Calculate the mean activity coefficient ln(gamma) of salt in
    seawater from salinity, temperature, and pressure.
    
    :arg float salt: Absolute salinity in kg/kg.
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg bool useext: If False (default) then the polynomial is
        calculated from _GSCOEFFS; if True, from _GSCOEFFS_EXT.
    :returns: Activity coefficient, unitless.
    :raises ValueError: If temp or pres are nonpositive, or if salt is
        not between 0 and 1.
    :raises RuntimeWarning: If salt, temp, or pres are outside the
        recommended bounds and chkbnd is True.
    
    :Examples:
    
    >>> actcoeff(0.035,300.,1e6)
    -0.527003008913
    """
    _chksalbnds(salt,temp,pres,chkbnd=chkbnd)
    # Treat boundary cases separately
    if salt == 1:
        errmsg = 'Activity coefficient is undefined for a salinity of 1'
        raise ValueError(errmsg)
    if salt == 0:
        lng = 0.
        return lng
    
    g = sal_g(0,0,0,salt,temp,pres,useext=useext)
    g_s = sal_g(1,0,0,salt,temp,pres,useext=useext)
    g1 = _sal_g_term(1,0,0,temp,pres,useext=useext)
    g2 = _sal_g_term(2,0,0,temp,pres,useext=useext)
    lng = (g + (1-salt)*g_s - g1-g2)/g1 + numpy.log(salt**(-1)-1)
    return lng

def activityw(salt,temp,pres,chkbnd=False,useext=False):
    """Calculate water activity in seawater.
    
    Calculate the activity of water in seawater from salinity,
    temperature, and pressure.
    
    :arg float salt: Absolute salinity in kg/kg.
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg bool useext: If False (default) then the polynomial is
        calculated from _GSCOEFFS; if True, from _GSCOEFFS_EXT.
    :returns: Activity of water, unitless.
    :raises ValueError: If temp or pres are nonpositive, or if salt is
        not between 0 and 1.
    :raises RuntimeWarning: If salt, temp, or pres are outside the
        recommended bounds and chkbnd is True.
    
    :Examples:
    
    >>> activityw(0.035,300.,1e6)
    0.981388410188
    """
    _chksalbnds(salt,temp,pres,chkbnd=chkbnd)
    # Treat boundary cases separately
    if salt == 0:
        awat = 1.
        return awat
    if salt == 1:
        g1 = _sal_g_term(1,0,0,temp,pres,useext=useext)
        gis = [_sal_g_term(iTerm,0,0,temp,pres,useext=useext)
            for iTerm in range(3,_NSALTERMS+1)]
        earg = salt
        for (ind,gi) in enumerate(gis):
            iTerm = ind + 3
            earg += gi/g1 * (.5*iTerm-1) * salt**(.5*iTerm)
        earg *= -_MWAT/_MSAL
        awat = numpy.exp(earg)
        return awat
    
    phi = osmcoeff(salt,temp,pres,useext=useext)
    awat = numpy.exp(-phi * salt/(1-salt) * _MWAT/_MSAL)
    return awat

def actpotential(salt,temp,pres,chkbnd=False,useext=False):
    """Calculate salt activity potential.
    
    Calculate the activity potential of salt in seawater from salinity,
    temperature, and pressure.
    
    :arg float salt: Absolute salinity in kg/kg.
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg bool useext: If False (default) then the polynomial is
        calculated from _GSCOEFFS; if True, from _GSCOEFFS_EXT.
    :returns: Activity potential, unitless.
    :raises ValueError: If temp or pres are nonpositive, or if salt is
        not between 0 and 1.
    :raises RuntimeWarning: If salt, temp, or pres are outside the
        recommended bounds and chkbnd is True.
    
    :Examples:
    
    >>> actpotential(0.035,300.,1e6)
    -0.429940465498
    """
    _chksalbnds(salt,temp,pres,chkbnd=chkbnd)
    # Treat boundary cases separately
    if salt == 1:
        errmsg = 'Activity potential is undefined for a salinity of 1'
        raise ValueError(errmsg)
    if salt == 0:
        psi = 0.
        return psi
    
    g = sal_g(0,0,0,salt,temp,pres,useext=useext)
    g1 = _sal_g_term(1,0,0,temp,pres,useext=useext)
    g2 = _sal_g_term(2,0,0,temp,pres,useext=useext)
    psi = (g/salt - g2)/g1 + numpy.log(salt**(-1)-1)
    return psi

def chemcoeff(salt,temp,pres,chkbnd=False,useext=False):
    """Calculate salt chemical coefficient.
    
    Calculate the chemical coefficient of salt in seawater from
    salinity, temperature, and pressure.
    
    :arg float salt: Absolute salinity in kg/kg.
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg bool useext: If False (default) then the polynomial is
        calculated from _GSCOEFFS; if True, from _GSCOEFFS_EXT.
    :returns: Chemical coefficient in J/kg.
    :raises ValueError: If temp or pres are nonpositive, or if salt is
        not between 0 and 1.
    :raises RuntimeWarning: If salt, temp, or pres are outside the
        recommended bounds and chkbnd is True.
    
    :Examples:
    
    >>> chemcoeff(0.035,300.,1e6)
    2754.04566958
    """
    _chksalbnds(salt,temp,pres,chkbnd=chkbnd)
    d = dilution(salt,temp,pres,useext=useext)
    chem = salt * d
    return chem

def dilution(salt,temp,pres,chkbnd=False,useext=False):
    """Calculate salt dilution coefficient.
    
    Calculate the dilution coefficient of salt in seawater from
    salinity, temperature, and pressure.
    
    :arg float salt: Absolute salinity in kg/kg.
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg bool useext: If False (default) then the polynomial is
        calculated from _GSCOEFFS; if True, from _GSCOEFFS_EXT.
    :returns: Dilution coefficient in J/kg.
    :raises ValueError: If temp or pres are nonpositive, or if salt is
        not between 0 and 1.
    :raises RuntimeWarning: If salt, temp, or pres are outside the
        recommended bounds and chkbnd is True.
    
    :Examples:
    
    >>> dilution(0.035,300.,1e6)
    78687.0191309
    """
    _chksalbnds(salt,temp,pres,chkbnd=chkbnd)
    # Treat boundary case separately
    if salt == 0:
        dil = _sal_g_term(1,0,0,temp,pres,useext=useext)
        return dil
    g_ss = sal_g(2,0,0,salt,temp,pres,useext=useext)
    dil = salt*g_ss
    return dil

def liqpot(salt,temp,pres,chkbnd=False,useext=False):
    """Calculate water chemical potential.
    
    Calculate the contribution of salinity to the chemical potential of
    water in seawater from salinity, temperature, and pressure.
    
    :arg float salt: Absolute salinity in kg/kg.
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg bool useext: If False (default) then the polynomial is
        calculated from _GSCOEFFS; if True, from _GSCOEFFS_EXT.
    :returns: Chemical potential in J/kg.
    :raises ValueError: If temp or pres are nonpositive, or if salt is
        not between 0 and 1.
    :raises RuntimeWarning: If salt, temp, or pres are outside the
        recommended bounds and chkbnd is True.
    
    :Examples:
    
    >>> liqpot(0.035,300.,1e6)
    -2601.18871107
    """
    _chksalbnds(salt,temp,pres,chkbnd=chkbnd)
    g = sal_g(0,0,0,salt,temp,pres,useext=useext)
    g_s = sal_g(1,0,0,salt,temp,pres,useext=useext)
    gliq = g - salt*g_s
    return gliq

def osmcoeff(salt,temp,pres,chkbnd=False,useext=False):
    """Calculate salt osmotic coefficient.
    
    Calculate the osmotic coefficient phi of salt in seawater from
    salinity, temperature, and pressure.
    
    :arg float salt: Absolute salinity in kg/kg.
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg bool useext: If False (default) then the polynomial is
        calculated from _GSCOEFFS; if True, from _GSCOEFFS_EXT.
    :returns: Osmotic coefficient, unitless.
    :raises ValueError: If temp or pres are nonpositive, or if salt is
        not between 0 and 1.
    :raises RuntimeWarning: If salt, temp, or pres are outside the
        recommended bounds and chkbnd is True.
    
    :Examples:
    
    >>> osmcoeff(0.035,300.,1e6)
    0.902937456585
    """
    _chksalbnds(salt,temp,pres,chkbnd=chkbnd)
    # Treat boundary case separately
    if salt == 0:
        phi = 1.
        return phi
    
    g = sal_g(0,0,0,salt,temp,pres,useext=useext)
    g_s = sal_g(1,0,0,salt,temp,pres,useext=useext)
    g1 = _sal_g_term(1,0,0,temp,pres,useext=useext)
    
    # Calculate phi
    phi = (g_s - g/salt) * (1-salt) / g1
    return phi

def salpot(salt,temp,pres,chkbnd=False,useext=False):
    """Calculate salt chemical potential.
    
    Calculate the chemical potential of salt in seawater from salinity,
    temperature, and pressure.
    
    :arg float salt: Absolute salinity in kg/kg.
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg bool useext: If False (default) then the polynomial is
        calculated from _GSCOEFFS; if True, from _GSCOEFFS_EXT.
    :returns: Chemical potential in J/kg.
    :raises ValueError: If temp or pres are nonpositive, or if salt is
        not between 0 and 1.
    :raises RuntimeWarning: If salt, temp, or pres are outside the
        recommended bounds and chkbnd is True.
    
    :Examples:
    
    >>> salpot(0.035,300.,1e6)
    77949.2100395
    """
    _chksalbnds(salt,temp,pres,chkbnd=chkbnd)
    gsal = sal_g(1,0,0,salt,temp,pres,useext=useext)
    return gsal

def saltenthalpy(salt,temp,pres,chkbnd=False,useext=False):
    """Calculate seawater enthalpy salt term.
    
    Calculate the contribution of salt to the specific enthalpy of
    seawater from salinity, temperature, and pressure.
    
    :arg float salt: Absolute salinity in kg/kg.
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg bool useext: If False (default) then the polynomial is
        calculated from _GSCOEFFS; if True, from _GSCOEFFS_EXT.
    :returns: Enthalpy in J/kg/(kg/kg).
    :raises ValueError: If temp or pres are nonpositive, or if salt is
        not between 0 and 1.
    :raises RuntimeWarning: If salt, temp, or pres are outside the
        recommended bounds and chkbnd is True.
    
    :Examples:
    
    >>> saltenthalpy(0.035,300.,1e6)
    -156107.959196
    """
    _chksalbnds(salt,temp,pres,chkbnd=chkbnd)
    # Treat boundary case separately
    if salt == 0:
        errmsg = 'Salt enthalpy is undefined for a salinity of 0'
        raise ValueError(errmsg)
    
    g = sal_g(0,0,0,salt,temp,pres,useext=useext)
    g_t = sal_g(0,1,0,salt,temp,pres,useext=useext)
    hsal = (g - temp*g_t) / salt
    return hsal

def saltentropy(salt,temp,pres,chkbnd=False,useext=False):
    """Calculate seawater entropy salt term.
    
    Calculate the contribution of salt to the specific entropy of
    seawater from salinity, temperature, and pressure.
    
    :arg float salt: Absolute salinity in kg/kg.
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg bool useext: If False (default) then the polynomial is
        calculated from _GSCOEFFS; if True, from _GSCOEFFS_EXT.
    :returns: Entropy in J/kg/K/(kg/kg).
    :raises ValueError: If temp or pres are nonpositive, or if salt is
        not between 0 and 1.
    :raises RuntimeWarning: If salt, temp, or pres are outside the
        recommended bounds and chkbnd is True.
    
    :Examples:
    
    >>> saltentropy(0.035,300.,1e6)
    -532.458305922
    """
    _chksalbnds(salt,temp,pres,chkbnd=chkbnd)
    # Treat boundary case separately
    if salt == 0:
        errmsg = 'Salt entropy is undefined for a salinity of 0'
        raise ValueError(errmsg)
    
    g_t = sal_g(0,1,0,salt,temp,pres,useext=useext)
    ssal = -g_t / salt
    return ssal

def saltvolume(salt,temp,pres,chkbnd=False,useext=False):
    """Calculate seawater volume salt term.
    
    Calculate the contribution of salt to the specific volume of
    seawater from salinity, temperature, and pressure.
    
    :arg float salt: Absolute salinity in kg/kg.
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg bool useext: If False (default) then the polynomial is
        calculated from _GSCOEFFS; if True, from _GSCOEFFS_EXT.
    :returns: Specific volume in m3/kg/(kg/kg).
    :raises ValueError: If temp or pres are nonpositive, or if salt is
        not between 0 and 1.
    :raises RuntimeWarning: If salt, temp, or pres are outside the
        recommended bounds and chkbnd is True.
    
    :Examples:
    
    >>> saltvolume(0.035,300.,1e6)
    -7.30285943768e-04
    """
    _chksalbnds(salt,temp,pres,chkbnd=chkbnd)
    # Treat boundary case separately
    if salt == 0:
        g1_p = _sal_g_term(1,0,1,temp,pres,useext=useext)
        v = g1_p
        return v
    
    g_p = sal_g(0,0,1,salt,temp,pres,useext=useext)
    vsal = g_p / salt
    return vsal


## Functions for mixing of parcels
def mixenthalpy(salt1,salt2,wgt1,temp,pres,chkbnd=False,useext=False):
    """Calculate salt mixing enthalpy.
    
    Calculate the excess enthalpy of mixing. This is the change in
    specific enthalpy that occurs when two seawater parcels with
    salinities (salt1,salt2) and mass ratio wgt1 with common temperature
    and pressure are mixed.
    
    :arg float salt1: Absolute salinity of parcel 1 in kg/kg.
    :arg float salt2: Absolute salinity of parcel 2 in kg/kg.
    :arg float wgt1: Ratio of the mass of parcel 1 to the total mass.
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg bool useext: If False (default) then the polynomial is
        calculated from _GSCOEFFS; if True, from _GSCOEFFS_EXT.
    :returns: Enthalpy of mixing in J/kg.
    :raises ValueError: If temp or pres are nonpositive, or if salt1 or
        salt2 are not between 0 and 1.
    :raises RuntimeWarning: If salt1, salt2, temp, or pres are outside
        the recommended bounds and chkbnd is True.
    :raises ValueError: If wgt1 is not between 0 and 1.
    
    :Examples:
    
    >>> mixenthalpy(0.01,0.035,0.6,300.,1e6)
    16.1539911284
    """
    _chksalbnds(salt1,temp,pres,chkbnd=chkbnd)
    _chksalbnds(salt2,temp,pres,chkbnd=chkbnd)
    if wgt1 < 0 or wgt1 > 1:
        errmsg = 'Mass fraction {0} is not between 0 and 1'.format(wgt1)
        raise ValueError(errmsg)
    
    h1 = salt1 * saltenthalpy(salt1,temp,pres,useext=useext)
    h2 = salt2 * saltenthalpy(salt2,temp,pres,useext=useext)
    salt0 = wgt1*salt1 + (1-wgt1)*salt2
    h0 = salt0 * saltenthalpy(salt0,temp,pres,useext=useext)
    hmix = h0 - wgt1*h1 - (1-wgt1)*h2
    return hmix

def mixentropy(salt1,salt2,wgt1,temp,pres,chkbnd=False,useext=False):
    """Calculate salt mixing entropy.
    
    Calculate the excess entropy of mixing. This is the change in
    specific entropy that occurs when two seawater parcels with
    salinities (salt1,salt2) and mass ratio wgt1 with common temperature
    and pressure are mixed.
    
    :arg float salt1: Absolute salinity of parcel 1 in kg/kg.
    :arg float salt2: Absolute salinity of parcel 2 in kg/kg.
    :arg float wgt1: Ratio of the mass of parcel 1 to the total mass.
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg bool useext: If False (default) then the polynomial is
        calculated from _GSCOEFFS; if True, from _GSCOEFFS_EXT.
    :returns: Entropy of mixing in J/kg/K.
    :raises ValueError: If temp or pres are nonpositive, or if salt1 or
        salt2 are not between 0 and 1.
    :raises RuntimeWarning: If salt1, salt2, temp, or pres are outside
        the recommended bounds and chkbnd is True.
    :raises ValueError: If wgt1 is not between 0 and 1.
    
    :Examples:
    
    >>> mixentropy(0.01,0.035,0.6,300.,1e6)
    0.966829422617
    """
    _chksalbnds(salt1,temp,pres,chkbnd=chkbnd)
    _chksalbnds(salt2,temp,pres,chkbnd=chkbnd)
    if wgt1 < 0 or wgt1 > 1:
        errmsg = 'Mass fraction {0} is not between 0 and 1'.format(wgt1)
        raise ValueError(errmsg)
    
    s1 = salt1 * saltentropy(salt1,temp,pres,useext=useext)
    s2 = salt2 * saltentropy(salt2,temp,pres,useext=useext)
    salt0 = wgt1*salt1 + (1-wgt1)*salt2
    s0 = salt0 * saltentropy(salt0,temp,pres,useext=useext)
    smix = s0 - wgt1*s1 - (1-wgt1)*s2
    return smix

def mixvolume(salt1,salt2,wgt1,temp,pres,chkbnd=False,useext=False):
    """Calculate salt mixing volume.
    
    Calculate the excess volume of mixing. This is the change in
    specific volume that occurs when two seawater parcels with
    salinities (salt1,salt2) and mass ratio wgt1 with common temperature
    and pressure are mixed.
    
    :arg float salt1: Absolute salinity of parcel 1 in kg/kg.
    :arg float salt2: Absolute salinity of parcel 2 in kg/kg.
    :arg float wgt1: Ratio of the mass of parcel 1 to the total mass.
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg bool useext: If False (default) then the polynomial is
        calculated from _GSCOEFFS; if True, from _GSCOEFFS_EXT.
    :returns: Volume of mixing in m3/kg.
    :raises ValueError: If temp or pres are nonpositive, or if salt1 or
        salt2 are not between 0 and 1.
    :raises RuntimeWarning: If salt1, salt2, temp, or pres are outside
        the recommended bounds and chkbnd is True.
    :raises ValueError: If wgt1 is not between 0 and 1.
    
    :Examples:
    
    >>> mixvolume(0.01,0.035,0.6,300.,1e6)
    -5.94174956892e-08
    """
    _chksalbnds(salt1,temp,pres,chkbnd=chkbnd)
    _chksalbnds(salt2,temp,pres,chkbnd=chkbnd)
    if wgt1 < 0 or wgt1 > 1:
        errmsg = 'Mass fraction {0} is not between 0 and 1'.format(wgt1)
        raise ValueError(errmsg)
    
    v1 = salt1 * saltvolume(salt1,temp,pres,useext=useext)
    v2 = salt2 * saltvolume(salt2,temp,pres,useext=useext)
    salt0 = wgt1*salt1 + (1-wgt1)*salt2
    v0 = salt0 * saltvolume(salt0,temp,pres,useext=useext)
    vmix = v0 - wgt1*v1 - (1-wgt1)*v2
    return vmix


## Equilibrium calculation auxiliary functions
def eq_liqpot(drvs,drvt,drvp,salt,temp,pres,chkbnd=False,useext=False):
    """Calculate water chemical potential with derivatives.
    
    Calculate the contribution of salt to the chemical potential of
    water in seawater or its derivatives with respect to salinity,
    temperature, and pressure. Only derivatives up to first order have
    been error checked.
    
    :arg int drvs: Number of salinity derivatives.
    :arg int drvt: Number of temperature derivatives.
    :arg int drvp: Number of pressure derivatives.
    :arg float salt: Absolute salinity in kg/kg.
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg bool useext: If False (default) then the polynomial is
        calculated from _GSCOEFFS; if True, from _GSCOEFFS_EXT.
    :returns: Chemical potential of water in units of
        (J/kg) / (kg/kg)^drvs / K^drvt / Pa^drvp.
    :raises ValueError: If temp or pres are nonpositive, or if salt is
        not between 0 and 1.
    :raises RuntimeWarning: If salt, temp, or pres are outside the
        recommended bounds and chkbnd is True.
    :raises ValueError: If any of (drvs,drvt,drvp) are negative, or if
        drvs>1.
    """
    _chksalbnds(salt,temp,pres,chkbnd=chkbnd)
    if drvs == 0:
        gs = sal_g(0,drvt,drvp,salt,temp,pres,useext=useext)
        gs_s = sal_g(1,drvt,drvp,salt,temp,pres,useext=useext)
        gliq = gs - salt*gs_s
    elif drvs == 1:
        gs_ss = sal_g(2,drvt,drvp,salt,temp,pres,useext=useext)
        gliq = -salt*gs_ss
    else:
        errmsg = 'Derivatives {0} not recognized'.format((drvs,drvt,drvp))
        raise ValueError(errmsg)
    return gliq

def eq_enthalpy(drvs,drvt,drvp,salt,temp,pres,chkbnd=False,useext=False):
    """Calculate salt enthalpy with derivatives.
    
    Calculate the contribution of salt to the specific enthalpy of
    seawater or its derivatives with respect to salinity, temperature,
    and pressure. Only derivatives up to first order have been error
    checked.
    
    :arg int drvs: Number of salinity derivatives.
    :arg int drvt: Number of temperature derivatives.
    :arg int drvp: Number of pressure derivatives.
    :arg float salt: Absolute salinity in kg/kg.
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg bool useext: If False (default) then the polynomial is
        calculated from _GSCOEFFS; if True, from _GSCOEFFS_EXT.
    :returns: Enthalpy in units of
        (J/kg) / (kg/kg)^drvs / K^drvt / Pa^drvp.
    :raises ValueError: If temp or pres are nonpositive, or if salt is
        not between 0 and 1.
    :raises RuntimeWarning: If salt, temp, or pres are outside the
        recommended bounds and chkbnd is True.
    :raises ValueError: If any of (drvs,drvt,drvp) are negative, or if
        drvt is nonzero and one of (drvs,drvp) is nonzero.
    """
    _chksalbnds(salt,temp,pres,chkbnd=chkbnd)
    # Temperature derivative is special
    if (drvs,drvt,drvp) == (0,1,0):
        g_tt = sal_g(0,2,0,salt,temp,pres,useext=useext)
        enth = -temp*g_tt
        return enth
    elif drvt != 0:
        errmsg = 'Derivatives {0} not recognized'.format((drvs,drvt,drvp))
        raise ValueError(errmsg)
    
    g = sal_g(drvs,0,drvp,salt,temp,pres,useext=useext)
    g_t = sal_g(drvs,1,drvp,salt,temp,pres,useext=useext)
    enth = g - temp*g_t
    return enth

def eq_entropy(drvs,drvt,drvp,salt,temp,pres,chkbnd=False,useext=False):
    """Calculate salt entropy with derivatives.
    
    Calculate the contribution of salt to the specific entropy of
    seawater or its derivatives with respect to salinity, temperature,
    and pressure. Only derivatives up to first order have been error
    checked.
    
    :arg int drvs: Number of salinity derivatives.
    :arg int drvt: Number of temperature derivatives.
    :arg int drvp: Number of pressure derivatives.
    :arg float salt: Absolute salinity in kg/kg.
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg bool useext: If False (default) then the polynomial is
        calculated from _GSCOEFFS; if True, from _GSCOEFFS_EXT.
    :returns: Entropy in units of
        (J/kg/K) / (kg/kg)^drvs / K^drvt / Pa^drvp.
    :raises ValueError: If temp or pres are nonpositive, or if salt is
        not between 0 and 1.
    :raises RuntimeWarning: If salt, temp, or pres are outside the
        recommended bounds and chkbnd is True.
    :raises ValueError: If any of (drvs,drvt,drvp) are negative.
    """
    _chksalbnds(salt,temp,pres,chkbnd=chkbnd)
    entr = -sal_g(drvs,drvt+1,drvp,salt,temp,pres,useext=useext)
    return entr

