"""Ice thermodynamic properties.

This module provides thermodynamic properties of ice derived from the
Gibbs free energy, with primary variables of temperature and pressure.
This module can also be called as a function,

    python ice_2.py

which will compare results from this module to reference values from
IAPWS 2006, table 6.

:Examples:

>>> cp(270.,1e5)
2073.47946211
>>> density(270.,1e5)
917.181167192
>>> kappat(270.,1e5)
1.17226047281e-10

:Functions:

* chempot: Chemical potential of ice.
* cp: Isobaric heat capacity of ice.
* density: Density of ice.
* enthalpy: Specific enthalpy of ice.
* entropy: Specific entropy of ice.
* expansion: Thermal expansion coefficient of ice.
* helmholtzenergy: Specific Helmholtz free energy of ice.
* internalenergy: Specific internal energy of ice.
* kappas: Isentropic compressibility of ice.
* kappat: Isothermal compressibility of ice.
* lapserate: Adiabatic lapse rate of ice.
* pcoefficient: Isochoric pressure coefficient of ice.
* specificvolume: Specific volume of ice.
* chkiapws06table6: Check module accuracy against IAPWS 2006 table 6.

"""

__all__ = ['chempot','cp','density','enthalpy','entropy','expansion',
    'helmholtzenergy','internalenergy','kappas','kappat','lapserate',
    'pcoefficient','specificvolume']

import constants0
import ice1

_CHKTOL = constants0.CHKTOL
_chkicebnds = constants0.chkicebnds
_ice_g = ice1.ice_g


## Thermodynamic properties
def chempot(temp,pres,chkbnd=False):
    """Calculate ice chemical potential.
    
    Calculate the chemical potential (specific Gibbs free energy) of ice
    from temperature and pressure.
    
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :returns: Chemical potential in J/kg.
    :raises ValueError: If temp or pres are nonpositive.
    :raises RuntimeWarning: If temp or pres are outside the recommended
        bounds and chkbnd is True.
    
    :Examples:
    
    >>> chempot(270.,1e5)
    -3786.74963128
    """
    _chkicebnds(temp,pres,chkbnd=chkbnd,stacklevel=2)
    g = _ice_g(0,0,temp,pres)
    return g

def cp(temp,pres,chkbnd=False):
    """Calculate ice isobaric heat capacity.
    
    Calculate the isobaric (constant pressure) heat capacity of ice from
    temperature and pressure.
    
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :returns: Heat capacity in J/kg/K.
    :raises ValueError: If temp or pres are nonpositive.
    :raises RuntimeWarning: If temp or pres are outside the recommended
        bounds and chkbnd is True.
    
    :Examples:
    
    >>> cp(270.,1e5)
    2073.47946211
    """
    _chkicebnds(temp,pres,chkbnd=chkbnd,stacklevel=2)
    g_tt = _ice_g(2,0,temp,pres)
    cp = -temp*g_tt
    return cp

def density(temp,pres,chkbnd=False):
    """Calculate ice density.
    
    Calculate the density of ice from temperature and pressure.
    
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :returns: Density in kg/m3.
    :raises ValueError: If temp or pres are nonpositive.
    :raises RuntimeWarning: If temp or pres are outside the recommended
        bounds and chkbnd is True.
    
    :Examples:
    
    >>> density(270.,1e5)
    917.181167192
    """
    _chkicebnds(temp,pres,chkbnd=chkbnd,stacklevel=2)
    g_p = _ice_g(0,1,temp,pres)
    dice = g_p**(-1)
    return dice

def enthalpy(temp,pres,chkbnd=False):
    """Calculate ice enthalpy.
    
    Calculate the specific enthalpy of ice from temperature and
    pressure.
    
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :returns: Pressure in Pa.
    :raises ValueError: If temp or pres are nonpositive.
    :raises RuntimeWarning: If temp or pres are outside the recommended
        bounds and chkbnd is True.
    
    :Examples:
    
    >>> enthalpy(270.,1e5)
    -339929.555499
    """
    _chkicebnds(temp,pres,chkbnd=chkbnd,stacklevel=2)
    g = _ice_g(0,0,temp,pres)
    g_t = _ice_g(1,0,temp,pres)
    h = g - temp*g_t
    return h

def entropy(temp,pres,chkbnd=False):
    """Calculate ice entropy.
    
    Calculate the specific entropy of ice from temperature and pressure.
    
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :returns: Entropy in J/kg/K.
    :raises ValueError: If temp or pres are nonpositive.
    :raises RuntimeWarning: If temp or pres are outside the recommended
        bounds and chkbnd is True.
    
    :Examples:
    
    >>> entropy(270.,1e5)
    -1244.97335506
    """
    _chkicebnds(temp,pres,chkbnd=chkbnd,stacklevel=2)
    g_t = _ice_g(1,0,temp,pres)
    s = -g_t
    return s

def expansion(temp,pres,chkbnd=False):
    """Calculate ice thermal expansion coefficient.
    
    Calculate the thermal expansion coefficient of ice from temperature
    and pressure.
    
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :returns: Expansion coefficient in 1/K.
    :raises ValueError: If temp or pres are nonpositive.
    :raises RuntimeWarning: If temp or pres are outside the recommended
        bounds and chkbnd is True.
    
    :Examples:
    
    >>> expansion(270.,1e5)
    1.58309329594e-04
    """
    _chkicebnds(temp,pres,chkbnd=chkbnd,stacklevel=2)
    g_p = _ice_g(0,1,temp,pres)
    g_tp = _ice_g(1,1,temp,pres)
    alpha = g_tp / g_p
    return alpha

def helmholtzenergy(temp,pres,chkbnd=False):
    """Calculate ice Helmholtz free energy.
    
    Calculate the specific Helmholtz free energy of ice from temperature
    and pressure.
    
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :returns: Helmholtz free energy in J/kg.
    :raises ValueError: If temp or pres are nonpositive.
    :raises RuntimeWarning: If temp or pres are outside the recommended
        bounds and chkbnd is True.
    
    :Examples:
    
    >>> helmholtzenergy(270.,1e5)
    -3895.77934490
    """
    _chkicebnds(temp,pres,chkbnd=chkbnd,stacklevel=2)
    g = _ice_g(0,0,temp,pres)
    g_p = _ice_g(0,1,temp,pres)
    f = g - pres*g_p
    return f

def internalenergy(temp,pres,chkbnd=False):
    """Calculate ice internal energy.
    
    Calculate the specific internal energy of ice from temperature and
    pressure.
    
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :returns: Internal energy in J/kg.
    :raises ValueError: If temp or pres are nonpositive.
    :raises RuntimeWarning: If temp or pres are outside the recommended
        bounds and chkbnd is True.
    
    :Examples:
    
    >>> internalenergy(270.,1e5)
    -340038.585212
    """
    _chkicebnds(temp,pres,chkbnd=chkbnd,stacklevel=2)
    g = _ice_g(0,0,temp,pres)
    g_t = _ice_g(1,0,temp,pres)
    g_p = _ice_g(0,1,temp,pres)
    u = g - pres*g_p - temp*g_t
    return u

def kappas(temp,pres,chkbnd=False):
    """Calculate ice isentropic compressibility.
    
    Calculate the isentropic compressibility of ice from temperature and
    pressure.
    
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :returns: Compressibility in 1/Pa.
    :raises ValueError: If temp or pres are nonpositive.
    :raises RuntimeWarning: If temp or pres are outside the recommended
        bounds and chkbnd is True.
    
    :Examples:
    
    >>> kappas(270.,1e5)
    1.13667916416e-10
    """
    _chkicebnds(temp,pres,chkbnd=chkbnd,stacklevel=2)
    g_p = _ice_g(0,1,temp,pres)
    g_tt = _ice_g(2,0,temp,pres)
    g_tp = _ice_g(1,1,temp,pres)
    g_pp = _ice_g(0,2,temp,pres)
    kappa = (g_tp**2 - g_tt*g_pp) / (g_p*g_tt)
    return kappa

def kappat(temp,pres,chkbnd=False):
    """Calculate ice isothermal compressibility.
    
    Calculate the isothermal compressibility of ice from temperature and
    pressure.
    
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :returns: Compressibility in 1/Pa.
    :raises ValueError: If temp or pres are nonpositive.
    :raises RuntimeWarning: If temp or pres are outside the recommended
        bounds and chkbnd is True.
    
    :Examples:
    
    >>> kappat(270.,1e5)
    1.17226047281e-10
    """
    _chkicebnds(temp,pres,chkbnd=chkbnd,stacklevel=2)
    g_p = _ice_g(0,1,temp,pres)
    g_pp = _ice_g(0,2,temp,pres)
    kappa = -g_pp / g_p
    return kappa

def lapserate(temp,pres,chkbnd=False):
    """Calculate ice adiabatic lapse rate.
    
    Calculate the adiabatic lapse rate of ice from temperature and
    pressure.
    
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :returns: Lapse rate in K/Pa.
    :raises ValueError: If temp or pres are nonpositive.
    :raises RuntimeWarning: If temp or pres are outside the recommended
        bounds and chkbnd is True.
    
    :Examples:
    
    >>> lapserate(270.,1e5)
    2.24758128545e-08
    """
    _chkicebnds(temp,pres,chkbnd=chkbnd,stacklevel=2)
    g_tt = _ice_g(2,0,temp,pres)
    g_tp = _ice_g(1,1,temp,pres)
    gamma = -g_tp / g_tt
    return gamma

def pcoefficient(temp,pres,chkbnd=False):
    """Calculate ice pressure coefficient.
    
    Calculate the isochoric pressure coefficient (rate of change of
    pressure with temperature at constant volume) of ice from
    temperature and pressure.
    
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :returns: Pressure coefficient in Pa/K.
    :raises ValueError: If temp or pres are nonpositive.
    :raises RuntimeWarning: If temp or pres are outside the recommended
        bounds and chkbnd is True.
    
    :Examples:
    
    >>> pcoefficient(270.,1e5)
    1350462.06254
    """
    _chkicebnds(temp,pres,chkbnd=chkbnd,stacklevel=2)
    g_tp = _ice_g(1,1,temp,pres)
    g_pp = _ice_g(0,2,temp,pres)
    pcoef = -g_tp / g_pp
    return pcoef

def specificvolume(temp,pres,chkbnd=False):
    """Calculate ice specific volume.
    
    Calculate the specific volume of ice from temperature and pressure.
    
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :returns: Specific volume in m3/kg.
    :raises ValueError: If temp or pres are nonpositive.
    :raises RuntimeWarning: If temp or pres are outside the recommended
        bounds and chkbnd is True.
    
    :Examples:
    
    >>> specificvolume(270.,1e5)
    1.09029713624e-03
    """
    _chkicebnds(temp,pres,chkbnd=chkbnd,stacklevel=2)
    v = _ice_g(0,1,temp,pres)
    return v


'''
### Functions to check values
def chk_iapws06_table6(printresult=True,tol=_CHKTOL):
    """Check accuracy against IAPWS 2006 table 6.
    
    Evaluate the functions in this module and compare to reference values of
    thermodynamic properties (e.g. heat capacity, lapse rate) in IAPWS 2006,
    table 6.
    
    Args:
        printresult (boolean, optional): If True (default) and any results are
            outside of the given tolerance, then the function name, reference
            value, resulting value, and relative error are printed.
        tol (float, optional): Tolerance to use when choosing to print results;
            default 1e-8.
    
    Returns:
        chkders (list): Dictionaries containing the comparisons for derivatives
            of the Gibbs free energy. Each dictionary contains the values of the
            derivatives 'ders'; the physical variables 'args'; the reference
            values 'refs'; and the results from this module 'results'.
        chkfuns (list): Dictionaries containing the comparisons for
            thermodynamic properties. Each dictionary contains the names of the
            functions 'names'; the physical variables 'args'; the reference
            values 'refs'; and the results from this module 'results'.
    """
    
    from values_common import runcheck
    
    FUNS = (ice_enthalpy,ice_helmholtz_energy,ice_internal_energy,
        ice_entropy,ice_cp,ice_density,ice_expansion,ice_p_coefficient,
        ice_kappa_t,ice_kappa_s)
    NAMES = ('enthalpy','helmholtz_energy','internal_energy','entropy',
        'cp','density','expansion','p_coefficient','kappa_t','kappa_s')
    DERS = ((0,0),(1,0),(0,1),(2,0),(1,1),(0,2))
    
    # Create dictionaries with the reference data
    checkDer1 = {'modname': 'ice1',
        'type': 'der',
        'args': (273.16,611.657),
        'funs': _ice_g,
        'names': '_ice_g',
        'ders': DERS,
        'refs': (0.611784135,0.122069433940e+4,0.109085812737e-2,
            -0.767602985875e+1,0.174387964700e-6,-0.128495941571e-12)}
    
    checkDer2 = {'modname': 'ice1',
        'type': 'der',
        'args': (273.152519,101325.),
        'funs': _ice_g,
        'names': '_ice_g',
        'ders': DERS,
        'refs': (0.10134274069e+3,0.122076932550e+4,0.109084388214e-2,
            -0.767598233365e+1,0.174362219972e-6,-0.128485364928e-12)}
    
    checkDer3 = {'modname': 'ice1',
        'type': 'der',
        'args': (100.,100e6),
        'funs': _ice_g,
        'names': '_ice_g',
        'ders': DERS,
        'refs': (-0.222296513088e+6,0.261195122589e+4,0.106193389260e-2,
            -0.866333195517e+1,0.274505162488e-7,-0.941807981761e-13)}
    
    checkName1 = {'modname': 'ice_2',
        'type': 'fun',
        'args': (273.16,611.657),
        'funs': FUNS,
        'names': NAMES,
        'refs': (-0.333444253966e+6,-0.55446875e-1,-0.333444921197e+6,
            -0.122069433940e+4,0.209678431622e+4,0.916709492200e+3,
            0.159863102566e-3,0.135714764659e+7,0.117793449348e-9,
            0.114161597779e-9)}
    
    checkName2 = {'modname': 'ice_2',
        'type': 'fun',
        'args': (273.152519,101325.),
        'funs': FUNS,
        'names': NAMES,
        'refs': (-0.333354873637e+6,-0.918701567e+1,-0.333465403393e+6,
            -0.122076932550e+4,0.209671391024e+4,0.916721463419e+3,
            0.159841589458e-3,0.135705899321e+7,0.117785291765e-9,
            0.114154442556e-9)}
    
    checkName3 = {'modname': 'ice_2',
        'type': 'fun',
        'args': (100.,100e6),
        'funs': FUNS,
        'names': NAMES,
        'refs': (-0.483491635676e+6,-0.328489902347e+6,
            -0.589685024936e+6,-0.261195122589e+4,0.866333195517e+3,
            0.941678203297e+3,0.258495528207e-4,0.291466166994e+6,
            0.886880048115e-10,0.886060982687e-10)}
    
    # Check each dictionary
    runcheck(checkDer1,printresult=printresult,tol=tol)
    runcheck(checkDer2,printresult=printresult,tol=tol)
    runcheck(checkDer3,printresult=printresult,tol=tol)
    
    runcheck(checkName1,printresult=printresult,tol=tol)
    runcheck(checkName2,printresult=printresult,tol=tol)
    runcheck(checkName3,printresult=printresult,tol=tol)
    
    # Return the dictionaries
    chkders = [checkDer1, checkDer2, checkDer3]
    chkfuns = [checkName1, checkName2, checkName3]
    return chkders, chkfuns



### Main function: Check tables
if __name__ == '__main__':
    chkders, chkfuns = chk_iapws06_table6();
'''

