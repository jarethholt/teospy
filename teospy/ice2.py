"""Ice thermodynamic properties.

This module provides thermodynamic properties of ice derived from the
Gibbs free energy, with primary variables of temperature and pressure.

This module can also be called as a function::

    python ice_2.py

which will compare results from this module to reference values from
IAPWS 2006, table 6.

:Examples:

>>> cp(270.,1e5)
2073.47946211
>>> density(270.,1e5)
917.181167192
>>> kappa_t(270.,1e5)
1.17226047281e-10

:Functions:

* :func:`chempot`: Chemical potential of ice.
* :func:`cp`: Isobaric heat capacity of ice.
* :func:`density`: Density of ice.
* :func:`enthalpy`: Specific enthalpy of ice.
* :func:`entropy`: Specific entropy of ice.
* :func:`expansion`: Thermal expansion coefficient of ice.
* :func:`helmholtzenergy`: Specific Helmholtz free energy of ice.
* :func:`internalenergy`: Specific internal energy of ice.
* :func:`kappa_s`: Isentropic compressibility of ice.
* :func:`kappa_t`: Isothermal compressibility of ice.
* :func:`lapserate`: Adiabatic lapse rate of ice.
* :func:`pcoefficient`: Isochoric pressure coefficient of ice.
* :func:`specificvolume`: Specific volume of ice.
* :func:`chkiapws06table6`: Check module accuracy against IAPWS 2006
  table 6.

"""

__all__ = ['chempot','cp','density','enthalpy','entropy','expansion',
    'helmholtzenergy','internalenergy','kappa_s','kappa_t','lapserate',
    'pcoefficient','specificvolume',
    'chkiapws06table6']

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
    _chkicebnds(temp,pres,chkbnd=chkbnd)
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
    _chkicebnds(temp,pres,chkbnd=chkbnd)
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
    _chkicebnds(temp,pres,chkbnd=chkbnd)
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
    _chkicebnds(temp,pres,chkbnd=chkbnd)
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
    _chkicebnds(temp,pres,chkbnd=chkbnd)
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
    _chkicebnds(temp,pres,chkbnd=chkbnd)
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
    _chkicebnds(temp,pres,chkbnd=chkbnd)
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
    _chkicebnds(temp,pres,chkbnd=chkbnd)
    g = _ice_g(0,0,temp,pres)
    g_t = _ice_g(1,0,temp,pres)
    g_p = _ice_g(0,1,temp,pres)
    u = g - pres*g_p - temp*g_t
    return u

def kappa_s(temp,pres,chkbnd=False):
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
    
    >>> kappa_s(270.,1e5)
    1.13667916416e-10
    """
    _chkicebnds(temp,pres,chkbnd=chkbnd)
    g_p = _ice_g(0,1,temp,pres)
    g_tt = _ice_g(2,0,temp,pres)
    g_tp = _ice_g(1,1,temp,pres)
    g_pp = _ice_g(0,2,temp,pres)
    kappa = (g_tp**2 - g_tt*g_pp) / (g_p*g_tt)
    return kappa

def kappa_t(temp,pres,chkbnd=False):
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
    
    >>> kappa_t(270.,1e5)
    1.17226047281e-10
    """
    _chkicebnds(temp,pres,chkbnd=chkbnd)
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
    _chkicebnds(temp,pres,chkbnd=chkbnd)
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
    _chkicebnds(temp,pres,chkbnd=chkbnd)
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
    _chkicebnds(temp,pres,chkbnd=chkbnd)
    v = _ice_g(0,1,temp,pres)
    return v



## Functions to check values
def chkiapws06table6(printresult=True,chktol=_CHKTOL):
    """Check accuracy against IAPWS 2006 table 6.
    
    Evaluate the functions in this module and compare to reference
    values of thermodynamic properties (e.g. heat capacity, lapse rate)
    in IAPWS 2006, table 6.
    
    :arg bool printresult: If True (default) and any results are outside
        of the given tolerance, then the function name, reference value,
        result value, and relative error are printed.
    :arg float chktol: Tolerance to use when choosing to print results
        (default _CHKTOL).
    :returns: :class:`~tester.Tester` instances containing the
        functions, arguments, reference values, results, and relative
        errors from the tests. The first instance involves derivatives
        of ice_g whereas the second tests the other thermodynamic
        functions.
    :rtype: (Tester,Tester)
    """
    from tester import Tester
    fargs0 = (273.16,611.657)
    fargs1 = (273.152519,101325.)
    fargs2 = (100.,1e8)
    propfargs = [fargs0,fargs1,fargs2]
    ders = [(0,0), (1,0), (0,1), (2,0), (1,1), (0,2)]
    
    # Tester instance for derivatives of ice_g
    derfuns = _ice_g
    derfnames = 'ice_g'
    # Derivatives change before arguments do here
    derfargs = [(der+fargs) for fargs in propfargs for der in ders]
    derargfmt = '({0:1g},{1:1g},{2:7.3f},{3:7g})'
    derrefs = [0.611784135,0.122069433940e+4,0.109085812737e-2,
        -0.767602985875e+1,0.174387964700e-6,-0.128495941571e-12,
        0.10134274069e+3,0.122076932550e+4,0.109084388214e-2,-0.767598233365e+1,
        0.174362219972e-6,-0.128485364928e-12,-0.222296513088e+6,
        0.261195122589e+4,0.106193389260e-2,-0.866333195517e+1,
        0.274505162488e-7,-0.941807981761e-13]
    header = 'Ice Gibbs energy derivatives'
    dertest = Tester(derfuns,derfargs,derrefs,derfnames,derargfmt,header=header)
    
    # Tester instance for other ice properties
    propfuns = [enthalpy,helmholtzenergy,internalenergy,entropy,cp,density,
        expansion,pcoefficient,kappa_t,kappa_s]
    propfnames = ['enthalpy','helmholtzenergy','internalenergy','entropy','cp',
        'density','expansion','pcoefficient','kappa_t','kappa_s']
    propargfmt = '({0:7.3f},{1:7g})'
    proprefs = [
        [-0.333444253966e+6,-0.333354873637e+6,-0.483491635676e+6],
        [-0.55446875e-1,-0.918701567e+1,-0.328489902347e+6],
        [-0.333444921197e+6,-0.333465403393e+6,-0.589685024936e+6],
        [-0.122069433940e+4,-0.122076932550e+4,-0.261195122589e+4],
        [0.209678431622e+4,0.209671391024e+4,0.866333195517e+3],
        [0.916709492200e+3,0.916721463419e+3,0.941678203297e+3],
        [0.159863102566e-3,0.159841589458e-3,0.258495528207e-4],
        [0.135714764659e+7,0.135705899321e+7,0.291466166994e+6],
        [0.117793449348e-9,0.117785291765e-9,0.886880048115e-10],
        [0.114161597779e-9,0.114154442556e-9,0.886060982687e-10]
    ]
    header = 'Ice thermodynamic properties'
    proptest = Tester(propfuns,propfargs,proprefs,propfnames,propargfmt,
        header=header)
    
    # Run Tester instances and print results
    dertest.run()
    proptest.run()
    if printresult:
        dertest.printresults(chktol=chktol)
        proptest.printresults(chktol=chktol)
    return dertest, proptest


## Main function: Check tables
if __name__ == '__main__':
    dertest, proptest = chkiapws06table6();

