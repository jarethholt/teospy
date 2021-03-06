"""Fluid water thermodynamic properties.

This module provides thermodynamic properties of pure fluid water
(liquid water or water vapour) derived from the Helmholtz free energy,
with primary variables temperature and fluid density. This module can
also be called as a function::

    python flu_2.py

which will compare results from this module to the reference values in
IAPWS 1995, table 7.

:Examples:

>>> cv(300.,1000.)
4105.20614776
>>> lapserate(300.,1000.)
2.03694039833e-08
>>> pressure(300.,1000.)
7833001.35597

:Functions:

* :func:`cp`: Isobaric heat capacity of fluid water.
* :func:`cv`: Isochoric heat capacity of fluid water.
* :func:`enthalpy`: Specific enthalpy of fluid water.
* :func:`entropy`: Specific entropy of fluid water.
* :func:`expansion`: Thermal expansion coefficient of fluid water.
* :func:`gibbsenergy`: Specific Gibbs free energy of fluid water.
* :func:`internalenergy`: Specific internal energy of fluid water.
* :func:`kappa_s`: Isentropic compressibility of fluid water.
* :func:`kappa_t`: Isothermal compressibility of fluid water.
* :func:`lapserate`: Adiabatic lapse rate of fluid water.
* :func:`pressure`: Pressure of fluid water.
* :func:`soundspeed`: Speed of sound in fluid water.
* :func:`eq_chempot`: Fluid water chemical potential with derivatives.
* :func:`eq_enthalpy`: Fluid water enthalpy with derivatives.
* :func:`eq_entropy`: Fluid water entropy with derivatives.
* :func:`eq_pressure`: Fluid water pressure with derivatives.
* :func:`chkiapws95table7`: Check module accuracy against IAPWS 1995
  table 7.

"""

__all__ = ['cp','cv','enthalpy','entropy','expansion','gibbsenergy',
    'internalenergy','kappa_s','kappa_t','lapserate','pressure','soundspeed',
    'eq_chempot','eq_enthalpy','eq_entropy','eq_pressure',
    'chkiapws95table7']

from teospy import constants0
from teospy import flu1

_CHKTOL = constants0.CHKTOL
_chkflubnds = constants0.chkflubnds
_flu_f = flu1.flu_f


## Thermodynamic properties
def cp(temp,dflu,chkbnd=False):
    """Calculate fluid water isobaric heat capacity.
    
    Calculate the isobaric (constant pressure) heat capacity of fluid
    water from the temperature and fluid density.
    
    :arg float temp: Temperature in K.
    :arg float dflu: Fluid water density in kg/m3.
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :returns: Heat capacity in J/kg/K.
    :raises ValueError: If temp or dflu are nonpositive.
    :raises RuntimeWarning: If temp or dflu are outside the recommended
        bounds and chkbnd is True.
    
    :Examples:
    
    >>> cp(300.,1000.)
    4159.37519963
    """
    _chkflubnds(temp,dflu,chkbnd=chkbnd)
    f_d = _flu_f(0,1,temp,dflu)
    f_tt = _flu_f(2,0,temp,dflu)
    f_td = _flu_f(1,1,temp,dflu)
    f_dd = _flu_f(0,2,temp,dflu)
    cp = f_td**2 - f_dd*f_tt - 2*f_d*f_tt/dflu
    cp /= f_dd + 2*f_d/dflu
    cp *= temp
    return cp

def cv(temp,dflu,chkbnd=False):
    """Calculate fluid water isochoric heat capacity.
    
    Calculate the isochoric (constant volume) heat capacity of fluid
    water from temperature and fluid density.
    
    :arg float temp: Temperature in K.
    :arg float dflu: Fluid water density in kg/m3.
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :returns: Heat capacity in J/kg/K.
    :raises ValueError: If temp or dflu are nonpositive.
    :raises RuntimeWarning: If temp or dflu are outside the recommended
        bounds and chkbnd is True.
    
    :Examples:
    
    >>> cv(300.,1000.)
    4105.20614776
    """
    _chkflubnds(temp,dflu,chkbnd=chkbnd)
    f_tt = _flu_f(2,0,temp,dflu)
    cv = -temp * f_tt
    return cv

def enthalpy(temp,dflu,chkbnd=False):
    """Calculate fluid water enthalpy.
    
    Calculate the specific enthalpy of fluid water from temperature and
    fluid density.
    
    :arg float temp: Temperature in K.
    :arg float dflu: Fluid water density in kg/m3.
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :returns: Enthalpy in J/kg.
    :raises ValueError: If temp or dflu are nonpositive.
    :raises RuntimeWarning: If temp or dflu are outside the recommended
        bounds and chkbnd is True.
    
    :Examples:
    
    >>> enthalpy(300.,1000.)
    119752.511434
    """
    _chkflubnds(temp,dflu,chkbnd=chkbnd)
    f = _flu_f(0,0,temp,dflu)
    f_t = _flu_f(1,0,temp,dflu)
    f_d = _flu_f(0,1,temp,dflu)
    h = f - temp*f_t + dflu*f_d
    return h

def entropy(temp,dflu,chkbnd=False):
    """Calculate fluid water entropy.
    
    Calculate the specific entropy of fluid water from temperature and
    fluid density.
    
    :arg float temp: Temperature in K.
    :arg float dflu: Fluid water density in kg/m3.
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :returns: Entropy in J/kg/K.
    :raises ValueError: If temp or dflu are nonpositive.
    :raises RuntimeWarning: If temp or dflu are outside the recommended
        bounds and chkbnd is True.
    
    :Examples:
    
    >>> entropy(300.,1000.)
    390.904170767
    """
    _chkflubnds(temp,dflu,chkbnd=chkbnd)
    f_t = _flu_f(1,0,temp,dflu)
    s = -f_t
    return s

def expansion(temp,dflu,chkbnd=False):
    """Calculate fluid water thermal expansion coefficient.
    
    Calculate the thermal expansion coefficient of fluid water from
    temperature and fluid density.
    
    :arg float temp: Temperature in K.
    :arg float dflu: Fluid water density in kg/m3.
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :returns: Expansion coefficient in 1/K.
    :raises ValueError: If temp or dflu are nonpositive.
    :raises RuntimeWarning: If temp or dflu are outside the recommended
        bounds and chkbnd is True.
    
    :Examples:
    
    >>> expansion(300.,1000.)
    2.82413312531e-04
    """
    _chkflubnds(temp,dflu,chkbnd=chkbnd)
    f_d = _flu_f(0,1,temp,dflu)
    f_td = _flu_f(1,1,temp,dflu)
    f_dd = _flu_f(0,2,temp,dflu)
    alpha = f_td / (2*f_d + dflu*f_dd)
    return alpha

def gibbsenergy(temp,dflu,chkbnd=False):
    """Calculate fluid water Gibbs free energy.
    
    Calculate the specific Gibbs free energy of fluid water from
    temperature and fluid density.
    
    :arg float temp: Temperature in K.
    :arg float dflu: Fluid water density in kg/m3.
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :returns: Gibbs energy in J/kg.
    :raises ValueError: If temp or dflu are nonpositive.
    :raises RuntimeWarning: If temp or dflu are outside the recommended
        bounds and chkbnd is True.
    
    :Examples:
    
    >>> gibbsenergy(300.,1000.)
    2481.2602039
    """
    _chkflubnds(temp,dflu,chkbnd=chkbnd)
    f = _flu_f(0,0,temp,dflu)
    f_d = _flu_f(0,1,temp,dflu)
    g = f + dflu*f_d
    return g

def internalenergy(temp,dflu,chkbnd=False):
    """Calculate fluid water internal energy.
    
    Calculate the specific internal energy of fluid water from
    temperature and fluid density.
    
    :arg float temp: Temperature in K.
    :arg float dflu: Fluid water density in kg/m3.
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :returns: Internal energy in J/kg.
    :raises ValueError: If temp or dflu are nonpositive.
    :raises RuntimeWarning: If temp or dflu are outside the recommended
        bounds and chkbnd is True.
    
    :Examples:
    
    >>> internalenergy(300.,1000.)
    111919.510078
    """
    _chkflubnds(temp,dflu,chkbnd=chkbnd)
    f = _flu_f(0,0,temp,dflu)
    f_t = _flu_f(1,0,temp,dflu)
    u = f - temp*f_t
    return u

def kappa_s(temp,dflu,chkbnd=False):
    """Calculate fluid water isentropic compressibility.
    
    Calculate the isentropic compressibility of fluid water from
    temperature and fluid density.
    
    :arg float temp: Temperature in K.
    :arg float dflu: Fluid water density in kg/m3.
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :returns: Compressibility in 1/Pa.
    :raises ValueError: If temp or dflu are nonpositive.
    :raises RuntimeWarning: If temp or dflu are outside the recommended
        bounds and chkbnd is True.
    
    :Examples:
    
    >>> kappa_s(300.,1000.)
    4.35960581171e-10
    """
    _chkflubnds(temp,dflu,chkbnd=chkbnd)
    f_d = _flu_f(0,1,temp,dflu)
    f_tt = _flu_f(2,0,temp,dflu)
    f_td = _flu_f(1,1,temp,dflu)
    f_dd = _flu_f(0,2,temp,dflu)
    denom = dflu**2 * (f_tt*(2*f_d + dflu*f_dd) - dflu*f_td**2)
    kappa = f_tt / denom
    return kappa

def kappa_t(temp,dflu,chkbnd=False):
    """Calculate fluid water isothermal compressibility.
    
    Calculate the isothermal compressibility of fluid water from
    temperature and fluid density.
    
    :arg float temp: Temperature in K.
    :arg float dflu: Fluid water density in kg/m3.
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :returns: Compressibility in 1/Pa.
    :raises ValueError: If temp or dflu are nonpositive.
    :raises RuntimeWarning: If temp or dflu are outside the recommended
        bounds and chkbnd is True.
    
    :Examples:
    
    >>> kappa_t(300.,1000.)
    4.41713172024e-10
    """
    _chkflubnds(temp,dflu,chkbnd=chkbnd)
    f_d = _flu_f(0,1,temp,dflu)
    f_dd = _flu_f(0,2,temp,dflu)
    denom = dflu**2 * (2*f_d + dflu*f_dd)
    kappa = denom**(-1)
    return kappa

def lapserate(temp,dflu,chkbnd=False):
    """Calculate fluid water adiabatic lapse rate.
    
    Calculate the adiabatic lapse rate of fluid water from temperature
    and fluid density.
    
    :arg float temp: Temperature in K.
    :arg float dflu: Fluid water density in kg/m3.
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :returns: Lapse rate in K/Pa.
    :raises ValueError: If temp or dflu are nonpositive.
    :raises RuntimeWarning: If temp or dflu are outside the recommended
        bounds and chkbnd is True.
    
    :Examples:
    
    >>> lapserate(300.,1000.)
    2.03694039833e-08
    """
    _chkflubnds(temp,dflu,chkbnd=chkbnd)
    f_d = _flu_f(0,1,temp,dflu)
    f_tt = _flu_f(2,0,temp,dflu)
    f_td = _flu_f(1,1,temp,dflu)
    f_dd = _flu_f(0,2,temp,dflu)
    denom = f_tt*(2*f_d + dflu*f_dd) - dflu*f_td**2
    gamma = -f_td / denom / dflu
    return gamma

def pressure(temp,dflu,chkbnd=False):
    """Calculate fluid water pressure.
    
    Calculate the pressure in fluid water from temperature and fluid
    density.
    
    :arg float temp: Temperature in K.
    :arg float dflu: Fluid water density in kg/m3.
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :returns: Pressure in Pa.
    :raises ValueError: If temp or dflu are nonpositive.
    :raises RuntimeWarning: If temp or dflu are outside the recommended
        bounds and chkbnd is True.
    
    :Examples:
    
    >>> pressure(300.,1000.)
    7833001.35597
    """
    _chkflubnds(temp,dflu,chkbnd=chkbnd)
    f_d = _flu_f(0,1,temp,dflu)
    p = dflu**2 * f_d
    return p

def soundspeed(temp,dflu,chkbnd=False):
    """Calculate fluid water sound speed.
    
    Calculate the speed of sound in fluid water from temperature and
    fluid density.
    
    :arg float temp: Temperature in K.
    :arg float dflu: Fluid water density in kg/m3.
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :returns: Sound speed in m/s.
    :raises ValueError: If temp or dflu are nonpositive.
    :raises RuntimeWarning: If temp or dflu are outside the recommended
        bounds and chkbnd is True.
    
    :Examples:
    
    >>> soundspeed(300.,1000.)
    1514.52479780
    """
    _chkflubnds(temp,dflu,chkbnd=chkbnd)
    f_d = _flu_f(0,1,temp,dflu)
    f_tt = _flu_f(2,0,temp,dflu)
    f_td = _flu_f(1,1,temp,dflu)
    f_dd = _flu_f(0,2,temp,dflu)
    c_sq = dflu**2*(f_dd - f_td**2/f_tt) + 2*dflu*f_d
    c = c_sq**(.5)
    return c


## Equilibrium calculation auxiliary functions
def eq_chempot(drvt,drvd,temp,dflu,chkbnd=False):
    """Calculate fluid water potential with derivatives.
    
    Calculate the chemical potential of fluid water or its derivatives
    with respect to temperature and fluid density.
    
    :arg int drvt: Number of temperature derivatives.
    :arg int drvd: Number of density derivatives.
    :arg float temp: Temperature in K.
    :arg float dflu: Fluid water density in kg/m3.
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :returns: Chemical potential in units of
        (J/kg) / K^drvt / (kg/m3)^drvd.
    :raises ValueError: If temp or dflu are nonpositive.
    :raises RuntimeWarning: If temp or dflu are outside the recommended
        bounds and chkbnd is True.
    :raises ValueError: If either of (drvt,drvd) are <0 or if
        drvt+drvd>1.
    """
    _chkflubnds(temp,dflu,chkbnd=chkbnd)
    
    # Run through derivative cases
    if (drvt,drvd) == (0,0):
        f = _flu_f(0,0,temp,dflu)
        f_d = _flu_f(0,1,temp,dflu)
        gflu = f + dflu*f_d
    elif (drvt,drvd) == (1,0):
        f_t = _flu_f(1,0,temp,dflu)
        f_td = _flu_f(1,1,temp,dflu)
        gflu = f_t + dflu*f_td
    elif (drvt,drvd) == (0,1):
        f_d = _flu_f(0,1,temp,dflu)
        f_dd = _flu_f(0,2,temp,dflu)
        gflu = 2*f_d + dflu*f_dd
    else:
        errmsg = 'Derivatives {0} not recognized'.format((drvt,drvd))
        raise ValueError(errmsg)
    return gflu

def eq_enthalpy(drvt,drvd,temp,dflu,chkbnd=False):
    """Calculate fluid water enthalpy with derivatives.
    
    Calculate the specific enthalpy of fluid water or its derivatives with
    respect to temperature and fluid density.
    
    :arg int drvt: Number of temperature derivatives.
    :arg int drvd: Number of density derivatives.
    :arg float temp: Temperature in K.
    :arg float dflu: Fluid water density in kg/m3.
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :returns: Enthalpy in units of
        (J/kg) / K^drvt / (kg/m3)^drvd.
    :raises ValueError: If temp or dflu are nonpositive.
    :raises RuntimeWarning: If temp or dflu are outside the recommended
        bounds and chkbnd is True.
    :raises ValueError: If either of (drvt,drvd) are <0 or if
        drvt+drvd>1.
    """
    _chkflubnds(temp,dflu,chkbnd=chkbnd)
    
    # Run through derivative cases
    if (drvt,drvd) == (0,0):
        f = _flu_f(0,0,temp,dflu)
        f_t = _flu_f(1,0,temp,dflu)
        f_d = _flu_f(0,1,temp,dflu)
        hflu = f - temp*f_t + dflu*f_d
    elif (drvt,drvd) == (1,0):
        f_tt = _flu_f(2,0,temp,dflu)
        f_td = _flu_f(1,1,temp,dflu)
        hflu = -temp*f_tt + dflu*f_td
    elif (drvt,drvd) == (0,1):
        f_d = _flu_f(0,1,temp,dflu)
        f_td = _flu_f(1,1,temp,dflu)
        f_dd = _flu_f(0,2,temp,dflu)
        hflu = 2*f_d + dflu*f_dd - temp*f_td
    else:
        errmsg = 'Derivatives {0} not recognized'.format((drvt,drvd))
        raise ValueError(errmsg)
    return hflu

def eq_entropy(drvt,drvd,temp,dflu,chkbnd=False):
    """Calculate fluid water entropy with derivatives.
    
    Calculate the specific entropy of fluid water or its derivatives with
    respect to temperature and fluid density.
    
    :arg int drvt: Number of temperature derivatives.
    :arg int drvd: Number of density derivatives.
    :arg float temp: Temperature in K.
    :arg float dflu: Fluid water density in kg/m3.
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :returns: Entropy in units of
        (J/kg/K) / K^drvt / (kg/m3)^drvd.
    :raises ValueError: If temp or dflu are nonpositive.
    :raises RuntimeWarning: If temp or dflu are outside the recommended
        bounds and chkbnd is True.
    :raises ValueError: If either of (drvt,drvd) are <0 or if
        drvt+drvd>1.
    """
    _chkflubnds(temp,dflu,chkbnd=chkbnd)
    entr = -_flu_f(drvt+1,drvd,temp,dflu)
    return entr

def eq_pressure(drvt,drvd,temp,dflu,chkbnd=False):
    """Calculate fluid water pressure with derivatives.
    
    Calculate the pressure of fluid water or its derivatives with
    respect to temperature and fluid density.
    
    :arg int drvt: Number of temperature derivatives.
    :arg int drvd: Number of density derivatives.
    :arg float temp: Temperature in K.
    :arg float dflu: Fluid water density in kg/m3.
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :returns: Pressure in units of
        Pa / K^drvt / (kg/m3)^drvd.
    :raises ValueError: If temp or dflu are nonpositive.
    :raises RuntimeWarning: If temp or dflu are outside the recommended
        bounds and chkbnd is True.
    :raises ValueError: If either of (drvt,drvd) are <0 or if
        drvt+drvd>1.
    """
    _chkflubnds(temp,dflu,chkbnd=chkbnd)
    
    # Run through derivative cases
    if (drvt,drvd) == (0,0):
        f_d = _flu_f(0,1,temp,dflu)
        pres = dflu**2 * f_d
    elif (drvt,drvd) == (1,0):
        f_td = _flu_f(1,1,temp,dflu)
        pres = dflu**2 * f_td
    elif (drvt,drvd) == (0,1):
        f_d = _flu_f(0,1,temp,dflu)
        f_dd = _flu_f(0,2,temp,dflu)
        pres = 2*dflu*f_d + dflu**2 * f_dd
    else:
        errmsg = 'Derivatives {0} not recognized'.format((drvt,drvd))
        raise ValueError(errmsg)
    return pres


## Functions to check values
def chkiapws95table7(printresult=True,chktol=_CHKTOL):
    """Check accuracy against IAPWS 1995 table 7.
    
    Evaluate the functions in this module and compare to reference
    values of thermodynamic properties (e.g. heat capacity, pressure) in
    IAPWS 1995, table 7.
    
    :arg bool printresult: If True (default) and any results are outside
        of the given tolerance, then the function name, reference value,
        result value, and relative error are printed.
    :arg float chktol: Tolerance to use when choosing to print results
        (default _CHKTOL).
    :returns: :class:`~teospy.tests.tester.Tester` instance containing
        the functions, arguments, reference values, results, and
        relative errors from the tests.
    """
    from teospy.tests.tester import Tester
    funs = [pressure,cv,soundspeed,entropy]
    fnames = ['pressure','cv','soundspeed','entropy']
    fargs = [
        (300.,0.9965560e3), (300.,0.1005308e4), (300.,0.1188202e4),
        (500.,0.4350000e0), (500.,0.4532000e1), (500.,0.8380250e3),
        (500.,0.1084564e4), (647.,0.3580000e3), (900.,0.2410000e0),
        (900.,0.5261500e2), (900.,0.8707690e3)
    ]
    argfmt = '({0:3g},{1:10.4e})'
    refs_scaled = [
        [0.992418352e-1,0.200022515e2,0.700004704e3,0.999679423e-1,
            0.999938125,0.100003858e2,0.700000405e3,0.220384756e2,
            0.100062559,0.200000690e2,0.700000006e3],
        [0.413018112e1,0.406798347e1,0.346135580e1,0.150817541e1,
            0.166991025e1,0.322106219e1,0.307437693e1,0.618315728e1,
            0.175890657e1,0.193510526e1,0.266422350e1],
        [0.150151914e4,0.153492501e4,0.244357992e4,0.548314253e3,
            0.535739001e3,0.127128441e4,0.241200877e4,0.252145078e3,
            0.724027147e3,0.698445674e3,0.201933608e4],
        [0.393062643,0.387405401,0.132609616,0.794488271e1,
            0.682502725e1,0.256690919e1,0.203237509e1,0.432092307e1,
            0.916653194e1,0.659070225e1,0.417223802e1]
    ]
    scales = [1e6,1e3,1.,1e3]
    refs = [[r*scale for r in ref]
        for (ref,scale) in zip(refs_scaled,scales)]
    header = 'Fluid water functions'
    test = Tester(funs,fargs,refs,fnames,argfmt,header=header)
    test.run()
    if printresult:
        test.printresults(chktol=chktol)
    return test


## Main function: Check tables
if __name__ == '__main__':
    test = chkiapws95table7()

