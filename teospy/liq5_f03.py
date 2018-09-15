"""Liquid water Gibbs energy function from Feistel 2003.

This module provides the Gibbs free energy of liquid water and related
properties using the Feistel (2003) formulation. This module can also be
called as a function,

    python liq5_f03.py

which will compare the results of this module to reference values from
IAPWS 2009, table 6.

:Examples:

>>> liq_g(0,0,300.,1e5)
-5.26505056073e3
>>> liq_g(0,1,300.,1e5)
1.00345554745e-3
>>> liq_g(2,0,300.,1e5)
-13.9354762020
>>> cp(300.,1e5)
4180.64286060
>>> density(300.,1e5)
996.556352243
>>> soundspeed(300.,1e5)
1501.52808421

:Functions:

* :func:`liq_g`: Liquid water Gibbs energy with derivatives.
* :func:`cp`: Liquid water isobaric heat capacity.
* :func:`density`: Liquid water density.
* :func:`expansion`: Liquid water thermal expansion coefficient.
* :func:`kappa_t`: Liquid water isothermal compressibility
* :func:`soundspeed`: Liquid water sound speed.
* :func:`enthalpy`: Liquid water enthalpy.
* :func:`entropy`: Liquid water entropy.
* :func:`helmholtzenergy`: Liquid water Helmholtz free energy.
* :func:`internalenergy`: Liquid water internal energy.
* :func:`chkiapws09table6`: Check module against IAPWS 2009, table 6.

"""

__all__ = ['liq_g',
    'cp','density','expansion','kappa_t','soundspeed',
    'entropy','enthalpy','helmholtzenergy','internalenergy',
    'chkiapws09table6']

import constants0

_CHKTOL = constants0.CHKTOL
_PATM = constants0.PATM
_TCELS = constants0.TCELS
_C_F03 = ((40.,1e8),
    (0,0,101.34274313967416), (1,0,5.9057834791235253),
    (2,0,-12357.785933039),   (3,0,736.741204151612),
    (4,0,-148.185936433658),  (5,0,58.0259125842571),
    (6,0,-18.9843846514172),  (7,0,3.05081646487967),
    (0,1,100015.695367145),   (1,1,-270.983805184062),
    (2,1,1455.0364540468),    (3,1,-672.50778314507),
    (4,1,397.968445406972),   (5,1,-194.618310617595),
    (6,1,63.5113936641785),   (7,1,-9.63108119393062),
    (0,2,-2544.5765420363),   (1,2,776.153611613101),
    (2,2,-756.558385769359),  (3,2,499.360390819152),
    (4,2,-301.815380621876),  (5,2,120.520654902025),
    (6,2,-22.2897317140459),  (0,3,284.517778446287),
    (1,3,-196.51255088122),   (2,3,273.479662323528),
    (3,3,-239.545330654412),  (4,3,152.196371733841),
    (5,3,-55.2723052340152),  (6,3,8.17060541818112),
    (0,4,-33.3146754253611),  (1,4,28.9796526294175),
    (2,4,-55.5604063817218),  (3,4,48.8012518593872),
    (4,4,-26.3748377232802),  (5,4,6.48190668077221),
    (0,5,4.20263108803084),   (1,5,-2.13290083518327),
    (2,5,4.34420671917197),   (3,5,-1.66307106208905),
    (0,6,-0.546428511471039))


## Gibbs energy function
def liq_g(drvt,drvp,temp,pres):
    """Calculate liquid water Gibbs energy using F03.
    
    Calculate the specific Gibbs free energy of liquid water or its
    derivatives with respect to temperature and pressure using the
    Feistel (2003) polynomial formulation.
    
    :arg int drvt: Number of temperature derivatives.
    :arg int drvp: Number of pressure derivatives.
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :returns: Gibbs free energy in units of
        (J/kg) / K^drvt / Pa^drvp.
    :raises ValueError: If drvt or drvp are negative.
    
    :Examples:
    
    >>> liq_g(0,0,300.,1e5)
    -5.26505056073e3
    >>> liq_g(1,0,300.,1e5)
    -393.062597709
    >>> liq_g(0,1,300.,1e5)
    1.00345554745e-3
    >>> liq_g(2,0,300.,1e5)
    -13.9354762020
    >>> liq_g(1,1,300.,1e5)
    2.75754520492e-7
    >>> liq_g(0,2,300.,1e5)
    -4.52067557155e-13
    """
    if drvt < 0 or drvp < 0:
        errmsg = 'Derivatives {0} cannot be negative'.format((drvt,drvp))
        raise ValueError(errmsg)
    
    TRED, PRED = _C_F03[0]
    y = (temp - _TCELS)/TRED
    z = (pres - _PATM)/PRED
    g = 0.
    for (j,k,c) in _C_F03[1:]:
        if y==0:
            if j==drvt:
                pwrt = 1.
            else:
                pwrt = 0.
        else:
            pwrt = y**(j-drvt)
        for l in range(drvt):
            pwrt *= j-l
        if z==0:
            if k==drvp:
                pwrp = 1.
            else:
                pwrp = 0.
        else:
            pwrp = z**(k-drvp)
        for l in range(drvp):
            pwrp *= k-l
        g += c * pwrt * pwrp
    g /= TRED**drvt * PRED**drvp
    return g


## Thermodynamic properties
def cp(temp,pres):
    """Calculate liquid water isobaric heat capacity using F03.
    
    Calculate the isobaric (constant pressure) heat capacity of liquid
    water using the Feistel (2003) polynomial formulation.
    
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :returns: Heat capacity in J/kg/K.
    
    :Examples:
    
    >>> cp(300.,1e5)
    4180.64286060
    """
    g_tt = liq_g(2,0,temp,pres)
    cp = -temp * g_tt
    return cp

def density(temp,pres):
    """Calculate liquid water density using F03.
    
    Calculate the density of liquid water using the Feistel (2003)
    polynomial formulation.
    
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :returns: Density in kg/m3.
    
    :Examples:
    
    >>> density(300.,1e5)
    996.556352243
    """
    g_p = liq_g(0,1,temp,pres)
    dliq = g_p**(-1)
    return dliq

def expansion(temp,pres):
    """Calculate liquid water thermal expansion coefficient using F03.
    
    Calculate the thermal expansion coefficient of liquid water using
    the Feistel (2003) polynomial formulation.
    
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :returns: Expansion coefficient in 1/K.
    
    :Examples:
    
    >>> expansion(300.,1e5)
    2.74804919056e-4
    """
    g_p = liq_g(0,1,temp,pres)
    g_tp = liq_g(1,1,temp,pres)
    alpha = g_tp / g_p
    return alpha

def kappa_t(temp,pres):
    """Calculate liquid water isothermal compressibility using F03.
    
    Calculate the isothermal compressibility of liquid water using the
    Feistel (2003) polynomial formulation.
    
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :returns: Compressibility in 1/Pa.
    
    :Examples:
    
    >>> kappa_t(300.,1e5)
    4.50510795725e-10
    """
    g_p = liq_g(0,1,temp,pres)
    g_pp = liq_g(0,2,temp,pres)
    kappa = -g_pp/g_p
    return kappa

def soundspeed(temp,pres):
    """Calculate liquid water sound speed using F03.
    
    Calculate the speed of sound in liquid water using the Feistel
    (2003) polynomial formulation.
    
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :returns: Sound speed in m/s.
    
    :Examples:
    
    >>> soundspeed(300.,1e5)
    1501.52808421
    """
    g_p = liq_g(0,1,temp,pres)
    g_tt = liq_g(2,0,temp,pres)
    g_tp = liq_g(1,1,temp,pres)
    g_pp = liq_g(0,2,temp,pres)
    csqinv = (g_tp**2/g_tt - g_pp) / g_p**2
    c = csqinv**(-.5)
    return c


## Functions not in the original Fortran
def enthalpy(temp,pres):
    """Calculate liquid water enthalpy using F03.
    
    Calculate the specific enthalpy of liquid water using the Feistel
    (2003) polynomial formulation.
    
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :returns: Enthalpy in J/kg.
    """
    g = liq_g(0,0,temp,pres)
    g_t = liq_g(1,0,temp,pres)
    h = g - temp*g_t
    return h

def entropy(temp,pres):
    """Calculate liquid water entropy using F03.
    
    Calculate the specific entropy of liquid water using the Feistel
    (2003) polynomial formulation.
    
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :returns: Entropy in J/kg/K.
    """
    g_t = liq_g(1,0,temp,pres)
    s = -g_t
    return s

def helmholtzenergy(temp,pres):
    """Calculate liquid water Helmholtz energy using F03.
    
    Calculate the specific Helmholtz free energy of liquid water using
    the Feistel (2003) polynomial formulation.
    
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :returns: Helmholtz energy in J/kg.
    """
    g = liq_g(0,0,temp,pres)
    g_p = liq_g(0,1,temp,pres)
    f = g - pres*g_p
    return f

def internalenergy(temp,pres):
    """Calculate liquid water internal energy using F03.
    
    Calculate the specific internal energy of liquid water using the
    Feistel (2003) polynomial formulation.
    
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :returns: Internal energy in J/kg.
    """
    g = liq_g(0,0,temp,pres)
    g_t = liq_g(1,0,temp,pres)
    g_p = liq_g(0,1,temp,pres)
    u = g - temp*g_t - pres*g_p
    return u


## Check tables
def chkiapws09table6(printresult=True,chktol=_CHKTOL):
    """Check accuracy against IAPWS 2009 table 6.
    
    Evaluate the functions in this module and compare to reference
    values from IAPWS 2009, table 6.
    
    :arg bool printresult: If True (default) and any results are outside
        of the given tolerance, then the function name, reference value,
        result value, and relative error are printed.
    :arg float chktol: Tolerance to use when choosing to print results
        (default _CHKTOL).
    :returns: :class:`~tester.Tester` instances containing the
        functions, arguments, reference values, results, and relative
        errors from the tests. The first instance tests liq_g and the
        second instance tests thermodynamic properties.
    """
    from tester import Tester
    args1 = [(273.15,101325.), (273.15,1e8), (313.15,101325.)]
    DERS2 = ((0,0),(1,0),(0,1),(2,0),(1,1),(0,2))
    
    funs = liq_g
    fargs = [(der+args) for args in args1 for der in DERS2]
    refs = [
            0.101342743e3,0.147644587,0.100015695e-2,-0.154472324e2,
                -0.677459513e-7,-0.508915308e-12,
            0.977303868e5,0.851506346e1,0.956683354e-3,-0.142970174e2,
                0.199088060e-6,-0.371527164e-12,
            -0.116198898e5,-0.572365181e3,0.100784471e-2,-0.133463968e2,
                0.388499694e-6,-0.445841077e-12
    ]
    fnames = 'liq_g'
    argfmt = '({0:1g},{1:1g},{2:6.2f},{3:6g})'
    header = 'F03 liq_g derivatives'
    testder = Tester(funs,fargs,refs,fnames,argfmt,header=header)
    
    funs = [cp,density,soundspeed,
        enthalpy,entropy,helmholtzenergy,internalenergy]
    fargs = args1
    refs = [
        [0.421941153e4,0.390523030e4,0.417942416e4],
        [0.999843071e3,0.104527793e4,0.992216354e3],
        [0.140240099e4,0.157543089e4,0.152891242e4],
        [0.610136242e2,0.954044973e5,0.167616267e6],
        [-0.147644587,-0.851506346e1,0.572365181e3],
        [0.183980891e-2,0.206205140e4,-0.117220097e5],
        [-0.403272791e2,-0.263838183e3,0.167514147e6]
    ]
    fnames = ['cp','density','soundspeed',
        'enthalpy','entropy','helmholtzenergy','internalenergy']
    argfmt = '({0:6.2f},{1:6g})'
    header = 'F03 thermodynamic properties'
    testprop = Tester(funs,fargs,refs,fnames,argfmt,header=header)
    
    # Run Tester instances and print results
    testder.run()
    testprop.run()
    if printresult:
        testder.printresults(chktol=chktol)
        testprop.printresults(chktol=chktol)
    return (testder, testprop)


## Main function: Check tables
if __name__ == '__main__':
    testder, testprop = chkiapws09table6()

