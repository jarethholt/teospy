"""Dry air Helmholtz potential and air-water virial coefficients.

This module implements the Helmholtz free energy of dry air and its
derivatives with respect to temperature and density. This module also
includes the virial coefficients for dry air-water vapour mixtures.

:Examples:

>>> dry_f(0,0,300.,1e-3)
-696239.965190
>>> air_baw(0,300.)
-2.95672747e-05
>>> air_caaw(0,300.)
8.01977741e-10
>>> air_caww(0,300.)
-1.15552784e-07

:Functions:

* :func:`dry_f`: Dry air Helmholtz free energy with derivatives.
* :func:`air_baw`: Air-water virial coefficient with derivatives.
* :func:`air_caaw`: Air-air-water virial coefficient with derivatives.
* :func:`air_caww`: Air-water-water virial coefficient with derivatives.

"""

__all__ = ['dry_f','air_baw','air_caaw','air_caww']

import numpy
import constants0

# Single constants
_RDRY = constants0.RDRY_L2000
_MDRY = constants0.MDRY
_TRED_DRY = 132.6312  # Reducing temperature, K
_DRED_DRY = 10447.7 * _MDRY  # Reducing density, kg/m3
_TRED_AW = 100.0  # Reducing temperature (K) for air-water coefficient
_TRED_AAW = 1.0  # Reducing temperature (K) for air-air-water coefficient
_TRED_AWW = 1.0  # Reducing temperature (K) for air-water-water coefficient

# Constants used in empirical equations
_C_DRYF0 = (
    (9.7450251743948, 10.0986147428912),
    (6.057194e-8, -2.10274769e-5, -1.58860716e-4),
    (1.5, -1.9536342e-4),
    2.490888032,
    (
        (1., 1., -1., -25.36365,  0.791309509),
        (1., 1., -1., -16.90741,  0.212236768),
        (2., 3.,  1.,  87.31279, -0.197938904)
    )
)
_C_DRYFR = (
    (
        (1, 0.  ,  0.118160747229   ),
        (1, 0.33,  0.713116392079   ),
        (1, 1.01, -1.61824192067    ),
        (2, 0.  ,  0.0714140178971  ),
        (3, 0.  , -0.0865421396646  ),
        (3, 0.15,  0.134211176704   ),
        (4, 0.  ,  0.0112626704218  ),
        (4, 0.2 , -0.0420533228842  ),
        (4, 0.35,  0.0349008431982  ),
        (6, 1.35,  0.000164957183186)
    ),
    (
        ( 1, 1,  1.6 , -0.101365037912   ),
        ( 3, 1,  0.8 , -0.17381369097    ),
        ( 5, 1,  0.95, -0.0472103183731  ),
        ( 6, 1,  1.25, -0.0122523554253  ),
        ( 1, 2,  3.6 , -0.146629609713   ),
        ( 3, 2,  6.  , -0.0316055879821  ),
        (11, 2,  3.25,  0.000233594806142),
        ( 1, 3,  3.5 ,  0.0148287891978  ),
        ( 3, 3, 15.  , -0.00938782884667 )
    )
)
_C_AW = (
    (-0.237,   66.5687),
    (-1.048, -238.834 ),
    (-3.183, -176.755 )
)
_C_AAW = (4.82737e2,1.05678e5,-6.56394e7,2.94442e10,-3.19317e12)
_C_AWW = (-1.0728876e1,3.47802e3,-3.83383e5,3.3406e7)


## Dry air auxiliary functions
def _alpha_ideal(tau):
    """Calculate dry air potential ideal term.
    
    Calculate the temperature-dependent terms of the ideal gas component
    of the Helmholtz potential (scaled Helmholtz free energy) for dry
    air.
    
    :arg float tau: Reduced temperature _TRED_DRY/temp(K).
    :returns: Helmholtz potential, unitless.
    """
    alpha = 0.0
    n4, n5 = _C_DRYF0[0]
    alpha += n4 + n5*tau
    for (k,n) in enumerate(_C_DRYF0[1]):
        alpha += n * tau**(k-3)
    k, n = _C_DRYF0[2]
    alpha += n * tau**k
    alpha += _C_DRYF0[3] * numpy.log(tau)
    for (a1,a2,b,c,n) in _C_DRYF0[4]:
        alpha += n * numpy.log(a1/a2 + b*numpy.exp(c*tau))
    return alpha

def _alpha_ideal_t(tau):
    """Calculate dry air potential ideal term T-derivative.
    
    Calculate the derivative of the ideal gas component of the dry air
    Helmholtz potential (scaled Helmholtz free energy) with respect to
    reduced temperature.
    
    :arg float tau: Reduced temperature _TRED_DRY/temp(K).
    :returns: Helmholtz potential derivative, unitless.
    """
    alpha = 0.0
    n4, n5 = _C_DRYF0[0]
    alpha += n5
    for (k,n) in enumerate(_C_DRYF0[1]):
        alpha += n * (k-3) * tau**(k-4)
    k, n = _C_DRYF0[2]
    alpha += n * k*tau**(k-1)
    alpha += _C_DRYF0[3] / tau
    for (a1,a2,b,c,n) in _C_DRYF0[4]:
        eterm = numpy.exp(c*tau)
        alpha += n * b*c*eterm / (a1/a2 + b*eterm)
    return alpha

def _alpha_ideal_tt(tau):
    """Calculate dry air potential ideal term TT-derivative.
    
    Calculate the second derivative of the ideal gas component of the
    dry air Helmholtz potential (scaled Helmholtz free energy) with
    respect to reduced temperature.
    
    :arg float tau: Reduced temperature _TRED_DRY/temp(K).
    :returns: Helmholtz potential second derivative, unitless.
    """
    alpha = 0.0
    for (k,n) in enumerate(_C_DRYF0[1]):
        alpha += n * (k-3)*(k-4) * tau**(k-5)
    k, n = _C_DRYF0[2]
    alpha += n * k*(k-1)*tau**(k-2)
    alpha += -_C_DRYF0[3] / tau**2
    for (a1,a2,b,c,n) in _C_DRYF0[4]:
        eterm = numpy.exp(c*tau)
        denom = a1/a2 + b*eterm
        alpha += n * a1/a2 * b * c**2 * eterm / denom**2
    return alpha

def _alpha_res(drvt,drvd,tau,dta):
    """Calculate dry air potential residual term.
    
    Calculate the residual (non-ideal gas) component of the Helmholtz
    potential (scaled Helmholtz free energy) or its derivatives with
    respect to reduced temperature and density. Derivatives up to
    second-order are available.
    
    :arg int drvt: Number of reduced temperature derivatives.
    :arg int drvd: Number of reduced density derivatives.
    :arg float tau: Reduced temperature _TRED_DRY/temp(K).
    :arg float dta: Reduced density ddry(kg/m3)/_DRED_DRY.
    :returns: Helmholtz potential or derivative, unitless.
    :raises ValueError: If drvt<0, drvd<0, or drvt+drvd>2.
    """
    if (drvt < 0 or drvd < 0 or drvt+drvd > 2):
        errmsg = 'Derivative {0} not recognized'.format((drvt,drvd))
        raise ValueError(errmsg)
    alpha = 0.0
    
    # First part: dual power series
    for (j,k,n) in _C_DRYFR[0]:
        if drvt == 0:
            a_tau = tau**k
        elif drvt == 1:
            a_tau = k * tau**(k-1)
        elif drvt == 2:
            a_tau = k*(k-1) * tau**(k-2)
        
        if drvd == 0:
            a_dta = dta**j
        elif drvd == 1:
            a_dta = j * dta**(j-1)
        elif drvd == 2:
            a_dta = j*(j-1) * dta**(j-2)
        alpha += n * a_tau * a_dta
    
    # Second part: power series with exponential correction
    for (j,l,k,n) in _C_DRYFR[1]:
        if drvt == 0:
            a_tau = tau**k
        elif drvt == 1:
            a_tau = k * tau**(k-1)
        elif drvt == 2:
            a_tau = k*(k-1) * tau**(k-2)
        
        dtal = dta**l
        eterm = numpy.exp(-dtal)
        if drvd == 0:
            a_dta = dta**j * eterm
        elif drvd == 1:
            a_dta = dta**(j-1)*eterm * (j-l*dtal)
        elif drvd == 2:
            a_dta = dta**(j-2)*eterm * ((j-1-l*dtal)*(j-l*dtal) - l**2*dtal)
        alpha += n * a_tau * a_dta
    return alpha


### Public functions
def dry_f(drvt,drvd,temp,ddry,chkbnd=False):
    """Calculate dry air Helmholtz free energy.
    
    Calculate the specific Helmholtz free energy of dry air or its
    derivatives with respect to temperature and density. Derivatives up
    to second order are available.
    
    :arg int drvt: Number of temperature derivatives.
    :arg int drvd: Number of density derivatives.
    :arg float temp: Temperature in K.
    :arg float ddry: Dry air density in kg/m3.
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :returns: Helmholtz free energy in units of
        (J/kg) / K^drvt / (kg/m3)^drvd.
    :raises ValueError: If either temp or ddry are nonpositive.
    :raises RuntimeWarning: If temp or ddry are outside the recommended
        bounds and chkbnd is True.
    :raises ValueError: If drvt<0, drvd<0, or drvt+drvd>2.
    
    :Examples:
    
    >>> dry_f(0,0,300.,1e-3)
    -696239.965190
    >>> dry_f(1,0,300.,1e-3)
    -2124.55145456
    >>> dry_f(0,1,300.,1e-3)
    8.61147149596e+07
    >>> dry_f(2,0,300.,1e-3)
    -2.39242390806
    >>> dry_f(1,1,300.,1e-3)
    287049.624545
    >>> dry_f(0,2,300.,1e-3)
    -8.61147380361e+10
    """
    constants0.chkdrybnds(temp,ddry,chkbnd=chkbnd)
    tau = _TRED_DRY / temp
    dta = ddry / _DRED_DRY
    
    # Run through each derivative case
    if (drvt,drvd) == (0,0):
        alpha = numpy.log(dta) + _alpha_ideal(tau)
        alpha += _alpha_res(0,0,tau,dta)
        f = _RDRY * temp * alpha
    elif (drvt,drvd) == (0,1):
        alpha_d = 1./dta + _alpha_res(0,1,tau,dta)
        f = _RDRY * temp * alpha_d / _DRED_DRY
    elif (drvt,drvd) == (0,2):
        alpha_dd = -1./dta**2 + _alpha_res(0,2,tau,dta)
        f = _RDRY * temp * alpha_dd / _DRED_DRY**2
    elif (drvt,drvd) == (1,0):
        alpha = numpy.log(dta) + _alpha_ideal(tau)
        alpha += _alpha_res(0,0,tau,dta)
        alpha_t = _alpha_ideal_t(tau)
        alpha_t += _alpha_res(1,0,tau,dta)
        f = _RDRY * (alpha - tau*alpha_t)
    elif (drvt,drvd) == (1,1):
        alpha_d = 1./dta + _alpha_res(0,1,tau,dta)
        alpha_td = _alpha_res(1,1,tau,dta)
        f = _RDRY * (alpha_d - tau*alpha_td) / _DRED_DRY
    elif (drvt,drvd) == (2,0):
        alpha_tt = _alpha_ideal_tt(tau)
        alpha_tt += _alpha_res(2,0,tau,dta)
        f = _RDRY * tau**2 * alpha_tt / temp
    else:
        errmsg = 'Derivatives {0} not recognized'.format((drvt,drvd))
        raise ValueError(errmsg)
    return f

def air_baw(drvt,temp):
    """Calculate air-water virial coefficient.
    
    Calculate the first dry air-water vapour virial coefficient or its
    derivative with respect to temperature. Derivatives up to second
    order are available.
    
    :arg int drvt: Number of temperature derivatives.
    :arg float temp: Temperature in K.
    :returns: Air-water virial coefficient, in units of
        (m3/mol) / K^drvt.
    :raises ValueError: If drvt<0 or drvt>2.
    
    :Examples:
    
    >>> air_baw(0,300.)
    -2.95672747e-05
    >>> air_baw(1,300.)
    2.80097360e-07
    >>> air_baw(2,300.)
    -2.42599241e-09
    """
    baw = 0.0
    tau = temp / _TRED_AW
    for (k,n) in _C_AW:
        if drvt == 0: baw += n * tau**k
        elif drvt == 1: baw += n * k * tau**(k-1) / _TRED_AW
        elif drvt == 2: baw += n * k*(k-1) * tau**(k-2) / _TRED_AW**2
        else:
            errmsg = 'Derivative {0} not recognized'.format(drvt)
            raise ValueError(errmsg)
    # Convert from cm3/mol to m3/mol
    baw *= 1e-6
    return baw

def air_caaw(drvt,temp):
    """Calculate air-air-water virial coefficient.
    
    Calculate the second dry air-dry air-water vapour virial coefficient
    or its derivative with respect to temperature. Derivatives up to
    second order are available.
    
    :arg int drvt: Number of temperature derivatives.
    :arg float temp: Temperature in K.
    :returns: Air-air-water virial coefficient, in units of
        (m3/mol)^2 / K^drvt.
    :raises ValueError: If drvt<0 or drvt>2.
    
    :Examples:
    
    >>> air_caaw(0,300.)
    8.01977741e-10 
    >>> air_caaw(1,300.)
    -1.96103457e-12
    >>> air_caaw(2,300.)
    1.70055638e-14
    """
    caaw = 0.0
    tau = temp / _TRED_AAW
    for (k,n) in enumerate(_C_AAW):
        if drvt == 0: caaw += n * tau**(-k)
        elif drvt == 1: caaw += n * (-k) * tau**(-k-1)
        elif drvt == 2: caaw += n * k*(k+1) * tau**(-k-2)
        else:
            errmsg = 'Derivative {0} not recognized'.format(drvt)
            raise ValueError(errmsg)
    # Convert from cm6/mol2 to m6/mol2
    caaw *= 1e-12
    return caaw

def air_caww(drvt,temp):
    """Calculate air-water-water virial coefficient.
    
    Calculate the second dry air-water vapour-water vapour virial
    coefficient or its derivative with respect to temperature.
    Derivatives up to second order are available.
    
    :arg int drvt: Number of temperature derivatives.
    :arg float temp: Temperature in K.
    :returns: Air-water-water virial coefficient, in units of
        (m3/mol)^2 / K^drvt.
    :raises ValueError: If drvt<0 or drvt>2.
    
    :Examples:
    
    >>> air_caww(0,300.)
    -1.15552784e-07
    >>> air_caww(1,300.)
    2.61363278e-09 
    >>> air_caww(2,300.)
    -7.51334582e-11
    """
    caww = 0.0
    tau = temp / _TRED_AWW
    if not (drvt in (0,1,2)):
        errmsg = 'Derivative {0} not recognized'.format(drvt)
        raise ValueError(errmsg)
    
    # Calculate leading exponential term
    earg = 0.0
    for (k,n) in enumerate(_C_AWW):
        earg += n * tau**(-k)
    caww = -1e-6 * numpy.exp(earg)
    
    # Multiply by derivatives of the exponent
    if drvt == 0:
        der = 1.0
    elif drvt == 1:
        der = 0.0
        for (k,n) in enumerate(_C_AWW):
            der += n * (-k) * tau**(-k-1)
    elif drvt == 2:
        der1, der2 = 0.0, 0.0
        for (k,n) in enumerate(_C_AWW):
            der1 += n * (-k) * tau**(-k-1)
            der2 += n * k*(k+1) * tau**(-k-2)
        der = der1**2 + der2
    caww *= der
    return caww


