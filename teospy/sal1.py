"""Seawater salt Gibbs free energy coefficients.

This module implements the coefficients of the power expansion of the
Gibbs free energy of salt in seawater and its derivatives with respect
to salinity, temperature, and pressure. The number of terms in the
expansion is provided by this module as `NSALTERMS`.

There are two sets of coefficients used in this module, `_GSCOEFFS` and
`_GSCOEFFS_EXT`. The `_EXT` version is an extension to the high-
temperature, high-salinity regime. Any function involving sea salt
supports the `useext` keyword. If True, the coefficients for the
extension are used. The default is False.

:Functions:

* :func:`sal_g_term`: Calculate a term in the Gibbs free energy of salt.

"""

__all__ = ['NSALTERMS','sal_g_term']

import constants0

# Single constants
_TCELS = constants0.TCELS
_PATM = constants0.PATM
_TRED = 40.
_PRED = 1e8

# Constants used in empirical equations
_GSCOEFFS = (
	(1,0,0, 7.23191662157061e+04), (1,1,0, 1.05903959312767e+04),
	(2,0,0, 2.67686938386265e+05), (2,0,1,-8.23738604560763e+04),
	(2,0,2, 9.57470498700625e+03), (2,0,3,-2.40198441433317e+03),
	(2,0,4, 3.94161790150021e+02), (2,0,5,-6.53120647714656e+01),
	(2,1,0, 3.82214384164908e+04), (2,1,1, 1.81423642207762e+04),
	(2,1,2,-8.55856527083052e+03), (2,1,3, 3.10256186902438e+03),
	(2,1,4,-7.87709720728709e+02), (2,1,5, 1.75337907450416e+02),
	(2,2,0, 2.18975275976525e+04), (2,2,1,-2.14181120172472e+04),
	(2,2,2, 8.39564917274930e+03), (2,2,3,-4.43694182294925e+03),
	(2,2,4, 1.09991432832012e+03), (2,2,5,-1.97071112050554e+02),
	(2,3,0,-5.60526002829731e+03), (2,3,1, 1.72746629142640e+04),
	(2,3,2,-5.09820084718488e+03), (2,3,3, 2.82571797474195e+03),
	(2,3,4,-2.76901128324072e+02), (2,4,0, 2.27492387832447e+03),
	(2,4,1,-7.40828530946612e+03), (2,4,2, 1.85938572788234e+03),
	(2,4,3,-9.07901154460196e+02), (2,5,0,-5.38966643478572e+02),
	(2,6,0, 5.30043046990093e+01), (3,0,0,-3.01880568417681e+05),
	(3,0,1, 2.47571333745467e+04), (3,0,2,-6.49079695264810e+03),
	(3,0,3, 8.44575339642748e+03), (3,0,4,-4.23565447618004e+02),
	(3,1,0,-6.12422533100840e+04), (3,1,1,-2.17574304584572e+04),
	(3,1,2, 1.03259263132272e+04), (3,1,3,-3.65946860988548e+03),
	(3,2,0,-5.34545474801933e+03), (3,2,1, 4.75455655384292e+04),
	(3,2,2,-6.72633342258010e+03), (3,2,3, 3.18244469131541e+03),
	(3,3,0,-1.24403255093626e+03), (3,3,1,-5.71353886748868e+04),
	(3,4,0, 1.08680464781702e+02), (3,4,1, 2.91144750376762e+04),
	(4,0,0, 1.25426900174329e+06), (4,0,1,-3.39242567121826e+04),
	(4,0,2,-2.52732136578246e+03), (4,0,3,-1.86830816013266e+04),
	(4,1,0, 3.36714093681915e+05), (4,1,1,-1.40350478031969e+04),
	(4,2,0,-4.24470261644411e+04), (4,3,0, 3.05652943762545e+04),
	(4,4,0,-1.06120320689368e+04), (4,5,0, 1.54599190546279e+03),
	(5,0,0,-3.37157703069480e+06), (5,0,1, 1.11272425476394e+05),
	(5,1,0,-6.05426088990805e+05), (6,0,0, 5.77111887317487e+06),
	(6,1,0, 5.66282039731641e+05), (7,0,0,-3.73403338866189e+06) )
NSALTERMS = max(g[0] for g in _GSCOEFFS)
_GSCOEFFS_EXT = (
	(1,0,0, 7.23191662157061e+04), (1,1,0, 1.05903959312767e+04),
	(2,0,0, 2.67686938386265e+05), (2,0,1,-8.23738604560763e+04),
	(2,0,2, 9.57470498700625e+03), (2,0,3,-2.40198441433317e+03),
	(2,0,4, 3.94161790150021e+02), (2,0,5,-6.53120647714656e+01),
	(2,1,0, 3.82214384164908e+04), (2,1,1, 1.72774406106906e+04),
	(2,1,2,-8.55856527083052e+03), (2,1,3, 3.10256186902438e+03),
	(2,1,4,-7.87709720728709e+02), (2,1,5, 1.75337907450416e+02),
	(2,2,0, 2.18975275976525e+04), (2,2,1,-2.08505851518528e+04),
	(2,2,2, 8.39564917274930e+03), (2,2,3,-4.43694182294925e+03),
	(2,2,4, 1.09991432832012e+03), (2,2,5,-1.97071112050554e+02),
	(2,3,0,-5.60526002829731e+03), (2,3,1, 1.58794275914092e+04),
	(2,3,2,-5.09820084718488e+03), (2,3,3, 2.82571797474195e+03),
	(2,3,4,-2.76901128324072e+02), (2,4,0, 2.27492387832447e+03),
	(2,4,1,-7.40828530946612e+03), (2,4,2, 1.85938572788234e+03),
	(2,4,3,-9.07901154460196e+02), (2,5,0,-5.38966643478572e+02),
	(2,5,1, 1.60565070211913e+03), (2,6,0, 5.30043046990093e+01),
	(2,6,1,-1.07678765811425e+02), (3,0,0,-3.01880568417681e+05),
	(3,0,1, 2.47571333745467e+04), (3,0,2,-6.49079695264810e+03),
	(3,0,3, 8.44575339642748e+03), (3,0,4,-4.23565447618004e+02),
	(3,1,0,-6.12422533100840e+04), (3,1,1,-1.01499012805071e+04),
	(3,1,2, 1.03259263132272e+04), (3,1,3,-3.65946860988548e+03),
	(3,2,0,-5.34545474801933e+03), (3,2,1, 4.75455655384292e+04),
	(3,2,2,-6.72633342258010e+03), (3,2,3, 3.18244469131541e+03),
	(3,3,0,-1.24403255093626e+03), (3,3,1,-5.89168201312478e+04),
	(3,4,0, 1.08680464781702e+02), (3,4,1, 2.91144750376762e+04),
	(3,5,1,-5.00834767678767e+03), (4,0,0, 1.25426900174329e+06),
	(4,0,1,-3.39242567121826e+04), (4,0,2,-2.52732136578246e+03),
	(4,0,3,-1.86830816013266e+04), (4,1,0, 3.36714093681915e+05),
	(4,1,1,-5.13995561277463e+04), (4,2,0,-4.24470261644411e+04),
	(4,2,1,-1.50372111563156e+04), (4,3,0, 3.05652943762545e+04),
	(4,3,1, 5.73557723369509e+04), (4,4,0,-1.06120320689368e+04),
	(4,4,1,-2.57973168671926e+04), (4,5,0, 1.54599190546279e+03),
	(4,6,1, 1.27516274815843e+03), (5,0,0,-3.37157703069480e+06),
	(5,0,1, 1.11272425476394e+05), (5,1,0,-6.05426088990805e+05),
	(6,0,0, 5.77111887317487e+06), (6,1,0, 5.66282039731641e+05),
	(7,0,0,-3.73403338866189e+06) )


## Auxiliary functions
def _poly_gyz(term,drvy,drvz,y,z,useext=False):
    """Evaluate a nondimensional polynomial from GSCOEFFS.
    
    Evaluate a bivariate polynomial or its derivatives whose
    coefficients are given in _GSCOEFFS. The full polynomial for a given
    term is
    
        g(i;y,z) = sum_{(i=term,j,k,c) in _GSCOEFFS} c * y^j * z^k.
    
    The high-temperature/high-salinity extension defined by
    _GSCOEFFS_EXT can also be used by setting useext to True.
    
    :arg int term: Number of the term in the salinity expansion to
        calculate, between 1 and NSALTERMS.
    :arg int drvy: Number of y-derivatives to take.
    :arg int drvz: Number of z-derivatives to take.
    :arg float y: Value of the first variable.
    :arg float z: Value of the second variable.
    :arg bool useext: If False (default) then the polynomial is
        calculated from _GSCOEFFS; if True, from _GSCOEFFS_EXT.
    :returns: Value of the polynomial or its derivatives.
    :raises ValueError: If drvy or drvz are negative.
    """
    if drvy < 0 or drvz < 0:
        errmsg = 'Derivatives {0} not recognized'.format((drvy,drvz))
        raise ValueError(errmsg)
    gi = 0.
    if useext:
        gsi = _GSCOEFFS_EXT
    else:
        gsi = _GSCOEFFS
    
    # Loop over powers of y and z
    for (i,j,k,c) in gsi:
        if (i != term) or (j < drvy) or (k < drvz):
            continue
        gijk = c * y**(j-drvy) * z**(k-drvz)
        for l in range(drvy):
            gijk *= (j-l)
        for l in range(drvz):
            gijk *= (k-l)
        gi += gijk
    return gi

def _poly_gtp(term,drvt,drvp,temp,pres,useext=False):
    """Evaluate a dimensional polynomial from GSCOEFFS.
    
    Evaluate a bivariate polynomial or its derivatives whose
    coefficients are given in _GSCOEFFS. The full polynomial for a given
    term is
    
        g(i;y,z) = sum_{(i=term,j,k,c) in _GSCOEFFS} c * y^j * z^k
        y, z = (temp-_TCELS)/_TRED, (pres-_PATM)/_PRED
    
    The high-temperature/high-salinity extension defined by
    _GSCOEFFS_EXT can also be used by setting useext to True.
    
    :arg int term: Number of the term in the salinity expansion to
        calculate, between 1 and NSALTERMS.
    :arg int drvt: Number of temperature derivatives to take.
    :arg int drvp: Number of pressure derivatives to take.
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg bool useext: If False (default) then the polynomial is
        calculated from _GSCOEFFS; if True, from _GSCOEFFS_EXT.
    :returns: Value of the polynomial in K^(-drvt) Pa^(-drvp).
    :raises ValueError: If drvt or drvp are negative.
    """
    y = (temp - _TCELS) / _TRED
    z = (pres - _PATM) / _PRED
    giyz = _poly_gyz(term,drvt,drvp,y,z,useext=useext)
    gi = giyz / (_TRED**drvt * _PRED**drvp)
    return gi


## Public functions
def sal_g_term(term,drvt,drvp,temp,pres,chkbnd=False,useext=False):
    """Calculate a salt Gibbs free energy term.
    
    Calculate a term in the expansion of the Gibbs free energy of salt
    in seawater in terms of salinity.
    
    :arg int term: Number of the term in the salinity expansion to
        calculate, between 1 and NSALTERMS.
    :arg int drvt: Number of temperature derivatives to take.
    :arg int drvp: Number of pressure derivatives to take.
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg bool useext: If False (default) then the polynomial is
        calculated from _GSCOEFFS; if True, from _GSCOEFFS_EXT.
    :returns: Gibbs free energy in units of
        (J/kg) / K^drvt / Pa^drvp.
    :raises ValueError: If term is not between 1 and NSALTERMS.
    :raises ValueError: If temp or pres are nonpositive.
    :raises RuntimeWarning: If temp or pres are outside the recommended
        bounds and chkbnd is True.
    :raises ValueError: If drvt or drvp are negative.
    
    :Examples:
    
    >>> sal_g_term(1,0,0,300.,1e6)
    79427.9694846
    >>> sal_g_term(2,1,0,300.,1e6)
    1558.10730393
    >>> sal_g_term(4,0,2,300.,1e6)
    -6.06204383305e-13
    >>> sal_g_term(1,0,0,300.,1e6,useext=True)
    79427.9694846
    >>> sal_g_term(2,1,0,300.,1e6,useext=True)
    1558.10730392553
    >>> sal_g_term(4,0,2,300.,1e6,useext=True)
    -6.06204383305e-13
    """
    if term < 1 or term > NSALTERMS:
        errmsg = ('The value of term {0} is not between 1 and '
            '{1}').format(term,NSALTERMS)
        raise ValueError(errmsg)
    constants0.chksalbnds(0.,temp,pres,chkbnd=chkbnd)
    gi = _poly_gtp(term,drvt,drvp,temp,pres,useext=useext)
    return gi


