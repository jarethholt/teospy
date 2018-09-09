"""Fluid water Helmholtz function.

This module implements the Helmholtz free energy of fluid water (liquid
or vapour) and its derivatives with respect to temperature and density.

This module can also be called as a function::

    python flu1.py

which will compare results from this module to reference values from
IAPWS 1995, table 6.

:Examples:

>>> flu_f(0,0,300.,1000.)
-5351.74115204
>>> flu_f(0,2,300.,1000.)
2.24824656167

:Functions:

* flu_f: Fluid water Helmholtz free energy with derivatives.
* chkiapws95table6: Compare module results to reference values.

"""

__all__ = ['flu_f','chkiapws95table6']

import numpy
import constants0

# Constants
_CHKTOL = constants0.CHKTOL
_RWAT = constants0.RWAT_IAPWS95
_TCP = constants0.TCP
_DCP = constants0.DCP
_TLOW = 130.

# Constants used in empirical equations
_IDEALCOEFFS = (
    (-8.32044648374969,6.68321052759323,3.00632),
    (
        (1.28728967,0.012436), (3.53734222,0.97315 ), (7.74073708,1.2795  ),
        (9.24437796,0.96956 ), (27.5075105,0.24873 )
    )
)
_RESIDCOEFFS = (
    (
        (1,-0.5  , 0.012533547935523 ), (1, 0.875, 7.8957634722828   ),
        (1, 1.   ,-8.7803203303561   ), (2,  .5  , 0.31802509345418  ),
        (2,  .75 ,-0.26145533859358  ), (3,  .375,-7.8199751687981e-3),
        (4, 1.   , 8.8089493102134e-3)
    ),
    (
        (1, 1, 4,-0.66856572307965   ), (1, 1, 6, 0.20433810950965   ),
        (1, 1,12,-6.6212605039687e-5 ), (1, 2, 1,-0.19232721156002   ),
        (1, 2, 5,-0.25709043003438   ), (1, 3, 4, 0.16074868486251   ),
        (1, 4, 2,-0.040092828925807  ), (1, 4,13, 3.9343422603254e-7 ),
        (1, 5, 9,-7.5941377088144e-6 ), (1, 7, 3, 5.6250979351888e-4 ),
        (1, 9, 4,-1.5608652257135e-5 ), (1,10,11, 1.1537996422951e-9 ),
        (1,11, 4, 3.6582165144204e-7 ), (1,13,13,-1.3251180074668e-12),
        (1,15, 1,-6.2639586912454e-10), (2, 1, 7,-0.10793600908932   ),
        (2, 2, 1, 0.017611491008752  ), (2, 2, 9, 0.22132295167546   ),
        (2, 2,10,-0.40247669763528   ), (2, 3,10, 0.58083399985759   ),
        (2, 4, 3, 4.9969146990806e-3 ), (2, 4, 7,-0.031358700712549  ),
        (2, 4,10,-0.74315929710341   ), (2, 5,10, 0.4780732991548    ),
        (2, 6, 6, 0.020527940895948  ), (2, 6,10,-0.13636435110343   ),
        (2, 7,10, 0.014180634400617  ), (2, 9, 1, 8.3326504880713e-3 ),
        (2, 9, 2,-0.029052336009585  ), (2, 9, 3, 0.038615085574206  ),
        (2, 9, 4,-0.020393486513704  ), (2, 9, 8,-1.6554050063734e-3 ),
        (2,10, 6, 1.9955571979541e-3 ), (2,10, 9, 1.5870308324157e-4 ),
        (2,12, 8,-1.638856834253e-5  ), (3, 3,16, 0.043613615723811  ),
        (3, 4,22, 0.034994005463765  ), (3, 4,23,-0.076788197844621  ),
        (3, 5,23, 0.022446277332006  ), (4,14,10,-6.2689710414685e-5 ),
        (6, 3,50,-5.5711118565645e-10), (6, 6,44,-0.19905718354408   ),
        (6, 6,46, 0.31777497330738   ), (6, 6,50,-0.11841182425981   )
    ),
    (
        (3,0,20,1,150,1.21,-31.306260323435),
        (3,1,20,1,150,1.21, 31.546140237781),
        (3,4,20,1,250,1.25,-2521.3154341695)
    ),
    (
        (3.5,.85,.2,28,700,0.32,0.3,-0.14874640856724),
        (3.5,.95,.2,32,800,0.32,0.3, 0.31806110878444)
    )
)
_LOWCOEFFS = (0.278296458178592,-.5,-3.,4.5,.5,-4.5)


## Ideal gas auxiliary functions
def _phi0(tau,dta):
    """Calculate fluid water potential ideal term.
    
    Calculate the ideal gas component of the Helmholtz potential (scaled
    free energy) for fluid water.
    
    :arg float tau: Reduced temperature _TCP/temp(K).
    :arg float dta: Reduced density dflu(kg/m3)/_DCP.
    :returns: Helmholtz potential, unitless.
    """
    n0, n1, n2 = _IDEALCOEFFS[0]
    phi = numpy.log(dta) + n0 + n1*tau + n2*numpy.log(tau)
    for (g,n) in _IDEALCOEFFS[1]:
        eterm = numpy.exp(-g*tau)
        phi += n * numpy.log(1 - eterm)
    
    # Extension to low temperatures
    if tau > _TCP/_TLOW:
        COEFF, A0, A1, A2, A3, A4 = _LOWCOEFFS
        TAULOW = _TCP/_TLOW
        phi += COEFF * (A0 / tau
            + A1 * (tau+TAULOW)*numpy.log(tau/TAULOW)/TAULOW**2
            + A2 * tau/TAULOW**2
            + A3 * tau**2/TAULOW**3
            + A4 / TAULOW)
    return phi

def _phi0_t(tau,dta):
    """Calculate fluid water potential ideal term T-derivative.
    
    Calculate the derivative of the ideal gas component of the Helmholtz
    potential (scaled free energy) for fluid water with respect to
    reduced temperature.
    
    :arg float tau: Reduced temperature _TCP/temp(K).
    :arg float dta: Reduced density dflu(kg/m3)/_DCP.
    :returns: Helmholtz potential derivative, unitless.
    """
    n0, n1, n2 = _IDEALCOEFFS[0]
    phi = n1 + n2/tau
    for (g,n) in _IDEALCOEFFS[1]:
        eterm = numpy.exp(-g*tau)
        phi += n * g * eterm/(1-eterm)
    
    # Extension to low temperatures
    if tau > _TCP/_TLOW:
        COEFF, A0, A1, A2, A3, A4 = _LOWCOEFFS
        TAULOW = _TCP/_TLOW
        phi += COEFF * (-A0 / tau**2
            + A1 * (numpy.log(tau/TAULOW) + 1 + TAULOW/tau)/TAULOW**2
            + A2 / TAULOW**2
            + A3 * 2*tau/TAULOW**3)
    return phi

def _phi0_d(tau,dta):
    """Calculate fluid water potential ideal term D-derivative.
    
    Calculate the derivative of the ideal gas component of the Helmholtz
    potential (scaled free energy) for fluid water with respect to
    reduced density.
    
    :arg float tau: Reduced temperature _TCP/temp(K).
    :arg float dta: Reduced density dflu(kg/m3)/_DCP.
    :returns: Helmholtz potential derivative, unitless.
    """
    return 1./dta

def _phi0_tt(tau,dta):
    """Calculate fluid water potential ideal term TT-derivative.
    
    Calculate the second derivative of the ideal gas component of the
    Helmholtz potential (scaled free energy) for fluid water with
    respect to reduced temperature.
    
    :arg float tau: Reduced temperature _TCP/temp(K).
    :arg float dta: Reduced density dflu(kg/m3)/_DCP.
    :returns: Helmholtz potential derivative, unitless.
    """
    n0, n1, n2 = _IDEALCOEFFS[0]
    phi = -n2 / tau**2
    for (g,n) in _IDEALCOEFFS[1]:
        eterm = numpy.exp(-g*tau)
        denom = 1. - eterm
        phi += n * -g**2 * eterm / denom**2
    
    # Extension to low temperatures
    if tau > _TCP/_TLOW:
        COEFF, A0, A1, A2, A3, A4 = _LOWCOEFFS
        TAULOW = _TCP/_TLOW
        phi += COEFF * (A0 * 2. / tau**3
            + A1 * (1./tau + 1 - TAULOW/tau**2)/TAULOW**2
            + A3 * 2./TAULOW**3)
    return phi

def _phi0_td(tau,dta):
    """Calculate fluid water potential ideal term TD-derivative.
    
    Calculate the mixed derivative of the ideal gas component of the
    Helmholtz potential (scaled free energy) for fluid water with
    respect to reduced temperature and density.
    
    :arg float tau: Reduced temperature _TCP/temp(K).
    :arg float dta: Reduced density dflu(kg/m3)/_DCP.
    :returns: Helmholtz potential derivative, unitless.
    """
    return 0.

def _phi0_dd(tau,dta):
    """Calculate fluid water potential ideal term DD-derivative.
    
    Calculate the second derivative of the ideal gas component of the
    Helmholtz potential (scaled free energy) for fluid water with
    respect to reduced density.
    
    :arg float tau: Reduced temperature _TCP/temp(K).
    :arg float dta: Reduced density dflu(kg/m3)/_DCP.
    :returns: Helmholtz potential derivative, unitless.
    """
    return -dta**-2


### Auxiliary residual (non-ideal gas) functions
def _theta(i,tau,dta):
    """Calculate _theta term.
    
    Calculate the value of the i(th) _theta parameter (IAPWS-95,
    table 5). This parameter is used in calculating the residual (non-
    ideal gas) component of the Helmholtz free energy of fluid water.
    
    :arg int i: Which of the two _theta terms to compute, either 0 or 1.
    :arg float tau: Reduced temperature _TCP/temp(K).
    :arg float dta: Reduced density dflu(kg/m3)/_DCP.
    :returns: Theta value.
    """
    ai, bi, bb, cc, dd, aa, bet, n = _RESIDCOEFFS[3][i]
    adta1 = abs(dta-1)
    tht = 1. - tau + aa*adta1**(1./bet)
    return tht

def _theta_t(i,tau,dta):
    """Calculate _theta term T-derivative.
    
    Calculate the derivative of the i(th) _theta parameter (IAPWS-95,
    table 5) with respect to reduced temperature. This parameter is used
    in calculating the residual (non-ideal gas) component of the
    Helmholtz free energy of fluid water.
    
    :arg int i: Which of the two _theta terms to compute, either 0 or 1.
    :arg float tau: Reduced temperature _TCP/temp(K).
    :arg float dta: Reduced density dflu(kg/m3)/_DCP.
    :returns: Theta value derivative.
    """
    return -1.

def _theta_d(i,tau,dta):
    """Calculate _theta term D-derivative.
    
    Calculate the derivative of the i(th) _theta parameter (IAPWS-95,
    table 5) with respect to reduced density. This parameter is used in
    calculating the residual (non-ideal gas) component of the Helmholtz
    free energy of fluid water.
    
    :arg int i: Which of the two _theta terms to compute, either 0 or 1.
    :arg float tau: Reduced temperature _TCP/temp(K).
    :arg float dta: Reduced density dflu(kg/m3)/_DCP.
    :returns: Theta value derivative.
    """
    ai, bi, bb, cc, dd, aa, bet, n = _RESIDCOEFFS[3][i]
    betinv = 1./bet
    adta1 = abs(dta-1)
    tht_d = aa * betinv * (dta-1) * adta1**(betinv-2)
    return tht_d

def _theta_tt(i,tau,dta):
    """Calculate _theta term TT-derivative.
    
    Calculate the second derivative of the i(th) _theta parameter
    (IAPWS-95, table 5) with respect to reduced density. This parameter
    is used in calculating the residual (non-ideal gas) component of the
    Helmholtz free energy of fluid water.
    
    :arg int i: Which of the two _theta terms to compute, either 0 or 1.
    :arg float tau: Reduced temperature _TCP/temp(K).
    :arg float dta: Reduced density dflu(kg/m3)/_DCP.
    :returns: Theta value derivative.
    """
    return 0.

def _theta_td(i,tau,dta):
    """Calculate _theta term TD-derivative.
    
    Calculate the mixed derivative of the i(th) _theta parameter
    (IAPWS-95, table 5) with respect to reduced temperature and density.
    This parameter is used in calculating the residual (non-ideal gas)
    component of the Helmholtz free energy of fluid water.
    
    :arg int i: Which of the two _theta terms to compute, either 0 or 1.
    :arg float tau: Reduced temperature _TCP/temp(K).
    :arg float dta: Reduced density dflu(kg/m3)/_DCP.
    :returns: Theta value derivative.
    """
    return 0.

def _theta_dd(i,tau,dta):
    """Calculate _theta term DD-derivative.
    
    Calculate the second derivative of the i(th) _theta parameter
    (IAPWS-95, table 5) with respect to reduced density. This parameter
    is used in calculating the residual (non-ideal gas) component of the
    Helmholtz free energy of fluid water.
    
    :arg int i: Which of the two _theta terms to compute, either 0 or 1.
    :arg float tau: Reduced temperature _TCP/temp(K).
    :arg float dta: Reduced density dflu(kg/m3)/_DCP.
    :returns: Theta value derivative.
    """
    ai, bi, bb, cc, dd, aa, bet, n = _RESIDCOEFFS[3][i]
    betinv = 1./bet
    adta1 = abs(dta-1)
    tht_dd = aa * betinv * (betinv-1) * adta1**(betinv-2)
    return tht_dd

def _delta(i,tau,dta):
    """Calculate _delta term.
    
    Calculate the i(th) _delta parameter (IAPWS-95, table 5). This
    parameter is used in calculating the residual (non-ideal gas)
    component of the Helmholtz free energy of fluid water.
    
    :arg int i: Which of the two _theta terms to compute, either 0 or 1.
    :arg float tau: Reduced temperature _TCP/temp(K).
    :arg float dta: Reduced density dflu(kg/m3)/_DCP.
    :returns: Delta value.
    """
    ai, bi, bb, cc, dd, aa, bet, n = _RESIDCOEFFS[3][i]
    adta1 = abs(dta-1)
    eta = _theta(i,tau,dta)**2 + bb * adta1**(2*ai)
    return eta

def _delta_t(i,tau,dta):
    """Calculate _delta term T-derivative.
    
    Calculate the derivative of the i(th) _delta parameter (IAPWS-95,
    table 5) with respect to reduced temperature. This parameter is used
    in calculating the residual (non-ideal gas) component of the
    Helmholtz free energy of fluid water.
    
    :arg int i: Which of the two _theta terms to compute, either 0 or 1.
    :arg float tau: Reduced temperature _TCP/temp(K).
    :arg float dta: Reduced density dflu(kg/m3)/_DCP.
    :returns: Delta value derivative.
    """
    eta = 2 * _theta(i,tau,dta) * _theta_t(i,tau,dta)
    return eta

def _delta_d(i,tau,dta):
    """Calculate _delta term D-derivative.
    
    Calculate the derivative of the i(th) _delta parameter (IAPWS-95,
    table 5) with respect to reduced density. This parameter is used in
    calculating the residual (non-ideal gas) component of the Helmholtz
    free energy of fluid water.
    
    :arg int i: Which of the two _theta terms to compute, either 0 or 1.
    :arg float tau: Reduced temperature _TCP/temp(K).
    :arg float dta: Reduced density dflu(kg/m3)/_DCP.
    :returns: Delta value derivative.
    """
    ai, bi, bb, cc, dd, aa, bet, n = _RESIDCOEFFS[3][i]
    adta1 = abs(dta-1)
    eta = (2 * _theta(i,tau,dta) * _theta_d(i,tau,dta)
        + bb * 2*ai * (dta-1) * adta1**(2*ai-2))
    return eta

def _delta_tt(i,tau,dta):
    """Calculate _delta term TT-derivative.
    
    Calculate the second derivative of the i(th) _delta parameter
    (IAPWS-95, table 5) with respect to reduced temperature. This
    parameter is used in calculating the residual (non-ideal gas)
    component of the Helmholtz free energy of fluid water.
    
    :arg int i: Which of the two _theta terms to compute, either 0 or 1.
    :arg float tau: Reduced temperature _TCP/temp(K).
    :arg float dta: Reduced density dflu(kg/m3)/_DCP.
    :returns: Delta value derivative.
    """
    eta = 2 * (_theta_t(i,tau,dta)**2 + _theta(i,tau,dta)*_theta_tt(i,tau,dta))
    return eta

def _delta_td(i,tau,dta):
    """Calculate _delta term TD-derivative.
    
    Calculate the mixed derivative of the i(th) _delta parameter
    (IAPWS-95, table 5) with respect to reduced temperature and density.
    This parameter is used in calculating the residual (non-ideal gas)
    component of the Helmholtz free energy of fluid water.
    
    :arg int i: Which of the two _theta terms to compute, either 0 or 1.
    :arg float tau: Reduced temperature _TCP/temp(K).
    :arg float dta: Reduced density dflu(kg/m3)/_DCP.
    :returns: Delta value derivative.
    """
    eta = 2 * (_theta_d(i,tau,dta)*_theta_t(i,tau,dta)
        + _theta(i,tau,dta)*_theta_td(i,tau,dta))
    return eta

def _delta_dd(i,tau,dta):
    """Calculate _delta term DD-derivative.
    
    Calculate the second derivative of the i(th) _delta parameter
    (IAPWS-95, table 5) with respect to reduced density. This parameter
    is used in calculating the residual (non-ideal gas) component of the
    Helmholtz free energy of fluid water.
    
    :arg int i: Which of the two _theta terms to compute, either 0 or 1.
    :arg float tau: Reduced temperature _TCP/temp(K).
    :arg float dta: Reduced density dflu(kg/m3)/_DCP.
    :returns: Delta value derivative.
    """
    ai, bi, bb, cc, dd, aa, bet, n = _RESIDCOEFFS[3][i]
    adta1 = abs(dta-1)
    eta_dd = (2 * _theta(i,tau,dta)*_theta_dd(i,tau,dta)
        + 2 * _theta_d(i,tau,dta)**2
        + bb * 2*ai * (2*ai-1) * adta1**(2*ai-2))
    return eta_dd

def _deltab_t(i,tau,dta):
    """Calculate _deltab term T-derivative.
    
    Calculate the derivative of the i(th) values of _delta^b
    (IAPWS-95, table 5) with respect to reduced temperature. This
    parameter is used in calculating the residual (non-ideal gas)
    component of the Helmholtz free energy of fluid water.
    
    :arg int i: Which of the two _theta terms to compute, either 0 or 1.
    :arg float tau: Reduced temperature _TCP/temp(K).
    :arg float dta: Reduced density dflu(kg/m3)/_DCP.
    :returns: Deltab value derivative.
    """
    ai, bi, bb, cc, dd, aa, bet, n = _RESIDCOEFFS[3][i]
    eta = _delta(i,tau,dta)
    etab = bi * eta**(bi-1) * _delta_t(i,tau,dta)
    return etab

def _deltab_d(i,tau,dta):
    """Calculate _deltab term D-derivative.
    
    Calculate the derivative of the i(th) values of _delta^b
    (IAPWS-95, table 5) with respect to reduced density. This parameter
    is used in calculating the residual (non-ideal gas) component of the
    Helmholtz free energy of fluid water.
    
    :arg int i: Which of the two _theta terms to compute, either 0 or 1.
    :arg float tau: Reduced temperature _TCP/temp(K).
    :arg float dta: Reduced density dflu(kg/m3)/_DCP.
    :returns: Deltab value derivative.
    """
    ai, bi, bb, cc, dd, aa, bet, n = _RESIDCOEFFS[3][i]
    eta = _delta(i,tau,dta)
    etab = bi * eta**(bi-1) * _delta_d(i,tau,dta)
    return etab

def _deltab_tt(i,tau,dta):
    """Calculate _deltab term TT-derivative.
    
    Calculate the second derivative of the i(th) values of _delta^b
    (IAPWS-95, table 5) with respect to reduced temperature. This
    parameter is used in calculating the residual (non-ideal gas)
    component of the Helmholtz free energy of fluid water.
    
    :arg int i: Which of the two _theta terms to compute, either 0 or 1.
    :arg float tau: Reduced temperature _TCP/temp(K).
    :arg float dta: Reduced density dflu(kg/m3)/_DCP.
    :returns: Deltab value derivative.
    """
    ai, bi, bb, cc, dd, aa, bet, n = _RESIDCOEFFS[3][i]
    eta = _delta(i,tau,dta)
    etab = bi * eta**(bi-2) * ((bi-1)*_delta_t(i,tau,dta)**2
        + eta*_delta_tt(i,tau,dta))
    return etab

def _deltab_td(i,tau,dta):
    """Calculate _deltab term TD-derivative.
    
    Calculate the mixed derivative of the i(th) values of _delta^b
    (IAPWS-95, table 5) with respect to reduced temperature and density.
    This parameter is used in calculating the residual (non-ideal gas)
    component of the Helmholtz free energy of fluid water.
    
    :arg int i: Which of the two _theta terms to compute, either 0 or 1.
    :arg float tau: Reduced temperature _TCP/temp(K).
    :arg float dta: Reduced density dflu(kg/m3)/_DCP.
    :returns: Deltab value derivative.
    """
    ai, bi, bb, cc, dd, aa, bet, n = _RESIDCOEFFS[3][i]
    eta = _delta(i,tau,dta)
    etab = bi * eta**(bi-2) * ((bi-1)*_delta_t(i,tau,dta)*_delta_d(i,tau,dta)
        + eta*_delta_td(i,tau,dta))
    return etab

def _deltab_dd(i,tau,dta):
    """Calculate _deltab term DD-derivative.
    
    Calculate the second derivative of the i(th) values of _delta^b
    (IAPWS-95, table 5) with respect to reduced density. This parameter
    is used in calculating the residual (non-ideal gas) component of the
    Helmholtz free energy of fluid water.
    
    :arg int i: Which of the two _theta terms to compute, either 0 or 1.
    :arg float tau: Reduced temperature _TCP/temp(K).
    :arg float dta: Reduced density dflu(kg/m3)/_DCP.
    :returns: Deltab value derivative.
    """
    ai, bi, bb, cc, dd, aa, bet, n = _RESIDCOEFFS[3][i]
    eta = _delta(i,tau,dta)
    etab = bi * eta**(bi-2) * ((bi-1)*_delta_d(i,tau,dta)**2
        + eta*_delta_dd(i,tau,dta))
    return etab

def _psi(i,tau,dta):
    """Calculate _psi term.
    
    Calculate the value of the i(th) _psi parameter (IAPWS-95, table 5).
    This parameter is used in calculating the residual (non-ideal gas)
    component of the Helmholtz free energy of fluid water.
    
    :arg int i: Which of the two _theta terms to compute, either 0 or 1.
    :arg float tau: Reduced temperature _TCP/temp(K).
    :arg float dta: Reduced density dflu(kg/m3)/_DCP.
    :returns: Psi value.
    """
    ai, bi, bb, cc, dd, aa, bet, n = _RESIDCOEFFS[3][i]
    ps = numpy.exp(-cc*(dta-1)**2 - dd*(tau-1)**2)
    return ps

def _psi_t(i,tau,dta):
    """Calculate _psi term T-derivative.
    
    Calculate the derivative of the value of the i(th) _psi parameter
    (IAPWS-95, table 5) with respect to reduced temperature. This
    parameter is used in calculating the residual (non-ideal gas)
    component of the Helmholtz free energy of fluid water.
    
    :arg int i: Which of the two _theta terms to compute, either 0 or 1.
    :arg float tau: Reduced temperature _TCP/temp(K).
    :arg float dta: Reduced density dflu(kg/m3)/_DCP.
    :returns: Psi value derivative.
    """
    ai, bi, bb, cc, dd, aa, bet, n = _RESIDCOEFFS[3][i]
    ps = -2 * dd * (tau-1) * _psi(i,tau,dta)
    return ps

def _psi_d(i,tau,dta):
    """Calculate _psi term D-derivative.
    
    Calculate the derivative of the value of the i(th) _psi parameter
    (IAPWS-95, table 5) with respect to reduced density. This parameter
    is used in calculating the residual (non-ideal gas) component of the
    Helmholtz free energy of fluid water.
    
    :arg int i: Which of the two _theta terms to compute, either 0 or 1.
    :arg float tau: Reduced temperature _TCP/temp(K).
    :arg float dta: Reduced density dflu(kg/m3)/_DCP.
    :returns: Psi value derivative.
    """
    ai, bi, bb, cc, dd, aa, bet, n = _RESIDCOEFFS[3][i]
    ps = -2 * cc * (dta-1) * _psi(i,tau,dta)
    return ps

def _psi_tt(i,tau,dta):
    """Calculate _psi term TT-derivative.
    
    Calculate the second derivative of the value of the i(th) _psi
    parameter (IAPWS-95, table 5) with respect to reduced temperature.
    This parameter is used in calculating the residual (non-ideal gas)
    component of the Helmholtz free energy of fluid water.
    
    :arg int i: Which of the two _theta terms to compute, either 0 or 1.
    :arg float tau: Reduced temperature _TCP/temp(K).
    :arg float dta: Reduced density dflu(kg/m3)/_DCP.
    :returns: Psi value derivative.
    """
    ai, bi, bb, cc, dd, aa, bet, n = _RESIDCOEFFS[3][i]
    ps = 2 * dd * (2*dd * (tau-1)**2 - 1)
    ps *= _psi(i,tau,dta)
    return ps

def _psi_td(i,tau,dta):
    """Calculate _psi term TD-derivative.
    
    Calculate the mixed derivative of the value of the i(th) _psi
    parameter (IAPWS-95, table 5) with respect to reduced temperature
    and density. This parameter is used in calculating the residual
    (non-ideal gas) component of the Helmholtz free energy of fluid
    water.
    
    :arg int i: Which of the two _theta terms to compute, either 0 or 1.
    :arg float tau: Reduced temperature _TCP/temp(K).
    :arg float dta: Reduced density dflu(kg/m3)/_DCP.
    :returns: Psi value derivative.
    """
    ai, bi, bb, cc, dd, aa, bet, n = _RESIDCOEFFS[3][i]
    ps = 4 * cc * dd * (dta-1) * (tau-1) * _psi(i,tau,dta)
    return ps

def _psi_dd(i,tau,dta):
    """Calculate _psi term DD-derivative.
    
    Calculate the second derivative of the value of the i(th) _psi
    parameter (IAPWS-95, table 5) with respect to reduced density. This
    parameter is used in calculating the residual (non-ideal gas)
    component of the Helmholtz free energy of fluid water.
    
    :arg int i: Which of the two _theta terms to compute, either 0 or 1.
    :arg float tau: Reduced temperature _TCP/temp(K).
    :arg float dta: Reduced density dflu(kg/m3)/_DCP.
    :returns: Psi value derivative.
    """
    ai, bi, bb, cc, dd, aa, bet, n = _RESIDCOEFFS[3][i]
    ps = 2 * cc * (2 * cc * (dta-1)**2 - 1)
    ps *= _psi(i,tau,dta)
    return ps

def _phir(tau,dta):
    """Calculate fluid water potential residual component.
    
    Calculate the residual (non-ideal gas) component of the Helmholtz
    potential (scaled free energy) for fluid water.
    
    :arg float tau: Reduced temperature _TCP/temp(K).
    :arg float dta: Reduced density dflu(kg/m3)/_DCP.
    """
    phi = 0.0
    for (d,t,n) in _RESIDCOEFFS[0]:
        phi += n * dta**d * tau**t
    for (c,d,t,n) in _RESIDCOEFFS[1]:
        eterm = numpy.exp(-dta**c)
        phi += n * dta**d * eterm * tau**t
    for (d,t,alf,eps,bet,gam,n) in _RESIDCOEFFS[2]:
        eterm = numpy.exp(-alf * (dta-eps)**2 - bet*(tau-gam)**2)
        phi += n * dta**d * tau**t * eterm
    for (i,gam) in enumerate(_RESIDCOEFFS[3]):
        ai, bi, bb, cc, dd, aa, bet, n = gam
        ps = _psi(i,tau,dta)
        eta = _delta(i,tau,dta)
        phi += n * eta**bi * dta * ps
    return phi

def _phir_t(tau,dta):
    """Calculate fluid water potential residual component T-derivative.
    
    Calculate the derivative of the residual (non-ideal gas) component
    of the Helmholtz potential (scaled free energy) for fluid water with
    respect to reduced temperature.
    
    :arg float tau: Reduced temperature _TCP/temp(K).
    :arg float dta: Reduced density dflu(kg/m3)/_DCP.
    """
    phi = 0.0
    for (d,t,n) in _RESIDCOEFFS[0]:
        phi += n * dta**d * t*tau**(t-1)
    for (c,d,t,n) in _RESIDCOEFFS[1]:
        term = dta**d * numpy.exp(-dta**c) * t * tau**(t-1)
        phi += n * term
    for (d,t,alf,eps,bet,gam,n) in _RESIDCOEFFS[2]:
        eterm = numpy.exp(-alf * (dta-eps)**2 - bet*(tau-gam)**2)
        eder = -2*bet*(tau-gam)
        term = dta**d * tau**(t-1) * eterm * (t + tau * eder)
        phi += n * term
    for (i,gam) in enumerate(_RESIDCOEFFS[3]):
        ai, bi, bb, cc, dd, aa, bet, n = gam
        ps = _psi(i,tau,dta)
        ps_t = _psi_t(i,tau,dta)
        eta = _delta(i,tau,dta)
        etab_t = _deltab_t(i,tau,dta)
        term = dta * (etab_t*ps + eta**bi * ps_t)
        phi += n * term
    return phi

def _phir_d(tau,dta):
    """Calculate fluid water potential residual component D-derivative.
    
    Calculate the derivative of the residual (non-ideal gas) component
    of the Helmholtz potential (scaled free energy) for fluid water with
    respect to reduced density.
    
    :arg float tau: Reduced temperature _TCP/temp(K).
    :arg float dta: Reduced density dflu(kg/m3)/_DCP.
    """
    phi = 0.0
    for (d,t,n) in _RESIDCOEFFS[0]:
        phi += n * d*dta**(d-1) * tau**t
    for (c,d,t,n) in _RESIDCOEFFS[1]:
        dtac = dta**c
        eterm = numpy.exp(-dtac)
        term = tau**t * dta**(d-1) * eterm * (d - c*dtac)
        phi += n * term
    for (d,t,alf,eps,bet,gam,n) in _RESIDCOEFFS[2]:
        eterm = numpy.exp(-alf * (dta-eps)**2 - bet*(tau-gam)**2)
        eder = -2 * alf * (dta-eps)
        term = tau**t * dta**(d-1) * eterm * (d + dta * eder)
        phi += n * term
    for (i,gam) in enumerate(_RESIDCOEFFS[3]):
        ai, bi, bb, cc, dd, aa, bet, n = gam
        eta = _delta(i,tau,dta)
        etab = eta**bi
        etab_d = _deltab_d(i,tau,dta)
        ps = _psi(i,tau,dta)
        ps_d = _psi_d(i,tau,dta)
        _deltab = _delta(i,tau,dta)**bi
        term = etab_d*dta*ps + etab*ps + etab*dta*ps_d
        phi += n * term
    return phi

def _phir_tt(tau,dta):
    """Calculate fluid water potential residual component TT-derivative.
    
    Calculate the second derivative of the residual (non-ideal gas)
    component of the Helmholtz potential (scaled free energy) for fluid
    water with respect to reduced temperature.
    
    :arg float tau: Reduced temperature _TCP/temp(K).
    :arg float dta: Reduced density dflu(kg/m3)/_DCP.
    """
    phi = 0.0
    for (d,t,n) in _RESIDCOEFFS[0]:
        phi += n * dta**d * t*(t-1)*tau**(t-2)
    for (c,d,t,n) in _RESIDCOEFFS[1]:
        phi += n * numpy.exp(-dta**c) * dta**d * t*(t-1)*tau**(t-2)
    for (d,t,alf,eps,bet,gam,n) in _RESIDCOEFFS[2]:
        eterm = numpy.exp(-alf * (dta-eps)**2 - bet*(tau-gam)**2)
        eder1 = -2*bet*(tau-gam)
        eder2 = -2*bet
        term = dta**d * tau**(t-2) * eterm * (t*(t-1) + 2*t*tau*eder1
            + tau**2*eder1**2 + tau**2*eder2)
        phi += n * term
    for (i,gam) in enumerate(_RESIDCOEFFS[3]):
        ai, bi, bb, cc, dd, aa, bet, n = gam
        eta = _delta(i,tau,dta)
        etab_t = _deltab_t(i,tau,dta)
        etab_tt = _deltab_tt(i,tau,dta)
        ps = _psi(i,tau,dta)
        ps_t = _psi_t(i,tau,dta)
        ps_tt = _psi_tt(i,tau,dta)
        term = dta * (etab_tt*ps + 2*etab_t*ps_t + eta**bi * ps_tt)
        phi += n * term
    return phi

def _phir_td(tau,dta):
    """Calculate fluid water potential residual component TD-derivative.
    
    Calculate the mixed derivative of the residual (non-ideal gas)
    component of the Helmholtz potential (scaled free energy) for fluid
    water with respect to reduced temperature and density.
    
    :arg float tau: Reduced temperature _TCP/temp(K).
    :arg float dta: Reduced density dflu(kg/m3)/_DCP.
    """
    phi = 0.0
    for (d,t,n) in _RESIDCOEFFS[0]:
        phi += n * d*dta**(d-1) * t*tau**(t-1)
    for (c,d,t,n) in _RESIDCOEFFS[1]:
        dtac = dta**c
        eterm = numpy.exp(-dtac)
        term = t*tau**(t-1) * dta**(d-1) * eterm * (d - c*dtac)
        phi += n * term
    for (d,t,alf,eps,bet,gam,n) in _RESIDCOEFFS[2]:
        eterm = numpy.exp(-alf * (dta-eps)**2 - bet*(tau-gam)**2)
        edert = -2*bet*(tau-gam)
        ederd = -2*alf*(dta-eps)
        term = tau**(t-1) * dta**(d-1) * eterm
        term *= (t + tau*edert) * (d + dta*ederd)
        phi += n * term
    for (i,gam) in enumerate(_RESIDCOEFFS[3]):
        ai, bi, bb, cc, dd, aa, bet, n = gam
        eta = _delta(i,tau,dta)
        etab = eta**bi
        etab_t = _deltab_t(i,tau,dta)
        etab_d = _deltab_d(i,tau,dta)
        etab_td = _deltab_td(i,tau,dta)
        ps = _psi(i,tau,dta)
        ps_t = _psi_t(i,tau,dta)
        ps_d = _psi_d(i,tau,dta)
        ps_td = _psi_td(i,tau,dta)
        term = etab_td*dta*ps + etab_t*(ps + dta*ps_d)
        term += (etab_d*dta + etab)*ps_t + etab*dta*ps_td
        phi += n * term
    return phi

def _phir_dd(tau,dta):
    """Calculate fluid water potential residual component DD-derivative.
    
    Calculate the second derivative of the residual (non-ideal gas)
    component of the Helmholtz potential (scaled free energy) for fluid
    water with respect to reduced density.
    
    :arg float tau: Reduced temperature _TCP/temp(K).
    :arg float dta: Reduced density dflu(kg/m3)/_DCP.
    """
    phi = 0.0
    for (d,t,n) in _RESIDCOEFFS[0]:
        phi += n * d*(d-1)*dta**(d-2) * tau**t
    for (c,d,t,n) in _RESIDCOEFFS[1]:
        dtac = dta**c
        eterm = numpy.exp(-dtac)
        term = tau**t * dta**(d-2) * eterm * ((d-1-c*dtac)*(d-c*dtac)
            - c**2*dtac)
        phi += n * term
    for (d,t,alf,eps,bet,gam,n) in _RESIDCOEFFS[2]:
        eterm = numpy.exp(-alf*(dta-eps)**2 - bet*(tau-gam)**2)
        eder1 = -2*alf*(dta-eps)
        eder2 = -2*alf
        term = tau**t * dta**(d-2) * eterm * (d*(d-1) + 2*d*dta*eder1
            + dta**2*(eder1**2 + eder2))
        phi += n * term
    for (i,gam) in enumerate(_RESIDCOEFFS[3]):
        ai, bi, bb, cc, dd, aa, bet, n = gam
        eta = _delta(i,tau,dta)
        etab = eta**bi
        etab_d = _deltab_d(i,tau,dta)
        etab_dd = _deltab_dd(i,tau,dta)
        ps = _psi(i,tau,dta)
        ps_d = _psi_d(i,tau,dta)
        ps_dd = _psi_dd(i,tau,dta)
        term = (etab_dd*dta + 2*etab_d) * ps
        term += 2 * (etab_d*dta + etab) * ps_d
        term += etab * dta * ps_dd
        phi += n * term
    return phi


### Public function
def flu_f(drvt,drvd,temp,dflu,chkbnd=False):
    """Calculate fluid water Helmholtz free energy.
    
    Calculate the specific Helmholtz free energy of fluid water or its
    derivatives with respect to temperature and density. Derivatives up
    to second order are available.
    
    :arg int drvt: Number of temperature derivatives.
    :arg int drvd: Number of density derivatives.
    :arg float temp: Temperature in K.
    :arg float dflu: Fluid water density in kg/m3.
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :returns: Helmholtz free energy in units of
        (J/kg) / K^drvt / (kg/m3)^drvd.
    :raises ValueError: If either temp or dflu are nonpositive.
    :raises RuntimeWarning: If temp or dflu are outside the recommended
        bounds and chkbnd is True.
    :raises ValueError: If drvt<0, drvd<0, or drvt+drvd>2.
    
    :Examples:
    
    >>> flu_f(0,0,300.,1000.)
    -5351.74115204
    >>> flu_f(1,0,300.,1000.)
    -390.904170767
    >>> flu_f(0,1,300.,1000.)
    7.83300135597
    >>> flu_f(2,0,300.,1000.)
    -13.6840204925
    >>> flu_f(1,1,300.,1000.)
    0.639359046588
    >>> flu_f(0,2,300.,1000.)
    2.24824656167
    """
    constants0.chkflubnds(temp,dflu,chkbnd=chkbnd)
    tau = _TCP / temp
    dta = dflu / _DCP
    rt = _RWAT * temp
    
    # Run through each derivative case
    if (drvt,drvd) == (0,0):
        f = _RWAT * temp * (_phi0(tau,dta) + _phir(tau,dta))
    elif (drvt,drvd) == (1,0):
        phi = _phi0(tau,dta) + _phir(tau,dta)
        phi_t = _phi0_t(tau,dta) + _phir_t(tau,dta)
        f = _RWAT * (phi - tau*phi_t)
    elif (drvt,drvd) == (0,1):
        f = _RWAT * temp / _DCP * (_phi0_d(tau,dta) + _phir_d(tau,dta))
    elif (drvt,drvd) == (2,0):
        phi_tt = _phi0_tt(tau,dta) + _phir_tt(tau,dta)
        f = _RWAT / temp * tau**2 * phi_tt
    elif (drvt,drvd) == (1,1):
        phi_d = _phi0_d(tau,dta) + _phir_d(tau,dta)
        phi_td = _phi0_td(tau,dta) + _phir_td(tau,dta)
        f = _RWAT / _DCP * (phi_d - tau*phi_td)
    elif (drvt,drvd) == (0,2):
        phi_dd = _phi0_dd(tau,dta) + _phir_dd(tau,dta)
        f = _RWAT * temp / _DCP**2 * phi_dd
    else:
        errmsg = 'Derivatives {0} not recognized'.format((drvt,drvd))
        raise ValueError(errmsg)
    return f

def chkiapws95table6(printresult=True,chktol=_CHKTOL):
    """Check accuracy against IAPWS 1995 table 6.
    
    Evaluate the functions in this module and compare to reference
    values from IAPWS 1995, table 6. These values include separate
    contributions from the ideal gas and residual components of the
    Helmholtz potential.
    
    :arg bool printresult: If True (default) and any results are outside
        of the given tolerance, then the function name, reference value,
        result value, and relative error are printed.
    :arg float chktol: Tolerance to use when choosing to print results
        (default _CHKTOL).
    :returns: Tester instances containing the functions, arguments,
        reference values, results, and relative errors from the tests.
        The first instance is for the ideal gas functions and the second
        for residual functions.
    :rtype: (Tester,Tester)
    """
    from tester import Tester
    TCHK = 500.
    DCHK = 838.025
    TAUCHK = _TCP/TCHK
    DTACHK = DCHK/_DCP
    fargs = (TAUCHK,DTACHK)
    argfmt = '({0:6.4f},{1:6.4f})'
    
    # Initialize Tester instances
    idealfuns = [_phi0,_phi0_d,_phi0_dd,_phi0_t,_phi0_td,_phi0_tt]
    idealfnames = ['phi0','phi0_d','phi0_dd','phi0_t','phi0_td','phi0_tt']
    idealrefs = [0.20479773347960e1,0.38423674711375,-0.147637877832556,
        0.90461110617524e1,0.,-0.193249185013052e1]
    idealtest = Tester(idealfuns,fargs,idealrefs,idealfnames,argfmt)
    residfuns = [_phir,_phir_d,_phir_dd,_phir_t,_phir_td,_phir_tt]
    residfnames = ['phir','phir_d','phir_dd','phir_t','phir_td','phir_tt']
    residrefs = [-0.34269320568156e1,-0.36436665036388,0.85606370097461,
        -0.58140343523842e1,-0.11217691467031e1,-0.223440736884336e1]
    residtest = Tester(residfuns,fargs,residrefs,residfnames,argfmt)
    
    # Run tester instances and print results
    idealtest.run()
    residtest.run()
    if printresult:
        msg = 'Fluid water ideal gas component'
        print(msg)
        idealtest.printresults(chktol=chktol)
        msg = 'Fluid water residual component'
        print(msg)
        residtest.printresults(chktol=chktol)
    return idealtest, residtest


## Main function: Check tables
if __name__ == '__main__':
    idealtest, residtest = chkiapws95table6()

