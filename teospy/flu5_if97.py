"""Fluid water Gibbs energy functions from IAPWS-97.

This module provides the Gibbs free energies of both liquid water and
water vapour as defined in IAPWS 1997. This module can also be called as
a function,

    python flu5_if97.py

which will compare the results of this module to reference values from
IAPWS 1997 tables 5 and 15.

:Examples:

>>> liq_volume(300.,3e6)
1.00215168e-3
>>> liq_cp(300.,3e6)
4.17301218e3
>>> liq_soundspeed(300.,3e6)
1.50773921e3
>>> vap_volume(300.,3.5e3)
39.4913866e2
>>> vap_cp(300.,3.5e3)
1.91300162e3
>>> vap_soundspeed(300.,3.5e3)
427.920172

:Functions:

* :func:`liq_g`: Liquid water Gibbs free energy with derivatives.
* :func:`vap_g`: Water vapour Gibbs free energy with derivatives.
* :func:`liq_cp`: Liquid water isobaric heat capacity.
* :func:`liq_density`: Liquid water density.
* :func:`liq_enthalpy`: Liquid water enthalpy.
* :func:`liq_entropy`: Liquid water entropy.
* :func:`liq_internalenergy`: Liquid water internal energy.
* :func:`liq_soundspeed`: Liquid water sound speed.
* :func:`liq_volume`: Liquid water specific volume.
* :func:`vap_cp`: Water vapour isobaric heat capacity.
* :func:`vap_density`: Water vapour density.
* :func:`vap_enthalpy`: Water vapour enthalpy.
* :func:`vap_entropy`: Water vapour entropy.
* :func:`vap_internalenergy`: Water vapour internal energy.
* :func:`vap_soundspeed`: Water vapour sound speed.
* :func:`vap_volume`: Water vapour specific volume.
* :func:`chkiapws97table`: Check module against IAPWS 1997.
* :func:`chkiapws97table5`: Check module against IAPWS 1997, table 5
  (liquid water).
* :func:`chkiapws97table15`: Check module against IAPWS 1997, table 15
  (water vapour).

"""

__all__ = ['liq_g','vap_g',
    'liq_cp','liq_density','liq_enthalpy','liq_entropy','liq_internalenergy',
    'liq_soundspeed','liq_volume','vap_cp','vap_density','vap_enthalpy',
    'vap_entropy','vap_internalenergy','vap_soundspeed','vap_volume',
    'chkiapws97table','chkiapws97table5','chkiapws97table15']

import numpy
import constants0

_CHKTOL = constants0.CHKTOL
_C_LIQ = ((1386., 1.653e7, 461.526, 1.222, 7.1),
    (
        ( 0, -2, 0.14632971213167),    ( 0, -1,-0.84548187169114),
        ( 0,  0,-3.756360367204),      ( 0,  1, 3.3855169168385),
        ( 0,  2,-0.95791963387872),    ( 0,  3, 0.15772038513228),
        ( 0,  4,-0.016616417199501),   ( 0,  5, 8.1214629983568e-4),
        ( 1, -9, 2.8319080123804e-4),  ( 1, -7,-6.0706301565874e-4),
        ( 1, -1,-0.018990068218419),   ( 1,  0,-0.032529748770505),
        ( 1,  1,-0.021841717175414),   ( 1,  3,-5.283835796993e-5),
        ( 2, -3,-4.7184321073267e-4),  ( 2,  0,-3.0001780793026e-4),
        ( 2,  1, 4.7661393906987e-5),  ( 2,  3,-4.4141845330846e-6),
        ( 2, 17,-7.2694996297594e-16), ( 3, -4,-3.1679644845054e-5),
        ( 3,  0,-2.8270797985312e-6),  ( 3,  6,-8.5205128120103e-10),
        ( 4, -5,-2.2425281908e-6),     ( 4, -2,-6.5171222895601e-7),
        ( 4, 10,-1.4341729937924e-13), ( 5, -8,-4.0516996860117e-7),
        ( 8,-11,-1.2734301741641e-9),  ( 8, -6,-1.7424871230634e-10),
        (21,-29,-6.8762131295531e-19), (23,-31, 1.4478307828521e-20),
        (29,-38, 2.6335781662795e-23), (30,-39,-1.1947622640071e-23),
        (31,-40, 1.8228094581404e-24), (32,-41,-9.3537087292458e-26)
    ))
_C_VAP = ((540., 1e6, 461.526, 0.5),
    (
        ( 0,-9.6927686500217),   ( 1,10.086655968018),
        (-5,-0.005608791128302), (-4, 0.071452738081455),
        (-3,-0.40710498223928),  (-2, 1.4240819171444),
        (-1,-4.383951131945),    ( 2,-0.28408632460772),
        ( 3, 0.021268463753307)
    ),
    (
        ( 1, 0,-1.7731742473213e-3),  ( 1, 1,-0.017834862292358),
        ( 1, 2,-0.045996013696365),   ( 1, 3,-0.057581259083432),
        ( 1, 6,-0.05032527872793),    ( 2, 1,-3.3032641670203e-5),
        ( 2, 2,-1.8948987516315e-4),  ( 2, 4,-3.9392777243355e-3),
        ( 2, 7,-0.043797295650573),   ( 2,36,-2.6674547914087e-5),
        ( 3, 0, 2.0481737692309e-8),  ( 3, 1, 4.3870667284435e-7),
        ( 3, 3,-3.227767723857e-5),   ( 3, 6,-1.5033924542148e-3),
        ( 3,35,-0.040668253562649),   ( 4, 1,-7.8847309559367e-10),
        ( 4, 2, 1.2790717852285e-8),  ( 4, 3, 4.8225372718507e-7),
        ( 5, 7, 2.2922076337661e-6),  ( 6, 3,-1.6714766451061e-11),
        ( 6,16,-2.1171472321355e-3),  ( 6,35,-23.895741934104),
        ( 7, 0,-5.905956432427e-18),  ( 7,11,-1.2621808899101e-6),
        ( 7,25,-0.038946842435739),   ( 8, 8, 1.1256211360459e-11),
        ( 8,36,-8.2311340897998),     ( 9,13, 1.9809712802088e-8),
        (10, 4, 1.0406965210174e-19), (10,10,-1.0234747095929e-13),
        (10,14,-1.0018179379511e-9),  (16,29,-8.0882908646985e-11),
        (16,50, 0.10693031879409),    (18,57,-0.33662250574171),
        (20,20, 8.9185845355421e-25), (20,35, 3.0629316876232e-13),
        (20,48,-4.2002467698208e-6),  (21,21,-5.9056029685639e-26),
        (22,53, 3.7826947613457e-6),  (23,39,-1.2768608934681e-15),
        (24,26, 7.3087610595061e-29), (24,40, 5.5414715350778e-17),
        (24,58,-9.436970724121e-7)
    ))


## Liquid water Gibbs function
def _gammaliq(drvt,drvp,tau,psi):
    """Calculate liquid water Gibbs potential.
    
    Calculate the Gibbs potential (scaled Gibbs free energy) of liquid
    water or its derivatives with respect to reduced temperature and
    reduced pressure.
    
    :arg int drvt: Number of reduced temperature derivatives.
    :arg int drvp: Number of reduced pressure derivatives.
    :arg float tau: Reduced temperature 1386/temp(K).
    :arg float psi: Reduced pressure pres(Pa)/1.653e7.
    :returns: Gibbs potential, unitless.
    :raises ValueError: if tau < 1.222 or psi is nonpositive.
    :raises ValueError: If drvt or drvp are negative.
    """
    TAU0, PSI0 = _C_LIQ[0][-2:]
    if tau < TAU0:
        errmsg = 'Scaled temperature must be larger than {0}'.format(TAU0)
        raise ValueError(errmsg)
    if psi <= 0:
        errmsg = 'Scaled pressure must be positive'
        raise ValueError(errmsg)
    if drvt < 0 or drvp < 0:
        errmsg = 'Derivatives must be nonnegative'
        raise ValueError(errmsg)
    
    pp = PSI0 - psi
    tt = tau - TAU0
    gam = 0.
    for (i,j,n) in _C_LIQ[1]:
        pwrt = tt**(j-drvt)
        for k in range(drvt):
            pwrt *= (j-k)
        if pp == 0:
            if i == 1:
                pwrp = 1.
            else:
                pwrp = 0.
        else:
            pwrp = pp**(i-drvp)
        for k in range(drvp):
            pwrp *= -(i-k)
        gam += n * pwrt * pwrp
    return gam

def liq_g(drvt,drvp,temp,pres):
    """Calculate liquid water Gibbs energy.
    
    Calculate the specific Gibbs free energy of liquid water or its
    derivatives with respect to temperature and pressure using the
    IAPWS-97 formulation.
    
    :arg int drvt: Number of temperature derivatives.
    :arg int drvp: Number of pressure derivatives.
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :returns: Gibbs energy in units of
        (J/kg) / K^drvt / Pa^drvp.
    :raises ValueError: If temp or pres are nonpositive, or if
        temp > 1134 K.
    :raises ValueError: If drvt or drvp are negative, or if
        (drvt+drvp)>2.
    """
    if temp <= 0.:
        errmsg = 'Temperature must be positive'
        raise ValueError(errmsg)
    if pres <= 0.:
        errmsg = 'Pressure must be positive'
        raise ValueError(errmsg)
    TRED, PRED, RWAT, TAU0 = _C_LIQ[0][:4]
    TMAX = TRED/TAU0
    if temp > TMAX:
        errmsg = ('Temperature {0} must be less than {1} K').format(temp,TMAX)
        raise ValueError(errmsg)
    psi = pres / PRED
    tau = TRED / temp
    rt = RWAT * temp
    
    if (drvt,drvp) == (0,0):
        gam = _gammaliq(0,0,tau,psi)
        gliq = rt*gam
    elif (drvt,drvp) == (0,1):
        gam_p = _gammaliq(0,1,tau,psi)
        gliq = rt*gam_p/PRED
    elif (drvt,drvp) == (0,2):
        gam_pp = _gammaliq(0,2,tau,psi)
        gliq = rt*gam_pp/PRED**2
    elif (drvt,drvp) == (1,0):
        gam = _gammaliq(0,0,tau,psi)
        gam_t = _gammaliq(1,0,tau,psi)
        gliq = RWAT*(gam - tau*gam_t)
    elif (drvt,drvp) == (1,1):
        gam_p = _gammaliq(0,1,tau,psi)
        gam_tp = _gammaliq(1,1,tau,psi)
        gliq = RWAT*(gam_p - tau*gam_tp)/PRED
    elif (drvt,drvp) == (2,0):
        gam_tt = _gammaliq(2,0,tau,psi)
        gliq = RWAT*tau**2 * gam_tt/temp
    else:
        # Should not have made it this far!
        errmsg = 'Derivatives {0} not recognized'.format((drvt,drvp))
        raise ValueError(errmsg)
    return gliq


## Water vapour Gibbs function
def _gammaideal(drvt,drvp,tau,psi):
    """Calculate water vapour ideal Gibbs potential.
    
    Calculate the ideal gas component of the Gibbs potential (scaled
    Gibbs free energy) for water vapour or its derivatives with respect
    to reduced temperature and reduced pressure.
    
    :arg int drvt: Number of reduced temperature derivatives.
    :arg int drvp: Number of reduced pressure derivatives.
    :arg float tau: Reduced temperature 540/temp(K).
    :arg float psi: Reduced pressure pres(Pa)/1e6.
    :returns: Gibbs potential, unitless.
    :raises ValueError: If tau or psi are nonpositive.
    :raises ValueError: If drvt or drvp are negative.
    """
    if tau <= 0:
        errmsg = 'Scaled temperature {0} must be positive'.format(tau)
        raise ValueError(errmsg)
    if psi <= 0:
        errmsg = 'Scaled pressure {0} must be positive'.format(psi)
        raise ValueError(errmsg)
    if drvt < 0 or drvp < 0:
        errmsg = 'Derivatives {0} must be nonnegative'.format((drvt,drvp))
        raise ValueError(errmsg)
    
    if drvt == 0:
        if drvp == 0:
            gam = numpy.log(psi)
        else:
            gam = psi**(-drvp)
            for k in range(2,drvp+1):
                gam *= (1-k)
    else:
        gam = 0.
    if drvp == 0:
        for (j,n) in _C_VAP[1]:
            pwrt = tau**(j-drvt)
            for k in range(drvt):
                pwrt *= j-k
            gam += n*pwrt
    return gam

def _gammaresid(drvt,drvp,tau,psi):
    """Calculate water vapour residual Gibbs potential.
    
    Calculate the residual (non-ideal gas) component of the Gibbs
    potential (scaled Gibbs free energy) of water vapour or its
    derivatives with respect to reduced temperature and reduced
    pressure.
    
    :arg int drvt: Number of reduced temperature derivatives.
    :arg int drvp: Number of reduced pressure derivatives.
    :arg float tau: Reduced temperature 540/temp(K).
    :arg float psi: Reduced pressure pres(Pa)/1e6.
    :returns: Gibbs potential, unitless.
    :raises ValueError: If psi is nonpositive or tau < 0.5.
    :raises ValueError: If drvt or drvp are negative.
    """
    TAU0 = _C_VAP[0][-1]
    if tau <= TAU0:
        errmsg = 'Scaled temperature must be larger than {0}'.format(TAU0)
        raise ValueError(errmsg)
    if psi <= 0:
        errmsg = 'Scaled pressure {0} must be positive'.format(psi)
        raise ValueError(errmsg)
    if drvt < 0 or drvp < 0:
        errmsg = 'Derivatives {0} must be nonnegative'.format((drvt,drvp))
        raise ValueError(errmsg)
    
    tt = tau - TAU0
    gam = 0.
    for (i,j,n) in _C_VAP[2]:
        pwrt = tt**(j-drvt)
        for k in range(drvt):
            pwrt *= j-k
        pwrp = psi**(i-drvp)
        for k in range(drvp):
            pwrp *= i-k
        gam += n * pwrt * pwrp
    return gam

def vap_g(drvt,drvp,temp,pres):
    """Calculate water vapour Gibbs energy.
    
    Calculate the specific Gibbs free energy of water vapour or its
    derivatives with respect to temperature and pressure using the
    IAPWS-97 formulation.
    
    :arg int drvt: Number of temperature derivatives.
    :arg int drvp: Number of pressure derivatives.
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :returns: Gibbs energy in units of
        (J/kg) / K^drvt / Pa^drvp.
    :raises ValueError: If temp or pres are nonpositive, or if
        temp > 1080 K.
    :raises ValueError: If drvt or drvp are negative, or if
        (drvt+drvp)>2.
    """
    if temp <= 0.:
        errmsg = 'Temperature {0} must be positive'.format(temp)
        raise ValueError(errmsg)
    if pres <= 0.:
        errmsg = 'Pressure {0} must be positive'.format(pres)
        raise ValueError(errmsg)
    
    # Calculate reduced variables
    TRED, PRED, RWAT, TAU0 = _C_VAP[0]
    TMAX = TRED/TAU0
    if temp > TMAX:
        errmsg = ('Temperature {0} must be less than {1} K').format(temp,TMAX)
        raise ValueError(errmsg)
    psi = pres/PRED
    tau = TRED/temp
    rt = RWAT*temp
    
    if (drvt,drvp) == (0,0):
        gam0 = _gammaideal(0,0,tau,psi)
        gamr = _gammaresid(0,0,tau,psi)
        gvap = rt*(gam0 + gamr)
    elif (drvt,drvp) == (0,1):
        gam0_p = _gammaideal(0,1,tau,psi)
        gamr_p = _gammaresid(0,1,tau,psi)
        gvap = rt*(gam0_p + gamr_p)/PRED
    elif (drvt,drvp) == (0,2):
        gam0_pp = _gammaideal(0,2,tau,psi)
        gamr_pp = _gammaresid(0,2,tau,psi)
        gvap = rt*(gam0_pp + gamr_pp)/PRED**2
    elif (drvt,drvp) == (1,0):
        gam0 = _gammaideal(0,0,tau,psi)
        gam0_t = _gammaideal(1,0,tau,psi)
        gamr = _gammaresid(0,0,tau,psi)
        gamr_t = _gammaresid(1,0,tau,psi)
        gvap = RWAT*((gam0 + gamr) - tau*(gam0_t + gamr_t))
    elif (drvt,drvp) == (1,1):
        gam0_p = _gammaideal(0,1,tau,psi)
        gam0_tp = _gammaideal(1,1,tau,psi)
        gamr_p = _gammaresid(0,1,tau,psi)
        gamr_tp = _gammaresid(1,1,tau,psi)
        gvap = RWAT*((gam0_p + gamr_p) - tau*(gam0_tp + gamr_tp))/PRED
    elif (drvt,drvp) == (2,0):
        gam0_tt = _gammaideal(2,0,tau,psi)
        gamr_tt = _gammaresid(2,0,tau,psi)
        gvap = RWAT*tau**2 * (gam0_tt + gamr_tt)/temp
    else:
        # Should not have made it this far!
        errmsg = 'Derivatives {0} not recognized'.format((drvt,drvp))
        raise ValueError(errmsg)
    return gvap


## Thermodynamic properties
def liq_cp(temp,pres):
    """Calculate liquid water isobaric heat capacity.
    
    Calculate the isobaric (constant pressure) heat capacity of liquid
    water using the IAPWS-97 formulation.
    
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :returns: Heat capacity in J/kg/K.
    """
    g_tt = liq_g(2,0,temp,pres)
    cp = -temp * g_tt
    return cp

def liq_density(temp,pres):
    """Calculate liquid water density.
    
    Calculate the density of liquid water using the IAPWS-97
    formulation.
    
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :returns: Density in kg/m3.
    """
    vliq = liq_g(0,1,temp,pres)
    dliq = vliq**(-1)
    return dliq

def liq_enthalpy(temp,pres):
    """Calculate liquid water enthalpy.
    
    Calculate the specific enthalpy of liquid water using the IAPWS-97
    formulation.
    
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :returns: Enthalpy in J/kg.
    """
    g = liq_g(0,0,temp,pres)
    g_t = liq_g(1,0,temp,pres)
    h = g - temp*g_t
    return h

def liq_entropy(temp,pres):
    """Calculate liquid water entropy.
    
    Calculate the specific entropy of liquid water using the IAPWS-97
    formulation.
    
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :returns: Entropy in J/kg/K.
    """
    g_t = liq_g(1,0,temp,pres)
    s = -g_t
    return s

def liq_internalenergy(temp,pres):
    """Calculate liquid water internal energy.
    
    Calculate the specific internal energy of liquid water using the
    IAPWS-97 formulation.
    
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :returns: Internal energy in J/kg.
    """
    g = liq_g(0,0,temp,pres)
    g_t = liq_g(1,0,temp,pres)
    g_p = liq_g(0,1,temp,pres)
    u = g - temp*g_t - pres*g_p
    return u

def liq_soundspeed(temp,pres):
    """Calculate liquid water sound speed.
    
    Calculate the speed of sound in liquid water using the IAPWS-97
    formulation.
    
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :returns: Sound speed in m/s.
    """
    g_p = liq_g(0,1,temp,pres)
    g_tt = liq_g(2,0,temp,pres)
    g_tp = liq_g(1,1,temp,pres)
    g_pp = liq_g(0,2,temp,pres)
    csqinv = (g_tp**2/g_tt - g_pp) / g_p**2
    c = csqinv**(-.5)
    return c

def liq_volume(temp,pres):
    """Calculate liquid water specific volume.
    
    Calculate the specific volume of liquid water using the IAPWS-97
    formulation.
    
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :returns: Specific volume in m3/kg.
    """
    v = liq_g(0,1,temp,pres)
    return v

def vap_cp(temp,pres):
    """Calculate water vapour isobaric heat capacity.
    
    Calculate the isobaric (constant pressure) heat capacity of water
    vapour using the IAPWS-97 formulation.
    
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :returns: Heat capacity in J/kg/K.
    """
    g_tt = vap_g(2,0,temp,pres)
    cp = -temp * g_tt
    return cp

def vap_density(temp,pres):
    """Calculate water vapour density.
    
    Calculate the density of water vapour using the IAPWS-97
    formulation.
    
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :returns: Density in kg/m3.
    """
    vvap = vap_g(0,1,temp,pres)
    dvap = vvap**(-1)
    return dvap

def vap_enthalpy(temp,pres):
    """Calculate water vapour enthalpy.
    
    Calculate the specific enthalpy of water vapour using the IAPWS-97
    formulation.
    
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :returns: Enthalpy in J/kg.
    """
    g = vap_g(0,0,temp,pres)
    g_t = vap_g(1,0,temp,pres)
    h = g - temp*g_t
    return h

def vap_entropy(temp,pres):
    """Calculate water vapour entropy.
    
    Calculate the specific entropy of water vapour using the IAPWS-97
    formulation.
    
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :returns: Entropy in J/kg/K.
    """
    g_t = vap_g(1,0,temp,pres)
    s = -g_t
    return s

def vap_internalenergy(temp,pres):
    """Calculate water vapour internal energy.
    
    Calculate the specific internal energy of water vapour using the
    IAPWS-97 formulation.
    
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :returns: Internal energy in J/kg.
    """
    g = vap_g(0,0,temp,pres)
    g_t = vap_g(1,0,temp,pres)
    g_p = vap_g(0,1,temp,pres)
    u = g - temp*g_t - pres*g_p
    return u

def vap_soundspeed(temp,pres):
    """Calculate water vapour sound speed.
    
    Calculate the speed of sound in water vapour using the IAPWS-97
    formulation.
    
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :returns: Sound speed in m/s.
    """
    g_p = vap_g(0,1,temp,pres)
    g_tt = vap_g(2,0,temp,pres)
    g_tp = vap_g(1,1,temp,pres)
    g_pp = vap_g(0,2,temp,pres)
    csqinv = (g_tp**2/g_tt - g_pp) / g_p**2
    c = csqinv**(-.5)
    return c

def vap_volume(temp,pres):
    """Calculate water vapour specific volume.
    
    Calculate the specific volume of water vapour using the IAPWS-97
    formulation.
    
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :returns: Specific volume in m3/kg.
    """
    v = vap_g(0,1,temp,pres)
    return v


## Check tables
def chkiapws97table(number,printresult=True,chktol=_CHKTOL):
    """Check accuracy against IAPWS-97.
    
    Evaluate the functions in this module and compare to reference
    values from IAPWS 1997. This function is simply a wrapper for the
    chkiapws97table* functions.
    
    :arg int number: Either 5 for liquid water or 15 for water vapour.
    :arg bool printresult: If True (default) and any results are outside
        of the given tolerance, then the function name, reference value,
        result value, and relative error are printed.
    :arg float chktol: Tolerance to use when choosing to print results
        (default _CHKTOL).
    :returns: :class:`~tester.Tester` instances containing the
        functions, arguments, reference values, results, and relative
        errors from the tests. See the individual functions to find what
        each instance corresponds to.
    :raises ValueError: If `number` is not 5 or 15.
    """
    
    if number == 5:
        chkdicts = chk_iapws97_table5(printresult=printresult,tol=tol)
    elif number == 15:
        chkdicts = chk_iapws97_table15(printresult=printresult,tol=tol)
    else:
        errmsg = 'Only tables 5 and 15 are available, not {0}'.format(number)
        raise ValueError(errmsg)
    
    return chkdicts

def chkiapws97table5(printresult=True,chktol=_CHKTOL):
    """Check liquid water accuracy against IAPWS-97.
    
    Evaluate the functions in this module and compare to reference
    values from IAPWS 1997, table 5 (liquid water properties).
    
    :arg bool printresult: If True (default) and any results are outside
        of the given tolerance, then the function name, reference value,
        result value, and relative error are printed.
    :arg float chktol: Tolerance to use when choosing to print results
        (default _CHKTOL).
    :returns: :class:`~tester.Tester` instance containing the
        functions, arguments, reference values, results, and relative
        errors from the tests.
    """
    from tester import Tester
    funs = [liq_cp,liq_enthalpy,liq_entropy,liq_internalenergy,liq_soundspeed,
        liq_volume]
    fargs = [(300.,3e6), (300.,80e6), (500.,3e6)]
    refs_scaled = [
        [0.417301218e1,0.401008987e1,0.465580682e1],
        [0.115331273e3,0.184142828e3,0.975542239e3],
        [0.392294792,0.368563852,0.258041912e1],
        [0.112324818e3,0.106448356e3,0.971934985e3],
        [0.150773921e4,0.163469054e4,0.124071337e4],
        [0.100215168e-2,0.971180894e-3,0.120241800e-2]
    ]
    scales = [1e3,1e3,1e3,1e3,1.,1.]
    refs = [[r*scale for r in ref] for (ref,scale) in zip(refs_scaled,scales)]
    fnames = ['cp','enthalpy','entropy','internalenergy','soundspeed','volume']
    argfmt = '({0:3g},{1:3g})'
    header = 'IF97 liquid functions'
    testliq = Tester(funs,fargs,refs,fnames,argfmt,header=header)
    
    # Run Tester instances and print results
    testliq.run()
    if printresult:
        testliq.printresults(chktol=chktol)
    return testliq

def chkiapws97table15(printresult=True,chktol=_CHKTOL):
    """Check water vapour accuracy against IAPWS-97.
    
    Evaluate the functions in this module and compare to reference
    values from IAPWS 1997, table 5 (water vapour properties).
    
    :arg bool printresult: If True (default) and any results are outside
        of the given tolerance, then the function name, reference value,
        result value, and relative error are printed.
    :arg float chktol: Tolerance to use when choosing to print results
        (default _CHKTOL).
    :returns: :class:`~tester.Tester` instance containing the
        functions, arguments, reference values, results, and relative
        errors from the tests.
    """
    from tester import Tester
    funs = [vap_cp,vap_enthalpy,vap_entropy,vap_internalenergy,vap_soundspeed,
        vap_volume]
    fargs = [(300.,3.5e3), (700.,3.5e3), (700.,30e6)]
    refs_scaled = [
        [0.191300162e1,0.208141274e1,0.103505092e2],
        [0.254991145e4,0.333568375e4,0.263149474e4],
        [0.852238967e1,0.101749996e2,0.517540298e1],
        [0.241169160e4,0.301262819e4,0.246861076e4],
        [0.427920172e3,0.644289068e3,0.480386523e3],
        [0.394913866e2,0.923015898e2,0.542946619e-2],
    ]
    scales = [1e3,1e3,1e3,1e3,1.,1.]
    refs = [[r*scale for r in ref] for (ref,scale) in zip(refs_scaled,scales)]
    fnames = ['cp','enthalpy','entropy','internalenergy','soundspeed','volume']
    argfmt = '({0:3g},{1:5g})'
    header = 'IF97 vapour functions'
    testvap = Tester(funs,fargs,refs,fnames,argfmt,header=header)
    
    # Run Tester instances and print results
    testvap.run()
    if printresult:
        testvap.printresults(chktol=chktol)
    return testvap


## Main function: Check tables
if __name__ == '__main__':
    testliq = chkiapws97table5()
    testvap = chkiapws97table15()

