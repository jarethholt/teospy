"""Gibbs SeaWater library.

This module implements the Gibbs Seawater (GSW) library, thermodynamic
functions for seawater based on a polynomial expression for the Gibbs
free energy of liquid water (from :mod:`liq5_f03`).

Unlike most other modules, the inputs to most functions here are not in
SI units. The salinity is absolute salinity in g/kg; the temperature is
in degrees Celsius; and the pressure is the seawater (gauge) pressure in
decibars (absolute pressure minus 1 atm).

:Examples:

>>> gsw_g(0,0,0,35.,26.85,0.)
-5113.70064124
>>> gsw_g(0,0,1,35.,26.85,0.)
0.977858058750e-3
>>> gsw_g(0,1,1,35.,26.85,0.)
0.304607508052e-6
>>> cp(35.7,25.5,1023.)
3974.42541260
>>> density(35.7,25.5,1023.)
1027.95249316
>>> pottemp(35.7,25.5,1023.,0.)
25.2720983155
>>> potdens(35.7,25.5,1023.,0.)
1023.66254941
>>> alpha_t(35.7,25.5,1023.)
3.09837839319e-4
>>> beta_t(35.7,25.5,1023.)
7.25729797839e-4
>>> tconfromtpot(35.7,25.5)
25.4805463842
>>> cabb_tcon(35.,20.,1e3)
8.96907383083e-6
>>> thrmb_tcon(35.,20.,1e3)
1.72708365652e-8
>>> asalfrompsal(35.,180.,40.,2e3)
35.1890932889958041
>>> psalfromasal(35.7,201.,-21.,1023.)
35.527515065427778

:Functions:

* :func:`asalfrompsal`: Calculate absolute salinity from practical
  salinity.
* :func:`psalfromasal`: Calculate practical salinity from absolute
  salinity.
* :func:`gsw_g`: Seawater Gibbs free energy with derivatives.
* :func:`alpha_t`: Seawater thermal expansion coefficient with respect
  to in-situ temperature.
* :func:`beta_t`: Seawater haline contraction coefficient at constant
  in-situ temperature.
* :func:`cp`: Seawater isobaric heat capacity.
* :func:`density`: Seawater density.
* :func:`enthalpy`: Seawater specific enthalpy.
* :func:`entropy`: Seawater specific entropy.
* :func:`kappa`: Seawater (isentropic) compressibility.
* :func:`kappa_t`: Seawater isothermal compressibility.
* :func:`specvol`: Seawater specific volume.
* :func:`svel`: Speed of sound in seawater.
* :func:`pottemp`: Seawater potential temperature.
* :func:`potdens`: Seawater potential density.
* :func:`alpha_tpot`: Seawater thermal expansion coefficient with
  respect to potential temperature.
* :func:`beta_tpot`: Seawater haline contraction coefficient at constant
  potential temperature.
* :func:`cabb_tpot`: Seawater cabbeling coefficient with respect to
  potential temperature.
* :func:`cabb_tpot_alt`: Seawater cabbeling coefficient with respect to
  potential temperature. This version uses finite differences rather
  than exact derivatives.
* :func:`thrmb_tpot`: Seawater thermobaric coefficient with respect to
  potential temperature.
* :func:`thrmb_tpot_alt`: Seawater thermobaric coefficient with respect
  to potential temperature. This version uses finite differences rather
  than exact derivatives.
* :func:`tconfromtpot`: Calculate seawater conservative temperature from
  potential temperature.
* :func:`tpotfromtcon`: Calculate seawater potential temperature from
  conservative temperature.
* :func:`alpha_tcon`: Seawater thermal expansion coefficient with
  respect to conservative temperature.
* :func:`beta_tcon`: Seawater haline contraction coefficient at constant
  conservative temperature.
* :func:`cabb_tcon`: Seawater cabbeling coefficient with respect to
  conservative temperature.
* :func:`cabb_tcon_alt`: Seawater cabbeling coefficient with respect to
  conservative temperature. This version uses finite differences rather
  than exact derivatives.
* :func:`thrmb_tcon`: Seawater thermobaric coefficient with respect to
  conservative temperature.
* :func:`thrmb_tcon_alt`: Seawater thermobaric coefficient with respect
  to conservative temperature. This version uses finite differences
  rather than exact derivatives.

"""

__all__ = ['asalfrompsal','psalfromasal',
    'gsw_g','alpha_t','beta_t','cp','density','enthalpy','entropy','kappa',
    'kappa_t','specvol','svel',
    'pottemp','potdens','alpha_tpot','beta_tpot','cabb_tpot','cabb_tpot_alt',
    'thrmb_tpot','thrmb_tpot_alt',
    'tconfromtpot','tpotfromtcon','alpha_tcon','beta_tcon','cabb_tcon',
    'cabb_tcon_alt','thrmb_tcon','thrmb_tcon_alt']

import warnings
import numpy
import constants0
import convert0
import sal2
import maths3
import liq5_f03

_G2KG = 1e-3  # Conversion factor from g/kg to kg/kg
_DBAR2PA = 1e4  # Conversion factor from dbar to Pa
_CHKTOL = constants0.CHKTOL
_TCELS = constants0.TCELS
_PATM = constants0.PATM
_SAL0 = constants0.SAL0
_SAL1 = constants0.SAL1
_CSEA = constants0.CSEA
_DTEMP = 1e-3  # Default temperature step size (cabbeling/thermobaric)
_DSALT = 1e-2  # Default salinity step size
_DPRES = 1e-1  # Default pressure step size
_C_TC = (-1.446013646344788e-2, -3.305308995852924e-3, 1.062415929128982e-4,
    9.477566673794488e-1, 2.166591947736613e-3, 3.828842955039902e-3,
    1., 6.506097115635800e-4, 3.830289486850898e-3, 1.247811760368034e-6)
_NMAX = 5
_C_TP = (8.65483913395442e-6,-1.41636299744881e-6,-7.38286467135737e-9,
    -8.38241357039698e-6,2.83933368585534e-8,1.77803965218656e-8,
    1.71155619208233e-10)
_newton = maths3.newton
_liq_g = liq5_f03.liq_g


## Salinity conversion routines
def asalfrompsal(s_psu,lon0,lat0,p_dbar):
    """Convert practical salinity to absolute salinity.
    
    Convert salinity in practical salinity units (PSU) to absolute
    salinity.
    
    :arg float s_psu: Salinity in practical salinity units (psu).
    :arg float lon0: Reference longitude in degrees East.
    :arg float lat0: Reference latitude in degrees North.
    :arg float p_dbar: Seawater (gauge) pressure in decibar.
    :returns: Absolute salinity in g/kg.
    :raises RuntimeWarning: If the reference point is not over the
        ocean, as determined by the GSW_FNAME gridded data file.
        GSW_ERRVAL is returned.
    
    :Examples:
    
    >>> asalfrompsal(35.527515065427778,201.,-21.,1023.)
    35.7
    >>> asalfrompsal(35.,180.,40.,2e3)
    35.1890932889958041
    >>> asalfrompsal(8.,20.,57.,0.)
    8.10483771428571406
    """
    pres = _PATM + p_dbar*_DBAR2PA
    asal = convert0.sal_asalfrompsal(s_psu,lon0,lat0,pres)
    s_gkg = asal / _G2KG
    return s_gkg

def psalfromasal(s_gkg,lon0,lat0,p_dbar):
    """Convert absolute salinity to practical salinity.
    
    Convert salinity in practical salinity units (PSU) to absolute
    salinity.
    
    :arg float s_gkg: Absolute salinity in g/kg.
    :arg float lon0: Reference longitude in degrees East.
    :arg float lat0: Reference latitude in degrees North.
    :arg float p_dbar: Seawater (gauge) pressure in decibar.
    :returns: Salinity in practical salinity units (psu).
    :raises RuntimeWarning: If the reference point is not over the
        ocean, as determined by the GSW_FNAME gridded data file.
        GSW_ERRVAL is returned.
    
    :Examples:
    
    >>> psalfromasal(35.7,201.,-21.,1023.)
    35.527515065427778
    """
    pres = _PATM + p_dbar*_DBAR2PA
    salt = s_gkg * _G2KG
    s_psu = convert0.sal_psalfromasal(salt,lon0,lat0,pres)
    return s_psu


## Gibbs function
def gsw_g(drvs,drvt,drvp,s_gkg,t_cels,p_dbar,chkbnd=False,useext=False):
    """Calculate seawater Gibbs energy with derivatives.
    
    Calculate the specific Gibbs free energy of seawater or its
    derivatives with respect to salinity, temperature, and pressure.
    
    :arg int drvs: Number of salinity derivatives.
    :arg int drvt: Number of temperature derivatives.
    :arg int drvp: Number of pressure derivatives.
    :arg float s_gkg: Absolute salinity in g/kg.
    :arg float t_cels: Temperature in degrees Celsius.
    :arg float p_dbar: Seawater (gauge) pressure in decibar.
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg bool useext: If False (default) then the salt contribution is
        calculated from _GSCOEFFS; if True, from _GSCOEFFS_EXT.
    :returns: Gibbs free energy in units of
        (J/kg) / (g/kg)^drvs / (deg C)^drvt / Pa^drvp.
    
    :Examples:
    
    >>> gsw_g(0,0,0,35.,26.85,0.)
    -5113.70064124
    >>> gsw_g(1,0,0,35.,26.85,0.)
    78.5928261339
    >>> gsw_g(0,1,0,35.,26.85,0.)
    -374.452000830
    >>> gsw_g(0,0,1,35.,26.85,0.)
    0.977858058750e-3
    >>> gsw_g(2,0,0,35.,26.85,0.)
    2.24755137017
    >>> gsw_g(1,1,0,35.,26.85,0.)
    0.789935187192
    >>> gsw_g(1,0,1,35.,26.85,0.)
    -0.716680931996e-6
    >>> gsw_g(0,2,0,35.,26.85,0.)
    -13.3358337534
    >>> gsw_g(0,1,1,35.,26.85,0.)
    0.304607508052e-6
    >>> gsw_g(0,0,2,35.,26.85,0.)
    -0.410939723950e-12
    """
    salt = s_gkg * _G2KG
    temp = t_cels + _TCELS
    pres = p_dbar*_DBAR2PA + _PATM
    if drvs == 0:
        gw = _liq_g(drvt,drvp,temp,pres)
    else:
        gw = 0.
    gs = sal2.sal_g(drvs,drvt,drvp,salt,temp,pres,chkbnd=chkbnd,useext=useext)
    g = gw + gs
    g_gsw = g * (_G2KG)**drvs
    return g_gsw


## Primary thermodynamic properties
def alpha_t(s_gkg,t_cels,p_dbar,chkbnd=False,useext=False):
    """Calculate seawater thermal expansion coefficient.
    
    Calculate the thermal expansion coefficient of seawater with respect
    to in-situ temperature.
    
    :arg float s_gkg: Absolute salinity in g/kg.
    :arg float t_cels: Temperature in degrees Celsius.
    :arg float p_dbar: Seawater (gauge) pressure in decibar.
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg bool useext: If False (default) then the salt contribution is
        calculated from _GSCOEFFS; if True, from _GSCOEFFS_EXT.
    
    Returns:
        alpha (float): Expansion coefficient in 1/(deg C).
    
    Examples:
        >>> gsw_alpha_t(35.7,25.5,1023.)
        3.09837839319e-4
    """
    
    g_p = gsw_g(0,0,1,s_gkg,t_cels,p_dbar,chkbnd=chkbnd,useext=useext)
    g_tp = gsw_g(0,1,1,s_gkg,t_cels,p_dbar,chkbnd=chkbnd,useext=useext)
    alpha = g_tp / g_p
    return alpha

def beta_t(s_gkg,t_cels,p_dbar,chkbnd=False,useext=False):
    """Calculate seawater haline contraction coefficient.
    
    Calculate the haline contraction coefficient of seawater at constant
    in-situ temperature.
    
    :arg float s_gkg: Absolute salinity in g/kg.
    :arg float t_cels: Temperature in degrees Celsius.
    :arg float p_dbar: Seawater (gauge) pressure in decibar.
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg bool useext: If False (default) then the salt contribution is
        calculated from _GSCOEFFS; if True, from _GSCOEFFS_EXT.
    :returns: Contraction coefficient in 1/(g/kg).
    
    :Examples:
    
    >>> beta_t(35.7,25.5,1023.)
    7.25729797839e-4
    """
    g_p = gsw_g(0,0,1,s_gkg,t_cels,p_dbar,chkbnd=chkbnd,useext=useext)
    g_sp = gsw_g(1,0,1,s_gkg,t_cels,p_dbar,chkbnd=chkbnd,useext=useext)
    beta = -g_sp / g_p
    return beta

def cp(s_gkg,t_cels,p_dbar,chkbnd=False,useext=False):
    """Calculate seawater isobaric heat capacity.
    
    Calculate the isobaric heat capacity of seawater.
    
    :arg float s_gkg: Absolute salinity in g/kg.
    :arg float t_cels: Temperature in degrees Celsius.
    :arg float p_dbar: Seawater (gauge) pressure in decibar.
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg bool useext: If False (default) then the salt contribution is
        calculated from _GSCOEFFS; if True, from _GSCOEFFS_EXT.
    :returns: Heat capacity in J/kg/(deg C).
    
    :Examples:
    
    >>> cp(35.7,25.5,1023.)
    3974.42541260
    """
    g_tt = gsw_g(0,2,0,s_gkg,t_cels,p_dbar,chkbnd=chkbnd,useext=useext)
    temp = t_cels + _TCELS
    cp = -temp * g_tt
    return cp

def density(s_gkg,t_cels,p_dbar,chkbnd=False,useext=False):
    """Calculate seawater density.
    
    Calculate the density of seawater.
    
    :arg float s_gkg: Absolute salinity in g/kg.
    :arg float t_cels: Temperature in degrees Celsius.
    :arg float p_dbar: Seawater (gauge) pressure in decibar.
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg bool useext: If False (default) then the salt contribution is
        calculated from _GSCOEFFS; if True, from _GSCOEFFS_EXT.
    :returns: Density in kg/m3.
    
    :Examples:
    
    >>> density(35.7,25.5,1023.)
    1027.95249316
    """
    g_p = gsw_g(0,0,1,s_gkg,t_cels,p_dbar,chkbnd=chkbnd,useext=useext)
    rho = g_p**(-1)
    return rho

def enthalpy(s_gkg,t_cels,p_dbar,chkbnd=False,useext=False):
    """Calculate seawater enthalpy.
    
    Calculate the specific enthalpy of seawater.
    
    :arg float s_gkg: Absolute salinity in g/kg.
    :arg float t_cels: Temperature in degrees Celsius.
    :arg float p_dbar: Seawater (gauge) pressure in decibar.
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg bool useext: If False (default) then the salt contribution is
        calculated from _GSCOEFFS; if True, from _GSCOEFFS_EXT.
    :returns: Enthalpy in J/kg.
    
    :Examples:
    
    >>> enthalpy(35.7,25.5,1023.)
    110776.712409
    """
    g = gsw_g(0,0,0,s_gkg,t_cels,p_dbar,chkbnd=chkbnd,useext=useext)
    g_t = gsw_g(0,1,0,s_gkg,t_cels,p_dbar,chkbnd=chkbnd,useext=useext)
    temp = t_cels + _TCELS
    h = g - temp*g_t
    return h

def entropy(s_gkg,t_cels,p_dbar,chkbnd=False,useext=False):
    """Calculate seawater entropy.
    
    Calculate the specific entropy of seawater.
    
    :arg float s_gkg: Absolute salinity in g/kg.
    :arg float t_cels: Temperature in degrees Celsius.
    :arg float p_dbar: Seawater (gauge) pressure in decibar.
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg bool useext: If False (default) then the salt contribution is
        calculated from _GSCOEFFS; if True, from _GSCOEFFS_EXT.
    :returns: Entropy in J/kg/(deg C).
    
    :Examples:
    
    >>> entropy(35.7,25.5,1023.)
    352.818797715
    """
    g_t = gsw_g(0,1,0,s_gkg,t_cels,p_dbar,chkbnd=chkbnd,useext=useext)
    s = -g_t
    return s

def kappa(s_gkg,t_cels,p_dbar,chkbnd=False,useext=False):
    """Calculate seawater isentropic compressibility.
    
    Calculate the isentropic compressibility of seawater.
    
    :arg float s_gkg: Absolute salinity in g/kg.
    :arg float t_cels: Temperature in degrees Celsius.
    :arg float p_dbar: Seawater (gauge) pressure in decibar.
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg bool useext: If False (default) then the salt contribution is
        calculated from _GSCOEFFS; if True, from _GSCOEFFS_EXT.
    :returns: Compressibility in 1/dbar.
    
    :Examples:
    
    >>> kappa(35.7,25.5,1023.)
    4.03386268546e-6
    """
    g_p = gsw_g(0,0,1,s_gkg,t_cels,p_dbar,chkbnd=chkbnd,useext=useext)
    g_tt = gsw_g(0,2,0,s_gkg,t_cels,p_dbar,chkbnd=chkbnd,useext=useext)
    g_tp = gsw_g(0,1,1,s_gkg,t_cels,p_dbar,chkbnd=chkbnd,useext=useext)
    g_pp = gsw_g(0,0,2,s_gkg,t_cels,p_dbar,chkbnd=chkbnd,useext=useext)
    kappa = (g_tp**2/g_tt - g_pp) / g_p
    kap_dbar = kappa * _DBAR2PA
    return kap_dbar

def kappa_t(s_gkg,t_cels,p_dbar,chkbnd=False,useext=False):
    """Calculate seawater isothermal compressibility.
    
    Calculate the isothermal compressibility of seawater.
    
    :arg float s_gkg: Absolute salinity in g/kg.
    :arg float t_cels: Temperature in degrees Celsius.
    :arg float p_dbar: Seawater (gauge) pressure in decibar.
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg bool useext: If False (default) then the salt contribution is
        calculated from _GSCOEFFS; if True, from _GSCOEFFS_EXT.
    :returns: Compressibility in 1/dbar.
    
    :Examples:
    
    >>> kappa_t(35.7,25.5,1023.)
    4.10403794615e-6
    """
    g_p = gsw_g(0,0,1,s_gkg,t_cels,p_dbar,chkbnd=chkbnd,useext=useext)
    g_pp = gsw_g(0,0,2,s_gkg,t_cels,p_dbar,chkbnd=chkbnd,useext=useext)
    kappa = -g_pp / g_p
    kap_dbar = kappa * _DBAR2PA
    return kap_dbar

def specvol(s_gkg,t_cels,p_dbar,chkbnd=False,useext=False):
    """Calculate seawater specific volume.
    
    Calculate the specific volume of seawater.
    
    :arg float s_gkg: Absolute salinity in g/kg.
    :arg float t_cels: Temperature in degrees Celsius.
    :arg float p_dbar: Seawater (gauge) pressure in decibar.
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg bool useext: If False (default) then the salt contribution is
        calculated from _GSCOEFFS; if True, from _GSCOEFFS_EXT.
    :returns: Specific volume in m3/kg.
    
    :Examples:
    
    >>> specvol(35.7,25.5,1023.)
    9.72807602158e-4
    """
    g_p = gsw_g(0,0,1,s_gkg,t_cels,p_dbar,chkbnd=chkbnd,useext=useext)
    vol = g_p
    return vol

def svel(s_gkg,t_cels,p_dbar,chkbnd=False,useext=False):
    """Calculate seawater sound speed.
    
    Calculate the speed of sound in seawater.
    
    :arg float s_gkg: Absolute salinity in g/kg.
    :arg float t_cels: Temperature in degrees Celsius.
    :arg float p_dbar: Seawater (gauge) pressure in decibar.
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg bool useext: If False (default) then the salt contribution is
        calculated from _GSCOEFFS; if True, from _GSCOEFFS_EXT.
    :returns: Sound speed in m/s.
    
    :Examples:
    
    >>> svel(35.7,25.5,1023.)
    1552.93372863
    """
    g_p = gsw_g(0,0,1,s_gkg,t_cels,p_dbar,chkbnd=chkbnd,useext=useext)
    g_tt = gsw_g(0,2,0,s_gkg,t_cels,p_dbar,chkbnd=chkbnd,useext=useext)
    g_tp = gsw_g(0,1,1,s_gkg,t_cels,p_dbar,chkbnd=chkbnd,useext=useext)
    g_pp = gsw_g(0,0,2,s_gkg,t_cels,p_dbar,chkbnd=chkbnd,useext=useext)
    csqinv = g_tp**2/g_tt - g_pp
    c = g_p * csqinv**(-.5)
    return c


## Potential temperature functions
def _approx_stppp(s_gkg,t_cels,p_dbar,ppot_dbar):
    """Approximate seawater potential temperature.
    
    Approximate the potential temperature of seawater from the salinity,
    in-situ temperature, in-situ pressure, and potential pressure.
    
    :arg float s_gkg: Absolute salinity in g/kg.
    :arg float t_cels: In-situ temperature in degrees Celsius.
    :arg float p_dbar: In-situ gauge pressure in decibar.
    :arg float ppot_dbar: Potential gauge pressure in decibar.
    :returns: Potential temperature in degrees Celsius.
    """
    s_psu = s_gkg * _SAL0/_SAL1
    pminus = p_dbar - ppot_dbar
    pplus = p_dbar + ppot_dbar
    tpot = (t_cels + _C_TP[0]*pminus + _C_TP[1]*pminus*s_psu
        + _C_TP[2]*pminus*pplus + _C_TP[3]*pminus*t_cels
        + _C_TP[4]*pminus*t_cels*s_psu + _C_TP[5]*pminus*t_cels**2
        + _C_TP[6]*pminus*t_cels*pplus)
    return tpot

def _diff_stppp(tp,s_gkg,t_cels,p_dbar,ppot_dbar,chkbnd=False,
    useext=False):
    """Calculate seawater disequilibrium at STPPp.
    
    Calculate both sides of the equation
    
        in-situ entropy = potential entropy
    
    and its derivative with respect to potential temperature. Solving
    this equation gives the potential temperature of seawater at the
    given salinity, in-situ temperature, in-situ pressure, and potential
    pressure.
    
    :arg float tp: Potential temperature in degrees Celsius.
    :arg float s_gkg: Absolute salinity in g/kg.
    :arg float t_cels: In-situ temperature in degrees Celsius.
    :arg float p_dbar: In-situ gauge pressure in decibar.
    :arg float ppot_dbar: Potential gauge pressure in decibar.
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg bool useext: If False (default) then the salt contribution is
        calculated from _GSCOEFFS; if True, from _GSCOEFFS_EXT.
    :returns: Left-hand side of the equation, right-hand side,
        derivative of LHS, and derivative of RHS.
    :rtype: tuple(float)
    """
    s1 = -gsw_g(0,1,0,s_gkg,t_cels,p_dbar,chkbnd=chkbnd,useext=useext)
    s2 = -gsw_g(0,1,0,s_gkg,tp,ppot_dbar,chkbnd=chkbnd,useext=useext)
    lhs = s1
    rhs = s2
    
    s2_t = -gsw_g(0,2,0,s_gkg,tp,ppot_dbar,chkbnd=chkbnd,useext=useext)
    dlhs = 0.
    drhs = s2_t
    return lhs, rhs, dlhs, drhs

def pottemp(s_gkg,t_cels,p_dbar,ppot_dbar,tp_cels=None,chkvals=False,
    chktol=_CHKTOL,tp_cels0=None,chkbnd=False,useext=False,
    mathargs=None):
    """Calculate seawater potential temperature.
    
    Calculate the potential temperature of seawater.
    
    If the calculation has already been done, the results can be passed
    to avoid unnecessary repeat calculations. If enough values are
    passed, they will be checked for consistency if chkvals is True.
    
    :arg float s_gkg: Absolute salinity in g/kg.
    :arg float t_cels: In-situ temperature in degrees Celsius.
    :arg float p_dbar: In-situ gauge pressure in decibar.
    :arg float ppot_dbar: Potential gauge pressure in decibar.
    :arg tp_cels: Potential temperature in degrees Celsius. If unknown,
        pass None (default) and it will be calculated.
    :type tp_cels: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg tp_cels0: Initial guess for the potential temperature in
        degrees Celsius. If None (default) then `_approx_stppp` is used.
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg bool useext: If False (default) then the salt contribution is
        calculated from _GSCOEFFS; if True, from _GSCOEFFS_EXT.
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Potential temperature in degrees Celsius.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    """
    if tp_cels is None:
        fargs = (s_gkg,t_cels,p_dbar,ppot_dbar)
        fkwargs = {'useext': useext}
        if mathargs is None:
            mathargs = dict()
        tp_cels = _newton(_diff_stppp,tp_cels0,_approx_stppp,fargs=fargs,
            fkwargs=fkwargs,**mathargs)
    
    if not chkvals:
        return tp_cels
    
    lhs, rhs, __, __ = _diff_stppp(tp_cels,s_gkg,t_cels,p_dbar,ppot_dbar,
        useext=useext)
    errs = list()
    for (l,r) in zip(lhs,rhs):
        if abs(r) >= chktol:
            errs.append(abs(l/r-1))
        else:
            errs.append(abs(l-r))
    if max(errs) > chktol:
        warnmsg = ('Given values {0} and solutions {1} disagree to more than '
            'the tolerance {2}').format(lhs,rhs,chktol)
        warnings.warn(warnmsg,RuntimeWarning)
    return tp_cels

def potdens(s_gkg,t_cels,p_dbar,ppot_dbar,tp_cels=None,chkvals=False,
    chktol=_CHKTOL,tp_cels0=None,chkbnd=False,useext=False,
    mathargs=None):
    """Calculate seawater potential density.
    
    Calculate the potential density of seawater.
    
    :arg float s_gkg: Absolute salinity in g/kg.
    :arg float t_cels: In-situ temperature in degrees Celsius.
    :arg float p_dbar: In-situ gauge pressure in decibar.
    :arg float ppot_dbar: Potential gauge pressure in decibar.
    :arg tp_cels: Potential temperature in degrees Celsius. If unknown,
        pass None (default) and it will be calculated.
    :type tp_cels: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg tp_cels0: Initial guess for the potential temperature in
        degrees Celsius. If None (default) then `_approx_stppp` is used.
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg bool useext: If False (default) then the salt contribution is
        calculated from _GSCOEFFS; if True, from _GSCOEFFS_EXT.
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Potential density in kg/m3.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    """
    tp_cels = pottemp(s_gkg,t_cels,p_dbar,ppot_dbar,tp_cels=tp_cels,
        chkvals=chkvals,chktol=chktol,tp_cels0=tp_cels0,chkbnd=chkbnd,
        useext=useext,mathargs=mathargs)
    dpot = density(s_gkg,tp_cels,ppot_dbar,useext=useext)
    return dpot

def alpha_tpot(s_gkg,t_cels,p_dbar,tp_cels=None,chkvals=False,
    chktol=_CHKTOL,tp_cels0=None,chkbnd=False,useext=False,
    mathargs=None):
    """Calculate thermal expansion wrt potential temperature.
    
    Calculate the thermal expansion coefficient of seawater with respect
    to potential temperature.
    
    :arg float s_gkg: Absolute salinity in g/kg.
    :arg float t_cels: In-situ temperature in degrees Celsius.
    :arg float p_dbar: In-situ gauge pressure in decibar.
    :arg tp_cels: Potential temperature in degrees Celsius. If unknown,
        pass None (default) and it will be calculated.
    :type tp_cels: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg tp_cels0: Initial guess for the potential temperature in
        degrees Celsius. If None (default) then `_approx_stppp` is used.
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg bool useext: If False (default) then the salt contribution is
        calculated from _GSCOEFFS; if True, from _GSCOEFFS_EXT.
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Expansion coefficient in 1/(deg C).
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> alpha_tpot(35.,20.,1e3)
    2.69753733317e-04
    """
    tp_cels = pottemp(s_gkg,t_cels,p_dbar,0.,tp_cels=tp_cels,chkvals=chkvals,
        chktol=chktol,tp_cels0=tp_cels0,chkbnd=chkbnd,useext=useext,
        mathargs=mathargs)
    g_p = gsw_g(0,0,1,s_gkg,t_cels,p_dbar,useext=useext)
    g_tt = gsw_g(0,2,0,s_gkg,t_cels,p_dbar,useext=useext)
    g_tp = gsw_g(0,1,1,s_gkg,t_cels,p_dbar,useext=useext)
    h_p = g_p
    h_ep = g_tp / (-g_tt)
    g0_tt = gsw_g(0,2,0,s_gkg,tp_cels,0.,useext=useext)
    h0_ee = (-g0_tt)**(-1)
    alpha = h_ep / (h_p * h0_ee)
    return alpha

def beta_tpot(s_gkg,t_cels,p_dbar,tp_cels=None,chkvals=False,
    chktol=_CHKTOL,tp_cels0=None,chkbnd=False,useext=False,
    mathargs=None):
    """Calculate haline contraction wrt potential temperature.
    
    Calculate the haline contraction coefficient of seawater at constant
    potential temperature.
    
    :arg float s_gkg: Absolute salinity in g/kg.
    :arg float t_cels: In-situ temperature in degrees Celsius.
    :arg float p_dbar: In-situ gauge pressure in decibar.
    :arg tp_cels: Potential temperature in degrees Celsius. If unknown,
        pass None (default) and it will be calculated.
    :type tp_cels: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg tp_cels0: Initial guess for the potential temperature in
        degrees Celsius. If None (default) then `_approx_stppp` is used.
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg bool useext: If False (default) then the salt contribution is
        calculated from _GSCOEFFS; if True, from _GSCOEFFS_EXT.
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Contraction coefficient in 1/(g/kg).
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> beta_tpot(35.,20.,1e3)
    7.31582583383e-4
    """
    tp_cels = pottemp(s_gkg,t_cels,p_dbar,0.,tp_cels=tp_cels,chkvals=chkvals,
        chktol=chktol,tp_cels0=tp_cels0,chkbnd=chkbnd,useext=useext,
        mathargs=mathargs)
    
    g_p = gsw_g(0,0,1,s_gkg,t_cels,p_dbar,useext=useext)
    g_sp = gsw_g(1,0,1,s_gkg,t_cels,p_dbar,useext=useext)
    g_st = gsw_g(1,1,0,s_gkg,t_cels,p_dbar,useext=useext)
    g_tt = gsw_g(0,2,0,s_gkg,t_cels,p_dbar,useext=useext)
    g_tp = gsw_g(0,1,1,s_gkg,t_cels,p_dbar,useext=useext)
    h_p = g_p
    h_ep = g_tp / (-g_tt)
    h_sp = g_sp + g_tp * g_st/(-g_tt)
    g0_s = gsw_g(1,0,0,s_gkg,tp_cels,0.,useext=useext)
    g0_st = gsw_g(1,1,0,s_gkg,tp_cels,0.,useext=useext)
    g0_tt = gsw_g(0,2,0,s_gkg,tp_cels,0.,useext=useext)
    h0_s = g0_s
    h0_e = tp_cels + _TCELS
    h0_se = g0_st / (-g0_tt)
    h0_ee = (-g0_tt)**(-1)
    beta = (h0_se*h_ep - h0_ee*h_sp) / (h_p * h0_ee)
    return beta

def cabb_tpot(s_gkg,t_cels,p_dbar,tp_cels=None,chkvals=False,
    chktol=_CHKTOL,tp_cels0=None,chkbnd=False,useext=False,
    mathargs=None):
    """Calculate cabbeling wrt potential temperature.
    
    Calculate the cabbeling coefficient of seawater with respect to
    potential temperature. This version uses the exact derivatives of
    the polynomial expression for the free energy of seawater.
    
    :arg float s_gkg: Absolute salinity in g/kg.
    :arg float t_cels: In-situ temperature in degrees Celsius.
    :arg float p_dbar: In-situ gauge pressure in decibar.
    :arg tp_cels: Potential temperature in degrees Celsius. If unknown,
        pass None (default) and it will be calculated.
    :type tp_cels: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg tp_cels0: Initial guess for the potential temperature in
        degrees Celsius. If None (default) then `_approx_stppp` is used.
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg bool useext: If False (default) then the salt contribution is
        calculated from _GSCOEFFS; if True, from _GSCOEFFS_EXT.
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Cabbeling coefficient in 1/(deg C)^2.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> cabb_tpot(35.,20.,1e3)
    8.765978783080191e-6
    """
    tp_cels = pottemp(s_gkg,t_cels,p_dbar,0.,tp_cels=tp_cels,chkvals=chkvals,
        chktol=chktol,tp_cels0=tp_cels0,chkbnd=chkbnd,useext=useext,
        mathargs=mathargs)
    
    alpha0 = alpha_tpot(s_gkg,t_cels,p_dbar,tp_cels=tp_cels,useext=useext)
    beta0 = beta_tpot(s_gkg,t_cels,p_dbar,tp_cels=tp_cels,useext=useext)
    ratio = alpha0/beta0
    
    g_p = gsw_g(0,0,1,s_gkg,t_cels,p_dbar,useext=useext)
    g_st = gsw_g(1,1,0,s_gkg,t_cels,p_dbar,useext=useext)
    g_sp = gsw_g(1,0,1,s_gkg,t_cels,p_dbar,useext=useext)
    g_tt = gsw_g(0,2,0,s_gkg,t_cels,p_dbar,useext=useext)
    g_tp = gsw_g(0,1,1,s_gkg,t_cels,p_dbar,useext=useext)
    g_sst = gsw_g(2,1,0,s_gkg,t_cels,p_dbar,useext=useext)
    g_ssp = gsw_g(2,0,1,s_gkg,t_cels,p_dbar,useext=useext)
    g_stt = gsw_g(1,2,0,s_gkg,t_cels,p_dbar,useext=useext)
    g_stp = gsw_g(1,1,1,s_gkg,t_cels,p_dbar,useext=useext)
    g_ttt = gsw_g(0,3,0,s_gkg,t_cels,p_dbar,useext=useext)
    g_ttp = gsw_g(0,2,1,s_gkg,t_cels,p_dbar,useext=useext)
    g0_st = gsw_g(1,1,0,s_gkg,tp_cels,0.,useext=useext)
    g0_tt = gsw_g(0,2,0,s_gkg,tp_cels,0.,useext=useext)
    g0_sst = gsw_g(2,1,0,s_gkg,tp_cels,0.,useext=useext)
    g0_stt = gsw_g(1,2,0,s_gkg,tp_cels,0.,useext=useext)
    g0_ttt = gsw_g(0,3,0,s_gkg,tp_cels,0.,useext=useext)
    
    dp = g_p**(-1)
    dtt = g_tt**(-1)
    t_s = (g0_st - g_st) * dtt
    t_th = g0_tt * dtt
    t_p = -g_tp * dtt
    vol_s = g_sp*dp
    vol_t = g_tp*dp
    cp_s = g_stt*dtt
    cp_t = g_ttt*dtt
    k_ttp = g_ttp - g_tp*(vol_t + cp_t)
    k_stp = g_stp - g_tp*(vol_s + cp_s) + t_s*k_ttp
    
    dalpha_th = dp * (t_th**2*k_ttp - t_p*g0_ttt)
    dalpha_s = dp * (t_th*k_stp - t_p*g0_stt)
    dbeta_s = -dp * (g_ssp - g_sp*vol_s + t_s*(g_stp - g_sp*vol_t + k_stp)
        + t_p*(g_sst + t_s*g_stt - g0_sst))
    cabb = dalpha_th + 2*ratio*dalpha_s - ratio**2*dbeta_s
    return cabb

def cabb_tpot_alt(s_gkg,t_cels,p_dbar,dtemp=_DTEMP,dsalt=_DSALT,
    tp_cels=None,chkvals=False,chktol=_CHKTOL,tp_cels0=None,
    chkbnd=False,useext=False,mathargs=None):
    """Calculate cabbeling coefficient wrt potential temperature.
    
    Calculate the cabbeling coefficient of seawater with respect to
    potential temperature. This version uses finite differences to
    approximate derivatives.
    
    :arg float s_gkg: Absolute salinity in g/kg.
    :arg float t_cels: In-situ temperature in degrees Celsius.
    :arg float p_dbar: In-situ gauge pressure in decibar.
    :arg float dtemp: Amount to change the temperature by when
        calculating finite differences (default _DTEMP).
    :arg float dsalt: Amount to change the salinity by when calculating
        finite differences (default _DSALT).
    :arg tp_cels: Potential temperature in degrees Celsius. If unknown,
        pass None (default) and it will be calculated.
    :type tp_cels: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg tp_cels0: Initial guess for the potential temperature in
        degrees Celsius. If None (default) then `_approx_stppp` is used.
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg bool useext: If False (default) then the salt contribution is
        calculated from _GSCOEFFS; if True, from _GSCOEFFS_EXT.
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Cabbeling coefficient in 1/(deg C)^2.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> cabb_tpot_alt(35.,20.,1e3)
    8.75963154048e-6
    """
    tp_cels = pottemp(s_gkg,t_cels,p_dbar,0.,tp_cels=tp_cels,chkvals=chkvals,
        chktol=chktol,tp_cels0=tp_cels0,chkbnd=chkbnd,useext=useext,
        mathargs=mathargs)
    alpha0 = alpha_tpot(s_gkg,t_cels,p_dbar,tp_cels=tp_cels,useext=useext)
    beta0 = beta_tpot(s_gkg,t_cels,p_dbar,tp_cels=tp_cels,useext=useext)
    ratio = alpha0/beta0
    
    # Calculate finite difference wrt potential temperature
    tp_l = tp_cels - dtemp
    tp_u = tp_cels + dtemp
    t_l = pottemp(s_gkg,tp_l,0.,p_dbar,useext=useext,mathargs=mathargs)
    t_u = pottemp(s_gkg,tp_u,0.,p_dbar,useext=useext,mathargs=mathargs)
    alpha_l = alpha_tpot(s_gkg,t_l,p_dbar,tp_cels=tp_l,useext=useext)
    alpha_u = alpha_tpot(s_gkg,t_u,p_dbar,tp_cels=tp_u,useext=useext)
    dalpha_th = (alpha_u - alpha_l) / (tp_u - tp_l)
    
    # Calculate finite difference wrt salinity
    if s_gkg >= dsalt:
        s_l = s_gkg - dsalt
        s_u = s_gkg + dsalt
    else:
        s_l = 0.
        s_u = dsalt
    t_l = pottemp(s_l,tp_cels,0.,p_dbar,useext=useext,mathargs=mathargs)
    t_u = pottemp(s_u,tp_cels,0.,p_dbar,useext=useext,mathargs=mathargs)
    alpha_l = alpha_tpot(s_l,t_l,p_dbar,tp_cels=tp_cels,useext=useext)
    alpha_u = alpha_tpot(s_u,t_u,p_dbar,tp_cels=tp_cels,useext=useext)
    beta_l = beta_tpot(s_l,t_l,p_dbar,tp_cels=tp_cels,useext=useext)
    beta_u = beta_tpot(s_u,t_u,p_dbar,tp_cels=tp_cels,useext=useext)
    dalpha_s = (alpha_u - alpha_l) / (s_u - s_l)
    dbeta_s = (beta_u - beta_l) / (s_u - s_l)
    cabb = dalpha_th + 2*ratio*dalpha_s - ratio**2*dbeta_s
    return cabb

def thrmb_tpot(s_gkg,t_cels,p_dbar,tp_cels=None,chkvals=False,
    chktol=_CHKTOL,tp_cels0=None,chkbnd=False,useext=False,
    mathargs=None):
    """Calculate thermobaricity wrt potential temperature.
    
    Calculate the thermobaric coefficient of seawater with respect to
    potential temperature. This version uses the exact derivatives of
    the polynomial expression for the free energy of seawater.
    
    :arg float s_gkg: Absolute salinity in g/kg.
    :arg float t_cels: In-situ temperature in degrees Celsius.
    :arg float p_dbar: In-situ gauge pressure in decibar.
    :arg tp_cels: Potential temperature in degrees Celsius. If unknown,
        pass None (default) and it will be calculated.
    :type tp_cels: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg tp_cels0: Initial guess for the potential temperature in
        degrees Celsius. If None (default) then `_approx_stppp` is used.
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg bool useext: If False (default) then the salt contribution is
        calculated from _GSCOEFFS; if True, from _GSCOEFFS_EXT.
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Thermobaric coefficient in 1/dbar/(deg C).
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> thrmb_tpot(35.,20.,1e3)
    1.70965495578822e-08
    """
    tp_cels = pottemp(s_gkg,t_cels,p_dbar,0.,tp_cels=tp_cels,chkvals=chkvals,
        chktol=chktol,tp_cels0=tp_cels0,chkbnd=chkbnd,useext=useext,
        mathargs=mathargs)
    alpha0 = alpha_tpot(s_gkg,t_cels,p_dbar,tp_cels=tp_cels,useext=useext)
    beta0 = beta_tpot(s_gkg,t_cels,p_dbar,tp_cels=tp_cels,useext=useext)
    
    g_p = gsw_g(0,0,1,s_gkg,t_cels,p_dbar,useext=useext)
    g_st = gsw_g(1,1,0,s_gkg,t_cels,p_dbar,useext=useext)
    g_sp = gsw_g(1,0,1,s_gkg,t_cels,p_dbar,useext=useext)
    g_tt = gsw_g(0,2,0,s_gkg,t_cels,p_dbar,useext=useext)
    g_tp = gsw_g(0,1,1,s_gkg,t_cels,p_dbar,useext=useext)
    g_pp = gsw_g(0,0,2,s_gkg,t_cels,p_dbar,useext=useext)
    g_stt = gsw_g(1,2,0,s_gkg,t_cels,p_dbar,useext=useext)
    g_stp = gsw_g(1,1,1,s_gkg,t_cels,p_dbar,useext=useext)
    g_spp = gsw_g(1,0,2,s_gkg,t_cels,p_dbar,useext=useext)
    g_ttt = gsw_g(0,3,0,s_gkg,t_cels,p_dbar,useext=useext)
    g_ttp = gsw_g(0,2,1,s_gkg,t_cels,p_dbar,useext=useext)
    g_tpp = gsw_g(0,1,2,s_gkg,t_cels,p_dbar,useext=useext)
    g0_st = gsw_g(1,1,0,s_gkg,tp_cels,0.,useext=useext)
    g0_tt = gsw_g(0,2,0,s_gkg,tp_cels,0.,useext=useext)

    dp = g_p**(-1)
    dtt = g_tt**(-1)
    dptt = dp*dtt
    t_s = (g0_st - g_st) * dtt
    t_p = -g_tp * dtt
    vol_p = (g_pp + t_p*g_tp) * dp
    cp_p = (g_ttp + t_p*g_ttt) * dtt
    k_tpp = g_tpp + t_p*g_ttp - g_tp*(vol_p + cp_p)
    
    dalpha_p = g0_tt*dptt * k_tpp
    dbeta_p = -dp * (g_spp + 2*t_p*g_stp + t_p**2*g_stt - g_sp*vol_p
        + t_s*k_tpp)
    thermb = dalpha_p - alpha0/beta0*dbeta_p
    thermb *= _DBAR2PA
    return thermb

def thrmb_tpot_alt(s_gkg,t_cels,p_dbar,dpres=_DPRES,tp_cels=None,
    chkvals=False,chktol=_CHKTOL,tp_cels0=None,chkbnd=False,
    useext=False,mathargs=None):
    """Calculate thermobaricity wrt potential temperature.
    
    Calculate the thermobaric coefficient of seawater with respect to
    potential temperature. This version uses finite differences to
    approximate derivatives.
    
    :arg float s_gkg: Absolute salinity in g/kg.
    :arg float t_cels: In-situ temperature in degrees Celsius.
    :arg float p_dbar: In-situ gauge pressure in decibar.
    :arg float dpres: Amount to change the pressure by when calculating
        finite differences (default _DPRES).
    :arg tp_cels: Potential temperature in degrees Celsius. If unknown,
        pass None (default) and it will be calculated.
    :type tp_cels: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg tp_cels0: Initial guess for the potential temperature in
        degrees Celsius. If None (default) then `_approx_stppp` is used.
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg bool useext: If False (default) then the salt contribution is
        calculated from _GSCOEFFS; if True, from _GSCOEFFS_EXT.
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Thermobaric coefficient in 1/dbar/(deg C).
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> thrmb_tpot_alt(35.,20.,1e3)
    1.70945045984e-8
    """
    tp_cels = pottemp(s_gkg,t_cels,p_dbar,0.,tp_cels=tp_cels,chkvals=chkvals,
        chktol=chktol,tp_cels0=tp_cels0,chkbnd=chkbnd,useext=useext,
        mathargs=mathargs)
    beta0 = beta_tpot(s_gkg,t_cels,p_dbar,tp_cels=tp_cels,useext=useext)
    
    # Calculate finite difference wrt pressure
    if p_dbar >= dpres:
        p_l = p_dbar - dpres
        p_u = p_dbar + dpres
    else:
        p_l = 0.
        p_u = dpres
    t_l = pottemp(s_gkg,tp_cels,0.,p_l,useext=useext,mathargs=mathargs)
    t_u = pottemp(s_gkg,tp_cels,0.,p_u,useext=useext,mathargs=mathargs)
    alpha_l = alpha_tpot(s_gkg,t_l,p_l,tp_cels=tp_cels,useext=useext)
    alpha_u = alpha_tpot(s_gkg,t_u,p_u,tp_cels=tp_cels,useext=useext)
    beta_l = beta_tpot(s_gkg,t_l,p_l,tp_cels=tp_cels,useext=useext)
    beta_u = beta_tpot(s_gkg,t_u,p_u,tp_cels=tp_cels,useext=useext)
    thermb = beta0 * (alpha_u/beta_u - alpha_l/beta_l) / (p_u - p_l)
    return thermb


## Conservative temperature functions
def tconfromtpot(s_gkg,tp_cels,chkbnd=False,useext=False):
    """Calculate seawater conservative from potential temperature.
    
    Calculate the conservative temperature of seawater from the
    potential temperature.
    
    :arg float s_gkg: Absolute salinity in g/kg.
    :arg float tp_cels: Potential temperature in degrees Celsius.
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg bool useext: If False (default) then the salt contribution is
        calculated from _GSCOEFFS; if True, from _GSCOEFFS_EXT.
    :returns: Conservative temperature in degrees Celsius.
    
    :Examples:
    
    >>> tconfromtpot(35.7,25.5)
    25.4805463842
    """
    hpot = enthalpy(s_gkg,tp_cels,0.,chkbnd=chkbnd,useext=useext)
    tcon = hpot / _CSEA
    return tcon

def _approx_stc(s_gkg,tc_cels):
    """Approximate Tp from STc.
    
    Approximate the potential temperature of seawater from the salinity
    and conservative temperature.
    
    :arg float s_gkg: Absolute salinity in g/kg.
    :arg float tc_cels: Conservative temperature in degrees Celsius.
    :returns: Potential temperature in degrees Celsius.
    """
    s1 = s_gkg * _SAL0/_SAL1
    th0num = (_C_TC[0] + _C_TC[1]*s1 + _C_TC[2]*s1**2
        + _C_TC[3]*tc_cels + _C_TC[4]*s1*tc_cels + _C_TC[5]*tc_cels**2)
    th0den = _C_TC[6] + _C_TC[7]*s1 + _C_TC[8]*tc_cels + _C_TC[9]*tc_cels**2
    tpot = th0num / th0den
    return tpot

def _diff_stc(tp,s_gkg,tc_cels,chkbnd=False,useext=False):
    """Calculate seawater disequilibrium at STc.
    
    Calculate both sides of the equation
    
        _CSEA*(conservative temperature) = potential enthalpy
    
    and their derivatives with respect to potential temperature. Solving
    this equation gives the potential temperature at the given salinity
    and conservative temperature.
    
    :arg float tp: Potential temperature in degrees Celsius.
    :arg float s_gkg: Absolute salinity in g/kg.
    :arg float tc_cels: Conservative temperature in degrees Celsius.
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg bool useext: If False (default) then the salt contribution is
        calculated from _GSCOEFFS; if True, from _GSCOEFFS_EXT.
    :returns: Left-hand side of the equation, right-hand side,
        derivative of LHS, and derivative of RHS.
    :rtype: tuple(float)
    """
    hp = enthalpy(s_gkg,tp,0.,chkbnd=chkbnd,useext=useext)
    lhs = _CSEA * tc_cels
    rhs = hp
    
    g_tt = gsw_g(0,2,0,s_gkg,tp,0.,chkbnd=chkbnd,useext=useext)
    hp_t = -(_TCELS+tp) * g_tt
    dlhs = 0.
    drhs = hp_t
    return lhs, rhs, dlhs, drhs

def tpotfromtcon(s_gkg,tc_cels,tp_cels=None,chkvals=False,
    chktol=_CHKTOL,tp_cels0=None,chkbnd=False,useext=False,
    mathargs=None):
    """Calculate seawater potential from conservative temperature.
    
    Calculate the potential temperature of seawater from the conservative
    temperature. This is done by solving the equation
    
        _CSEA*(conservative temperature) = potential enthalpy.
    
    If the calculation has already been done, the results can be passed
    to avoid unnecessary repeat calculations. If enough values are
    passed, they will be checked for consistency if chkvals is True.
    
    :arg float s_gkg: Absolute salinity in g/kg.
    :arg float tc_cels: Conservative temperature in degrees Celsius.
    :arg tp_cels: Potential temperature in degrees Celsius. If unknown,
        pass None (default) and it will be calculated.
    :type tp_cels: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg tp_cels0: Initial guess for the potential temperature in
        degrees Celsius. If None (default) then `_approx_stc` is used.
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg bool useext: If False (default) then the salt contribution is
        calculated from _GSCOEFFS; if True, from _GSCOEFFS_EXT.
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Potential temperature in degrees Celsius.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> tpotfromtcon(35.7,25.4805463842239)
    25.5
    """
    if tp_cels is None:
        fargs = (s_gkg,tc_cels)
        fkwargs = {'useext': useext}
        if mathargs is None:
            mathargs = dict()
        tp_cels = _newton(_diff_stc,tp_cels0,_approx_stc,fargs=fargs,
            fkwargs=fkwargs,**mathargs)
    
    if not chkvals:
        return tp_cels
    
    lhs, rhs, __, __ = _diff_stc(tp_cels,s_gkg,tc_cels,useext=useext)
    errs = [abs(lhs/rhs-1)]
    if max(errs) > chktol:
        warnmsg = ('Given value {0} and solution {1} disagree to more than the '
            'tolerance {2}').format(lhs,rhs,chktol)
        warnings.warn(warnmsg,RuntimeWarning)
    return tp_cels

def alpha_tcon(s_gkg,t_cels,p_dbar,tp_cels=None,chkvals=False,
    chktol=_CHKTOL,tp_cels0=None,chkbnd=False,useext=False,
    mathargs=None):
    """Calculate thermal expansion wrt conservative temperature.
    
    Calculate the thermal expansion coefficient of seawater with respect
    to conservative temperature.
    
    :arg float s_gkg: Absolute salinity in g/kg.
    :arg float t_cels: In-situ temperature in degrees Celsius.
    :arg float p_dbar: In-situ gauge pressure in decibar.
    :arg tp_cels: Potential temperature in degrees Celsius. If unknown,
        pass None (default) and it will be calculated.
    :type tp_cels: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg tp_cels0: Initial guess for the potential temperature in
        degrees Celsius. If None (default) then `_approx_stppp` is used.
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg bool useext: If False (default) then the salt contribution is
        calculated from _GSCOEFFS; if True, from _GSCOEFFS_EXT.
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Expansion coefficient in 1/(deg C).
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> alpha_tcon(35.,20.,1e3)
    2.69418609861e-4
    """
    tp_cels = pottemp(s_gkg,t_cels,p_dbar,0.,tp_cels=tp_cels,chkvals=chkvals,
        chktol=chktol,tp_cels0=tp_cels0,chkbnd=chkbnd,useext=useext,
        mathargs=mathargs)
    
    g_p = gsw_g(0,0,1,s_gkg,t_cels,p_dbar,useext=useext)
    g_tt = gsw_g(0,2,0,s_gkg,t_cels,p_dbar,useext=useext)
    g_tp = gsw_g(0,1,1,s_gkg,t_cels,p_dbar,useext=useext)
    h_p = g_p
    h_ep = g_tp / -g_tt
    h0_e = tp_cels + _TCELS
    alpha_h = h_ep / (h_p * h0_e)
    alpha = alpha_h * _CSEA
    return alpha

def beta_tcon(s_gkg,t_cels,p_dbar,tp_cels=None,chkvals=False,
    chktol=_CHKTOL,tp_cels0=None,chkbnd=False,useext=False,
    mathargs=None):
    """Calculate haline contraction wrt conservative temperature.
    
    Calculate the haline contraction coefficient of seawater at constant
    conservative temperature.
    
    :arg float s_gkg: Absolute salinity in g/kg.
    :arg float t_cels: In-situ temperature in degrees Celsius.
    :arg float p_dbar: In-situ gauge pressure in decibar.
    :arg tp_cels: Potential temperature in degrees Celsius. If unknown,
        pass None (default) and it will be calculated.
    :type tp_cels: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg tp_cels0: Initial guess for the potential temperature in
        degrees Celsius. If None (default) then `_approx_stppp` is used.
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg bool useext: If False (default) then the salt contribution is
        calculated from _GSCOEFFS; if True, from _GSCOEFFS_EXT.
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Contraction coefficient in 1/(g/kg).
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> beta_tcon(35.,20.,1e3)
    7.23213672954e-4
    """
    tp_cels = pottemp(s_gkg,t_cels,p_dbar,0.,tp_cels=tp_cels,chkvals=chkvals,
        chktol=chktol,tp_cels0=tp_cels0,chkbnd=chkbnd,useext=useext,
        mathargs=mathargs)
    
    g_p = gsw_g(0,0,1,s_gkg,t_cels,p_dbar,useext=useext)
    g_sp = gsw_g(1,0,1,s_gkg,t_cels,p_dbar,useext=useext)
    g_st = gsw_g(1,1,0,s_gkg,t_cels,p_dbar,useext=useext)
    g_tt = gsw_g(0,2,0,s_gkg,t_cels,p_dbar,useext=useext)
    g_tp = gsw_g(0,1,1,s_gkg,t_cels,p_dbar,useext=useext)
    g0_s = gsw_g(1,0,0,s_gkg,tp_cels,0.,useext=useext)
    h_p = g_p
    h_ep = g_tp / (-g_tt)
    h_sp = g_sp + g_tp * g_st/(-g_tt)
    h0_s = g0_s
    h0_e = tp_cels + _TCELS
    beta = (h0_s*h_ep - h0_e*h_sp) / (h_p*h0_e)
    return beta

def cabb_tcon(s_gkg,t_cels,p_dbar,tp_cels=None,chkvals=False,
    chktol=_CHKTOL,tp_cels0=None,chkbnd=False,useext=False,
    mathargs=None):
    """Calculate cabbeling wrt conservative temperature.
    
    Calculate the cabbeling coefficient of seawater with respect to
    conservative temperature. This version uses the exact derivatives of
    the polynomial expression for the free energy of seawater.
    
    :arg float s_gkg: Absolute salinity in g/kg.
    :arg float t_cels: In-situ temperature in degrees Celsius.
    :arg float p_dbar: In-situ gauge pressure in decibar.
    :arg tp_cels: Potential temperature in degrees Celsius. If unknown,
        pass None (default) and it will be calculated.
    :type tp_cels: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg tp_cels0: Initial guess for the potential temperature in
        degrees Celsius. If None (default) then `_approx_stppp` is used.
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg bool useext: If False (default) then the salt contribution is
        calculated from _GSCOEFFS; if True, from _GSCOEFFS_EXT.
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Cabbeling coefficient in 1/(deg C)^2.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> cabb_tcon(35.,20.,1e3)
    9.163247496098528e-06
    """
    tp_cels = pottemp(s_gkg,t_cels,p_dbar,0.,tp_cels=tp_cels,chkvals=chkvals,
        chktol=chktol,tp_cels0=tp_cels0,chkbnd=chkbnd,useext=useext,
        mathargs=mathargs)
    alpha0 = alpha_tcon(s_gkg,t_cels,p_dbar,tp_cels=tp_cels,useext=useext)
    beta0 = beta_tcon(s_gkg,t_cels,p_dbar,tp_cels=tp_cels,useext=useext)
    ratio = alpha0/beta0
    
    g_p = gsw_g(0,0,1,s_gkg,t_cels,p_dbar,useext=useext)
    g_st = gsw_g(1,1,0,s_gkg,t_cels,p_dbar,useext=useext)
    g_sp = gsw_g(1,0,1,s_gkg,t_cels,p_dbar,useext=useext)
    g_tt = gsw_g(0,2,0,s_gkg,t_cels,p_dbar,useext=useext)
    g_tp = gsw_g(0,1,1,s_gkg,t_cels,p_dbar,useext=useext)
    g_sst = gsw_g(2,1,0,s_gkg,t_cels,p_dbar,useext=useext)
    g_ssp = gsw_g(2,0,1,s_gkg,t_cels,p_dbar,useext=useext)
    g_stt = gsw_g(1,2,0,s_gkg,t_cels,p_dbar,useext=useext)
    g_stp = gsw_g(1,1,1,s_gkg,t_cels,p_dbar,useext=useext)
    g_ttt = gsw_g(0,3,0,s_gkg,t_cels,p_dbar,useext=useext)
    g_ttp = gsw_g(0,2,1,s_gkg,t_cels,p_dbar,useext=useext)
    g0_s = gsw_g(1,0,0,s_gkg,tp_cels,0.,useext=useext)
    g0_ss = gsw_g(2,0,0,s_gkg,tp_cels,0.,useext=useext)
    g0_st = gsw_g(1,1,0,s_gkg,tp_cels,0.,useext=useext)
    g0_tt = gsw_g(0,2,0,s_gkg,tp_cels,0.,useext=useext)
    
    dp = g_p**(-1)
    dtt = g_tt**(-1)
    dptt = dp*dtt
    tpinv = (_TCELS + tp_cels)**(-1)
    t_s = (-g_st + tpinv*g0_s)*dtt
    t_tc = -_CSEA*tpinv*dtt
    t_p = -g_tp*dtt
    th_s = (g0_s*tpinv - g0_st)/g0_tt
    th_tc = -_CSEA*tpinv/g0_tt
    vol_s = g_sp*dp
    vol_t = g_tp*dp
    cp_s = g_stt*dtt
    cp_t = g_ttt*dtt
    k_ttp = g_ttp - g_tp*(vol_t + cp_t)
    k_stp = g_stp - g_tp*(vol_s + cp_s) + t_s*k_ttp
    
    dalpha_tc = dp*t_tc * (t_tc*k_ttp - th_tc*tpinv*g_tp)
    dalpha_s = dp*t_tc * (k_stp - th_s*tpinv*g_tp)
    dbeta_s = -dp * (g_ssp - g_sp*vol_s + t_s*(g_stp - g_sp*vol_t + k_stp)
        + t_p*(g_sst + t_s*g_stt - (g0_ss - th_s**2*g0_tt)*tpinv))
    cabb = dalpha_tc + 2*ratio*dalpha_s - ratio**2*dbeta_s
    return cabb

def cabb_tcon_alt(s_gkg,t_cels,p_dbar,dtemp=_DTEMP,dsalt=_DSALT,
    tp_cels=None,chkvals=False,chktol=_CHKTOL,tp_cels0=None,
    chkbnd=False,useext=False,mathargs=None):
    """Calculate cabbeling wrt conservative temperature.
    
    Calculate the cabbeling coefficient of seawater with respect to
    conservative temperature. This version uses finite differences to
    approximate derivatives.
    
    :arg float s_gkg: Absolute salinity in g/kg.
    :arg float t_cels: In-situ temperature in degrees Celsius.
    :arg float p_dbar: In-situ gauge pressure in decibar.
    :arg float dtemp: Amount to change the conservative temperature by
        when calculating finite differences (default _DTEMP).
    :arg float dsalt: Amount to change the salinity by when calculating
        finite differences (default _DSALT).
    :arg tp_cels: Potential temperature in degrees Celsius. If unknown,
        pass None (default) and it will be calculated.
    :type tp_cels: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg tp_cels0: Initial guess for the potential temperature in
        degrees Celsius. If None (default) then `_approx_stppp` is used.
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg bool useext: If False (default) then the salt contribution is
        calculated from _GSCOEFFS; if True, from _GSCOEFFS_EXT.
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Cabbeling coefficient in 1/(deg C)^2.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> cabb_tcon_alt(35.,20.,1e3)
    8.96907383083e-6
    """
    tp_cels = pottemp(s_gkg,t_cels,p_dbar,0.,tp_cels=tp_cels,chkvals=chkvals,
        chktol=chktol,tp_cels0=tp_cels0,chkbnd=chkbnd,useext=useext,
        mathargs=mathargs)
    alpha0 = alpha_tcon(s_gkg,t_cels,p_dbar,tp_cels=tp_cels,useext=useext)
    beta0 = beta_tcon(s_gkg,t_cels,p_dbar,tp_cels=tp_cels,useext=useext)
    ratio = alpha0/beta0
    
    # Calculate finite difference wrt conservative temperature
    tc_cels = tconfromtpot(s_gkg,tp_cels,useext=useext)
    tc_l = tc_cels - dtemp
    tc_u = tc_cels + dtemp
    tp_l = tpotfromtcon(s_gkg,tc_l,useext=useext,mathargs=mathargs)
    tp_u = tpotfromtcon(s_gkg,tc_u,useext=useext,mathargs=mathargs)
    t_l = pottemp(s_gkg,tp_l,0.,p_dbar,useext=useext,mathargs=mathargs)
    t_u = pottemp(s_gkg,tp_u,0.,p_dbar,useext=useext,mathargs=mathargs)
    alpha_l = alpha_tcon(s_gkg,t_l,p_dbar,tp_cels=tp_l,useext=useext)
    alpha_u = alpha_tcon(s_gkg,t_u,p_dbar,tp_cels=tp_u,useext=useext)
    dalpha_tc = (alpha_u - alpha_l) / (tc_u - tc_l)
    
    # Calculate finite difference wrt salinity
    if s_gkg >= dsalt:
        s_l = s_gkg - dsalt
        s_u = s_gkg + dsalt
    else:
        s_l = 0.
        s_u = dsalt
    tp_l = tpotfromtcon(s_l,tc_cels,useext=useext,mathargs=mathargs)
    tp_u = tpotfromtcon(s_u,tc_cels,useext=useext,mathargs=mathargs)
    t_l = pottemp(s_l,tp_l,0.,p_dbar,useext=useext,mathargs=mathargs)
    t_u = pottemp(s_u,tp_u,0.,p_dbar,useext=useext,mathargs=mathargs)
    alpha_l = alpha_tcon(s_l,t_l,p_dbar,tp_cels=tp_l,useext=useext)
    alpha_u = alpha_tcon(s_u,t_u,p_dbar,tp_cels=tp_u,useext=useext)
    beta_l = beta_tcon(s_l,t_l,p_dbar,tp_cels=tp_l,useext=useext)
    beta_u = beta_tcon(s_u,t_u,p_dbar,tp_cels=tp_u,useext=useext)
    dalpha_s = (alpha_u - alpha_l) / (s_u - s_l)
    dbeta_s = (beta_u - beta_l) / (s_u - s_l)
    cabb = dalpha_tc + 2*ratio*dalpha_s - ratio**2*dbeta_s
    return cabb

def thrmb_tcon(s_gkg,t_cels,p_dbar,tp_cels=None,chkvals=False,
    chktol=_CHKTOL,tp_cels0=None,chkbnd=False,useext=False,
    mathargs=None):
    """Calculate thermobaricity wrt conservative temperature.
    
    Calculate the thermobaric coefficient of seawater with respect to
    conservative temperature. This version uses the exact derivatives of
    the polynomial expression for the free energy of seawater.
    
    :arg float s_gkg: Absolute salinity in g/kg.
    :arg float t_cels: In-situ temperature in degrees Celsius.
    :arg float p_dbar: In-situ gauge pressure in decibar.
    :arg tp_cels: Potential temperature in degrees Celsius. If unknown,
        pass None (default) and it will be calculated.
    :type tp_cels: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg tp_cels0: Initial guess for the potential temperature in
        degrees Celsius. If None (default) then `_approx_stppp` is used.
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg bool useext: If False (default) then the salt contribution is
        calculated from _GSCOEFFS; if True, from _GSCOEFFS_EXT.
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Thermobaric coefficient in 1/dbar/(deg C).
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> thrmb_tcon(35.,20.,1e3)
    1.7272850266583742e-08
    """
    tp_cels = pottemp(s_gkg,t_cels,p_dbar,0.,tp_cels=tp_cels,chkvals=chkvals,
        chktol=chktol,tp_cels0=tp_cels0,chkbnd=chkbnd,useext=useext,
        mathargs=mathargs)
    alpha0 = alpha_tcon(s_gkg,t_cels,p_dbar,tp_cels=tp_cels,useext=useext)
    beta0 = beta_tcon(s_gkg,t_cels,p_dbar,tp_cels=tp_cels,useext=useext)
    
    g_p = gsw_g(0,0,1,s_gkg,t_cels,p_dbar,useext=useext)
    g_st = gsw_g(1,1,0,s_gkg,t_cels,p_dbar,useext=useext)
    g_sp = gsw_g(1,0,1,s_gkg,t_cels,p_dbar,useext=useext)
    g_tt = gsw_g(0,2,0,s_gkg,t_cels,p_dbar,useext=useext)
    g_tp = gsw_g(0,1,1,s_gkg,t_cels,p_dbar,useext=useext)
    g_pp = gsw_g(0,0,2,s_gkg,t_cels,p_dbar,useext=useext)
    g_stt = gsw_g(1,2,0,s_gkg,t_cels,p_dbar,useext=useext)
    g_stp = gsw_g(1,1,1,s_gkg,t_cels,p_dbar,useext=useext)
    g_spp = gsw_g(1,0,2,s_gkg,t_cels,p_dbar,useext=useext)
    g_ttt = gsw_g(0,3,0,s_gkg,t_cels,p_dbar,useext=useext)
    g_ttp = gsw_g(0,2,1,s_gkg,t_cels,p_dbar,useext=useext)
    g_tpp = gsw_g(0,1,2,s_gkg,t_cels,p_dbar,useext=useext)
    g0_s = gsw_g(1,0,0,s_gkg,tp_cels,0.,useext=useext)
    
    dp = g_p**(-1)
    dtt = g_tt**(-1)
    dptt = dp*dtt
    thinv = (_TCELS + tp_cels)**(-1)
    t_s = dtt * (-g_st + g0_s*thinv)
    t_p = -g_tp*dtt
    vol_p = dp * (g_pp + t_p*g_tp)
    cp_p = dtt * (g_ttp + t_p*g_ttt)
    k_tpp = g_tpp + t_p*g_ttp - g_tp*(vol_p + cp_p)
    
    dalpha_p = -_CSEA*thinv * dptt * k_tpp
    dbeta_p = -dp * (g_spp + 2*t_p*g_stp + t_p**2*g_stt - g_sp*vol_p
        + t_s*k_tpp)
    thermb = dalpha_p - alpha0/beta0*dbeta_p
    thermb *= _DBAR2PA  # Since g_p is given in (J/kg)/Pa
    return thermb

def thrmb_tcon_alt(s_gkg,t_cels,p_dbar,dpres=_DPRES,tp_cels=None,
    chkvals=False,chktol=_CHKTOL,tp_cels0=None,chkbnd=False,
    useext=False,mathargs=None):
    """Calculate thermobaricity wrt conservative temperature.
    
    Calculate the thermobaric coefficient of seawater with respect to
    conservative temperature. This version uses finite differences to
    approximate derivatives.
    
    :arg float s_gkg: Absolute salinity in g/kg.
    :arg float t_cels: In-situ temperature in degrees Celsius.
    :arg float p_dbar: In-situ gauge pressure in decibar.
    :arg float dpres: Amount to change the pressure by when calculating
        finite differences (default _DPRES).
    :arg tp_cels: Potential temperature in degrees Celsius. If unknown,
        pass None (default) and it will be calculated.
    :type tp_cels: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg tp_cels0: Initial guess for the potential temperature in
        degrees Celsius. If None (default) then `_approx_stppp` is used.
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg bool useext: If False (default) then the salt contribution is
        calculated from _GSCOEFFS; if True, from _GSCOEFFS_EXT.
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Thermobaric coefficient in 1/dbar/(deg C).
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> thrmb_tcon_alt(35.,20.,1e3)
    1.72708365652e-8
    """
    tp_cels = pottemp(s_gkg,t_cels,p_dbar,0.,tp_cels=tp_cels,chkvals=chkvals,
        chktol=chktol,tp_cels0=tp_cels0,chkbnd=chkbnd,useext=useext,
        mathargs=mathargs)
    beta0 = beta_tcon(s_gkg,t_cels,p_dbar,tp_cels=tp_cels,useext=useext)
    
    # Calculate finite difference wrt pressure
    if p_dbar >= dpres:
        p_l = p_dbar - dpres
        p_u = p_dbar + dpres
    else:
        p_l = 0.
        p_u = dpres
    t_l = pottemp(s_gkg,tp_cels,0.,p_l,useext=useext,mathargs=mathargs)
    t_u = pottemp(s_gkg,tp_cels,0.,p_u,useext=useext,mathargs=mathargs)
    alpha_l = alpha_tcon(s_gkg,t_l,p_l,tp_cels=tp_cels,useext=useext)
    alpha_u = alpha_tcon(s_gkg,t_u,p_u,tp_cels=tp_cels,useext=useext)
    beta_l = beta_tcon(s_gkg,t_l,p_l,tp_cels=tp_cels,useext=useext)
    beta_u = beta_tcon(s_gkg,t_u,p_u,tp_cels=tp_cels,useext=useext)
    thermb = beta0 * (alpha_u/beta_u - alpha_l/beta_l) / (p_u - p_l)
    return thermb

