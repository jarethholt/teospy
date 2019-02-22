"""Humid air properties in non-SI units.

This module provides the lapse rate of humid air in the common units of
degrees Celsius per 100 m.

:Examples:

>>> lapserate_c100m(50.,20.,1e3)
0.971588085046

:Functions:

* :func:`lapserate_c100m`: Calculate the dry adiabatic lapse rate as the
  change in temperature in degrees Celsius per 100 meters altitude.
"""

__all__ = ['lapserate_c100m']

from teospy import constants0
from teospy import air3a
from teospy import air3b
from teospy import liqair4a

_CHKTOL = constants0.CHKTOL
_TCELS = constants0.TCELS
_HPA2PA = 100.
_PCT2FRAC = .01
_DALT = 100.  # Common altitude difference in m
_GRAV = 9.81  # Standard gravity in m/s2

def lapserate_c100m(rh_pct,t_cels,p_hpa,dhum=None,asat=None,dhsat=None,
    dlsat=None,chkvals=False,chktol=_CHKTOL,dhum0=None,asat0=None,
    dhsat0=None,dlsat0=None,chkbnd=False,mathargs=None):
    """Calculate adiabatic lapse rate wrt altitude.
    
    Calculate the dry adiabatic lapse rate of humid air, in terms of an
    altitude difference of 100 m rather than a pressure difference of 1
    Pa. The relative humidity here uses the WMO definition,
    
        rh_wmo = [(1-airf)/airf] / [(1-asat)/asat]
    
    where airf is the total dry air mass fraction and asat is the mass
    fraction at saturation.
    
    :arg float rh_pct: Relative humidity as a percentage (0-100). Note
        that the moist lapse rate is not provided here.
    :arg float t_cels: Temperature in degrees Celsius.
    :arg float p_hpa: Pressure in hPa.
    :arg dhum: Humid air density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dhum: float or None
    :arg asat: Humid air dry fraction at saturation in kg/kg. If
        unknown, pass None (default) and it will be calculated.
    :type asat: float or None
    :arg dhsat: Humid air density at saturation in kg/m3. If unknown,
        pass None (default) and it will be calculated.
    :type dhsat: float or None
    :arg dlsat: Liquid water density at saturation in kg/m3. If unknown,
        pass None (default) and it will be calculated.
    :type dlsat: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg dhum0: Initial guess for the humid air density in kg/m3. If
        None (default) then `liqair4._approx_tp` is used.
    :type dhum0: float or None
    :arg asat0: Initial guess for the humid air dry fraction at
        saturation in kg/m3. If None (default) then `liqair4._approx_tp`
        is used.
    :type asat0: float or None
    :arg dhsat0: Initial guess for the humid air density at saturation
        in kg/m3. If None (default) then `_approx_tp` is used.
    :type dhsat0: float or None
    :arg dlsat0: Initial guess for the liquid water density at
        saturation in kg/m3. If None (default) then `liqair4._approx_tp`
        is used.
    :type dlsat0: float or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Lapse rate in degrees Celsius per 100 meters.
    
    :Examples:
    
    >>> lapserate_c100m(50.,20.,1e3)
    0.971588085046
    """
    rh_wmo = rh_pct * _PCT2FRAC
    temp = t_cels + _TCELS
    pres = p_hpa * _HPA2PA
    
    asat, __, __, dhsat, dlsat = liqair4a.eq_atpe(temp=temp,pres=pres,airf=asat,
        dhum=dhsat,dliq=dlsat,chkvals=chkvals,chktol=chktol,airf0=asat0,
        dhum0=dhsat0,dliq0=dlsat0,chkbnd=chkbnd,mathargs=mathargs)
    airf = liqair4a.airffromrh_wmo(rh_wmo,temp,pres,asat=asat,dhsat=dhsat,
        dlsat=dlsat)
    dhum = air3a.eq_atp(airf,temp,pres,dhum=dhum,chkvals=chkvals,chktol=chktol,
        dhum0=dhum0,chkbnd=chkbnd,mathargs=mathargs)
    gamma_kpa = air3b.lapserate(airf,temp,pres,dhum=dhum)
    gamma_c100m = gamma_kpa * dhum*_GRAV*_DALT
    return gamma_c100m

