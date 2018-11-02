"""Physical constants and recommended bounds.

This module defines constants and parameters that are common across all
modules, including physical constants. It also contains functions that
produce errors for invalid inputs and warnings for out-of-range inputs.

The recommended bounds are the constants ``[SUB]_[V][LIM]``, where
``SUB`` is the substance it applies to; ``V`` is the variable (``T`` for
temperature, ``D`` for density); and ``LIM`` is either ``MIN`` or
``MAX``. ``CHKTOL`` is included here as a recommended tolerance for
iteration and error checking. All other constants are physical
constants. Those with the ``SO`` prefix involve the standard seawater
values; the ``CP`` prefix is for the vapour-liquid critical point of
water; and the ``TP`` prefix is for the vapour-liquid-ice triple point
of water.

In addition, there are several heat capacity values ``C[SUB]`` which are
not in the original TEOS routines. They are used in various
approximations that are used to initialize iterative routines.

:Functions:

* :func:`chkdrybnds`: Check inputs to the dry air functions.
* :func:`chkflubnds`: Check inputs to the fluid water functions.
* :func:`chkhumbnds`: Check inputs to the humid air functions.
* :func:`chkicebnds`: Check inputs to the ice functions.
* :func:`chksalbnds`: Check inputs to the seawater functions.

"""

__all__ = ['DRY_TMIN','DRY_TMAX','DRY_DMAX',
           'HUM_TMIN','HUM_TMAX','HUM_DMAX',
           'FLU_TMIN','FLU_TMAX','FLU_DMAX',
           'ICE_TMAX','ICE_PMAX',
           'SAL_SMAX','SAL_TMIN','SAL_TMAX','SAL_PMIN','SAL_PMAX',
           'GAS_CONSTANT_MOLAR_SI','GAS_CONSTANT_MOLAR_L2000',
           'MOLAR_MASS_H2O_SI','MOLAR_MASS_SEASALT_SI','MOLAR_MASS_AIR_SI',
           'MOLAR_MASS_AIR_L2000','GAS_CONSTANT_H2O_IAPWS95',
           'GAS_CONSTANT_H2O_SI','GAS_CONSTANT_AIR_SI',
           'GAS_CONSTANT_AIR_L2000','SEALEVEL_PRESSURE_SI',
           'CELSIUS_TEMPERATURE_SI',
           'SO_SALINITY_SI','SO_TEMPERATURE_SI','SO_PRESSURE_SI',
           'CP_DENSITY_SI','CP_TEMPERATURE_SI','CP_PRESSURE_SI',
           'CP_CHEMPOT_SI',
           'TP_TEMPERATURE_SI','TP_PRESSURE_EXP_SI','TP_PRESSURE_IAPWS95_SI',
           'TP_DENSITY_VAP_IAPWS95_SI','TP_DENSITY_LIQ_IAPWS95_SI',
           'TP_DENSITY_ICE_IAPWS95_SI','TP_ENTHALPY_VAP_SI',
           'TP_ENTHALPY_ICE_SI',
           'CICE','CVAP','CLIQ','CDRY','CSEA','SAL0','CHKTOL',
           'RUNIV','MWAT','MSAL','MDRY','RWAT_IAPWS95','RDRY','RDRY_L2000',
           'PATM','TCELS','SAL1','DCP','TCP','PCP','TTP','PTP','DVTP','DLTP',
           'DITP','LLVTP','LILTP',
           'chkdrybnds','chkflubnds','chkhumbnds','chkicebnds','chksalbnds']

import warnings

# Recommended limits
DRY_TMIN = 60.  # Dry air minimum temperature (K)
DRY_TMAX = 873.  # Dry air maximum temperature (K)
DRY_DMAX = 1035.  # Dry air maximum density (kg/m3)
HUM_TMIN = 193.  # Moist air minimum temperature (K)
HUM_TMAX = 473.  # Moist air maximum temperature (K)
HUM_DMAX = 100.  # Moist air maximum density (kg/m3)
FLU_TMIN = 50.  # Fluid water minimum temperature (K)
FLU_TMAX = 1273.  # Fluid water maximum temperature (K)
FLU_DMAX = 1240.  # Fluid water maximum density (kg/m3)
ICE_TMAX = 273.16  # Ice maximum temperature (K)
ICE_PMAX = 2e8  # Ice maximum pressure (Pa)
SAL_SMAX = 0.12  # Seawater maximum salinity (kg/kg)
SAL_TMIN = 262.  # Seawater minimum temperature (K)
SAL_TMAX = 353.  # Seawater maximum temperature (K)
SAL_PMIN = 100.  # Seawater minimum pressure (Pa)
SAL_PMAX = 1e8  # Seawater maximum pressure (Pa)

# Physical constants
GAS_CONSTANT_MOLAR_SI = 8.314472  # Molar gas constant R (J/mol/K)
GAS_CONSTANT_MOLAR_L2000 = 8.31451  # Gas constant from Lemmon 2000
MOLAR_MASS_H2O_SI = 0.018015268  # Molar mass of H2O (kg/mol)
MOLAR_MASS_SEASALT_SI = 0.0314038218  # Molar mass of sea salt (kg/mol)
MOLAR_MASS_AIR_SI = 0.02896546  # Molar mass of dry air (kg/mol)
MOLAR_MASS_AIR_L2000 = MOLAR_MASS_AIR_SI
GAS_CONSTANT_H2O_IAPWS95 = 461.51805  # Specific gas constant of H2O (J/kg/K)
GAS_CONSTANT_H2O_SI = GAS_CONSTANT_MOLAR_SI / MOLAR_MASS_H2O_SI
GAS_CONSTANT_AIR_SI = GAS_CONSTANT_MOLAR_SI / MOLAR_MASS_AIR_SI
GAS_CONSTANT_AIR_L2000 = GAS_CONSTANT_MOLAR_L2000 / MOLAR_MASS_AIR_L2000
SEALEVEL_PRESSURE_SI = 101325.  # 1 atm in Pa
CELSIUS_TEMPERATURE_SI = 273.15  # 0 Celsius in K
SO_SALINITY_SI = 0.03516504  # Salinity of KCl-normalised water (kg/kg)
SO_TEMPERATURE_SI = CELSIUS_TEMPERATURE_SI
SO_PRESSURE_SI = SEALEVEL_PRESSURE_SI
CP_DENSITY_SI = 322.  # Water critical point density (kg/m3)
CP_TEMPERATURE_SI = 647.096  # Water critical point temperature (K)
CP_PRESSURE_SI = 22064000.  # Water critical point pressure (Pa)
CP_CHEMPOT_SI = -767471.156792841  # Water critical point potential (J/kg)
TP_TEMPERATURE_SI = 273.16  # Water triple point temperature (K)
TP_PRESSURE_EXP_SI = 611.657  # Water triple point pressure (Pa)
TP_PRESSURE_IAPWS95_SI = 611.654771007894
TP_DENSITY_VAP_IAPWS95_SI = 4.85457572477859e-03  # TP vapour density (kg/m3)
TP_DENSITY_LIQ_IAPWS95_SI = 999.792520031621  # TP liquid density (kg/m3)
TP_DENSITY_ICE_IAPWS95_SI = 916.709492199488  # TP ice density (kg/m3)
TP_ENTHALPY_VAP_SI = 2500915.1914657  # TP enthalpy of vaporization (J/kg)
TP_ENTHALPY_ICE_SI = -333444.253967839  # TP enthalpy of fusion (J/kg)

# Reference values NOT included in the original Constants_0.f
CICE = 2096.7843  # Water triple point ice heat capacity (J/kg/K)
CVAP = 1884.3520  # Water triple point vapour heat capacity (J/kg/K)
CLIQ = 4219.9115  # Water triple point liquid heat capacity (J/kg/K)
CDRY = 1005.6844  # Heat capacity of dry air at 0 Celsius and 1 atm (J/kg/K)
CSEA = 3991.86795711963  # Standard heat capacity of seawater (J/kg/K)
SAL0 = 0.035  # Typical reference salinity (kg/kg)
CHKTOL = 1e-8  # Chosen tolerance for error checking

# Shorter aliases for common variables
RUNIV = GAS_CONSTANT_MOLAR_SI
MWAT = MOLAR_MASS_H2O_SI
MSAL = MOLAR_MASS_SEASALT_SI
MDRY = MOLAR_MASS_AIR_SI
RWAT = GAS_CONSTANT_H2O_SI
RWAT_IAPWS95 = GAS_CONSTANT_H2O_IAPWS95
RDRY = GAS_CONSTANT_AIR_SI
RDRY_L2000 = GAS_CONSTANT_AIR_L2000
PATM = SEALEVEL_PRESSURE_SI
TCELS = CELSIUS_TEMPERATURE_SI
SAL1 = SO_SALINITY_SI
DCP = CP_DENSITY_SI
TCP = CP_TEMPERATURE_SI
PCP = CP_PRESSURE_SI
TTP = TP_TEMPERATURE_SI
PTPE = TP_PRESSURE_EXP_SI
PTPI = TP_PRESSURE_IAPWS95_SI
DVTP = TP_DENSITY_VAP_IAPWS95_SI
DLTP = TP_DENSITY_LIQ_IAPWS95_SI
DITP = TP_DENSITY_ICE_IAPWS95_SI
LLVTP = TP_ENTHALPY_VAP_SI
LILTP = -TP_ENTHALPY_ICE_SI


### Warning functions
def chkdrybnds(temp,ddry,chkbnd=True,stacklevel=2):
    """Check inputs to the dry air functions.
    
    Check whether the given temperature and density are within the
    recommended bounds for dry air. Raises errors for non-positive
    values; raises warnings for out-of-bound values. Warnings can be
    suppressed by setting chkbnd to False.
    
    :arg float temp: Temperature in K.
    :arg float ddry: Density in kg/m3.
    :arg bool chkbnd: If True (default) then warnings are raised when
        the given values are valid but outside the recommended bounds.
    :arg int stacklevel: Controls how many levels deep to raise the
        warning from (default 2).
    :returns: None
    :raises ValueError: If temp or ddry are nonpositive (<=0).
    :raises RuntimeWarning: If temp < DRY_TMIN, temp > DRY_TMAX, or
        ddry > DRY_DMAX.
    """
    
    # Nonpositive temperature and density are always bad
    if temp <= 0:
        errmsg = 'The given temperature {0} is not positive'.format(temp)
        raise ValueError(errmsg)
    if ddry <= 0:
        errmsg = 'The given density {0} is not positive'.format(ddry)
        raise ValueError(errmsg)
    
    # Stop here if the warnings are suppressed
    if not chkbnd:
        return None
    
    # Check the bounds
    if temp < DRY_TMIN or temp > DRY_TMAX:
        warnmsg = ('Temperature {0} exceeds the recommended bounds '
                   '{1}-{2} K').format(temp,DRY_TMIN,DRY_TMAX)
        warnings.warn(warnmsg,RuntimeWarning,stacklevel=stacklevel)
    if ddry > DRY_DMAX:
        warnmsg = ('Density {0} exceeds the recommended upper bound '
            '{1} kg/m3').format(ddry,DRY_DMAX)
        warnings.warn(warnmsg,RuntimeWarning,stacklevel=stacklevel)
    return None

def chkflubnds(temp,dflu,chkbnd=True,stacklevel=2):
    """Check inputs to the fluid water functions.
    
    Check whether the given temperature and density are within the
    recommended bounds for fluid water. Raises errors for non-positive
    values; raises warnings for out-of-bound values. Warnings can be
    suppressed by setting chkbnd to False.
    
    :arg float temp: Temperature in K.
    :arg float dflu: Density in kg/m3.
    :arg bool chkbnd: If True (default) then warnings are raised when
        the given values are valid but outside the recommended bounds.
    :arg int stacklevel: Controls how many levels deep to raise the
        warning from (default 2).
    :returns: None
    :raises ValueError: If temp or dflu are nonpositive (<=0).
    :raises RuntimeWarning: If temp < FLU_TMIN, temp > FLU_TMAX, or
        dflu > FLU_DMAX.
    """
    
    # Nonpositive temperature and dfluity are always bad
    if temp <= 0:
        errmsg = 'The given temperature {0} is not positive'.format(temp)
        raise ValueError(errmsg)
    if dflu <= 0:
        errmsg = 'The given density {0} is not positive'.format(dflu)
        raise ValueError(errmsg)
    
    # Stop here if the warnings are suppressed
    if not chkbnd:
        return None
    
    # Check the bounds
    if temp < FLU_TMIN or temp > FLU_TMAX:
        warnmsg = ('Temperature {0} exceeds the recommended bounds '
            '{1}-{2} K').format(temp,FLU_TMIN,FLU_TMAX)
        warnings.warn(warnmsg,RuntimeWarning,stacklevel=stacklevel)
    if dflu > FLU_DMAX:
        warnmsg = ('Density {0} exceeds the recommended upper bound '
            '{1} kg/m3').format(dflu,FLU_DMAX)
        warnings.warn(warnmsg,RuntimeWarning,stacklevel=stacklevel)
    return None

def chkhumbnds(airf,temp,dhum,chkbnd=True,stacklevel=2):
    """Check inputs to the humid air functions.
    
    Check whether the given dry air mass fraction, temperature, and
    density are within the recommended bounds for humid air. Raises
    errors for non- positive values; raises warnings for out-of-bound
    values. Warnings can be suppressed by setting chkbnd to False.
    
    :arg float airf: Dry air mass fraction in kg/kg.
    :arg float temp: Temperature in K.
    :arg float dhum: Density in kg/m3.
    :arg bool chkbnd: If True (default) then warnings are raised when
        the given values are valid but outside the recommended bounds.
    :arg int stacklevel: Controls how many levels deep to raise the
        warning from (default 2).
    :returns: None
    :raises ValueError: If temp or dhum are nonpositive (<=0), or if
        airf is not between 0 and 1.
    :raises RuntimeWarning: If temp < HUM_TMIN, temp > HUM_TMAX, or
        ddry > HUM_DMAX.
    """
    
    # Nonpositive temperature and density are always bad
    if temp <= 0:
        errmsg = 'The given temperature {0} is not positive'.format(temp)
        raise ValueError(errmsg)
    if dhum <= 0:
        errmsg = 'The given density {0} is not positive'.format(dhum)
        raise ValueError(errmsg)
    if airf < 0 or airf > 1:
        errmsg = ('The given dry air mass fraction {0} is not between '
            '0 and 1').format(airf)
        raise ValueError(errmsg)
    
    # Stop here if the warnings are suppressed
    if not chkbnd:
        return None
    
    # If the air is dry, use those bounds instead
    if airf == 1.:
        chkdrybnds(temp,dhum,chkbnd=chkbnd,stacklevel=stacklevel)
        return None
    
    # If the air is pure water vapour, use those bounds instead
    if airf == 0.:
        chkflubnds(temp,dhum,chkbnd=chkbnd,stacklevel=stacklevel)
        return None
    
    # Check the humid air bounds
    if temp < HUM_TMIN or temp > HUM_TMAX:
        warnmsg = ('Temperature {0} exceeds the recommended bounds '
            '{1}-{2} K').format(temp,HUM_TMIN,HUM_TMAX)
        warnings.warn(warnmsg,RuntimeWarning,stacklevel=stacklevel)
    if dhum > HUM_DMAX:
        warnmsg = ('Density {0} exceeds the recommended upper bound '
            '{1} kg/m3').format(dhum,HUM_DMAX)
        warnings.warn(warnmsg,RuntimeWarning,stacklevel=stacklevel)
    return None

def chkicebnds(temp,pres,chkbnd=True,stacklevel=2):
    """Check inputs to the ice functions.
    
    Check whether the given temperature and pressure are within the
    recommended bounds for ice. Raises errors for non-positive values;
    raises warnings for out-of-bound values. Warnings can be suppressed
    by setting chkbnd to False.
    
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg bool chkbnd: If True (default) then warnings are raised when
        the given values are valid but outside the recommended bounds.
    :arg int stacklevel: Controls how many levels deep to raise the
        warning from (default 2).
    :returns: None
    :raises ValueError: If temp or pres are nonpositive (<=0).
    :raises RuntimeWarning: If temp > ICE_TMAX or pres > ICE_PMAX.
    """
    
    # Nonpositive temperature and pressure are always bad
    if temp <= 0:
        errmsg = 'The given temperature {0} is not positive'.format(temp)
        raise ValueError(errmsg)
    if pres <= 0:
        errmsg = 'The given pressure {0} is not positive'.format(pres)
        raise ValueError(errmsg)
    
    # Stop here if the warnings are suppressed
    if not chkbnd:
        return None
    
    # Check the bounds
    if temp > ICE_TMAX:
        warnmsg = ('Temperature {0} exceeds the recommended upper '
            'bound {1} K').format(temp,ICE_TMAX)
        warnings.warn(warnmsg,RuntimeWarning,stacklevel=stacklevel)
    if pres > ICE_PMAX:
        warnmsg = ('Pressure {0} exceeds the recommended upper bound '
            '{1} Pa').format(pres,ICE_PMAX)
        warnings.warn(warnmsg,RuntimeWarning,stacklevel=stacklevel)
    return None

def chksalbnds(salt,temp,pres,chkbnd=True,stacklevel=2):
    """Check inputs to the salinity functions.
    
    Check whether the given salinity, temperature, and density are
    within the recommended bounds for salt in seawater. Raises errors
    for non-positive values; raises warnings for out-of-bound values.
    Warnings can be suppressed by setting chkbnd to False.
    
    :arg float salt: Salinity in kg/kg.
    :arg float temp: Temperature in K.
    :arg float ddry: Density in kg/m3.
    :arg bool chkbnd: If True (default) then warnings are raised when
        the given values are valid but outside the recommended bounds.
    :arg int stacklevel: Controls how many levels deep to raise the
        warning from (default 2).
    :returns: None
    :raises ValueError: If temp or pres are nonpositive (<=0), or if
        salt is not between 0 and 1.
    :raises RuntimeWarning: If salt > SAL_SMAX, temp < SAL_TMIN,
        temp > SAL_TMAX, pres < SAL_PMIN, or pres > SAL_PMAX.
    """
    
    # Nonpositive temperature and pressure are always bad
    if temp <= 0:
        errmsg = 'The given temperature {0} is not positive'.format(temp)
        raise ValueError(errmsg)
    if pres <= 0:
        errmsg = 'The given pressure {0} is not positive'.format(pres)
        raise ValueError(errmsg)
    if salt < 0 or salt >= 1:
        errmsg = 'The given salinity {0} is not between 0 and 1'.format(salt)
    
    # Stop here if the warnings are suppressed
    if not chkbnd:
        return None
    
    # Check the bounds
    if salt > SAL_SMAX:
        warnmsg = ('Salinity {0} exceeds the recommended upper bound '
            '{1} kg/kg').format(salt,SAL_SMAX)
        warnings.warn(warnmsg,RuntimeWarning,stacklevel=stacklevel)
    if temp < SAL_TMIN or temp > SAL_TMAX:
        warnmsg = ('Temperature {0} exceeds the recommended bounds '
            '{1}-{2} K').format(temp,SAL_TMIN,SAL_TMAX)
        warnings.warn(warnmsg,RuntimeWarning,stacklevel=stacklevel)
    if pres < SAL_PMIN or pres > SAL_PMAX:
        warnmsg = ('Pressure {0} exceeds the recommended bounds '
            '{1}-{2} Pa').format(temp,SAL_PMIN,SAL_PMAX)
        warnings.warn(warnmsg,RuntimeWarning,stacklevel=stacklevel)
    return None


