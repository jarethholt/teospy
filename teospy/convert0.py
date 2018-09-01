"""Conversion between common units and SI units.

This module provides functions for converting between mass fractions and
mole fractions for both air and seawater, and from practical salinity
units to the IAPWS-recommended absolute salinity.

Practical salinity is calculated in part from a lookup table. The
datafile name of this table is GSW_FNAME. When practical salinity cannot
be calculated from the table, the value GSW_ERRVAL is returned.

Functions:

* air_molarmass: Average molar mass of humid air.
* air_molfractionvap: Mole fraction of water vapour in humid air.
* air_molfractiondry: Mole fraction of dry air in humid air.
* air_massfractionvap: Mass fraction of water vapour in humid air.
* air_massfractiondry: Mass fraction of dry air in humid air.
* sal_molality: Molality of salt in seawater.
* gsw_safromsp: Convert practical to absolute salinity with inputs in
  GSW units.
* gsw_spfromsa: Convert absolute to practical salinity with inputs in
  GSW units.
* sal_asalfrompsal: Convert practical to absolute salinity with inputs
  in SI units.
* sal_psalfromasal: Convert absolute to practical salinity with inputs
  in SI units.

"""

__all__ = ['GSW_FNAME','GSW_ERRVAL',
    'air_molarmass','air_molfractionvap','air_molfractiondry',
    'air_massfractionvap','air_massfractiondry','sal_molality',
    'gsw_safromsp','gsw_spfromsa','sal_asalfrompsal','sal_psalfromasal']

import warnings
import constants0

# Physical constants
_MWAT = constants0.MWAT
_MDRY = constants0.MDRY
_MSAL = constants0.MSAL
_PATM = constants0.PATM
_SAL0 = constants0.SAL0
_SAL1 = constants0.SAL1
_KG2G = 1e3
_DBAR2PA = 1e4

# Constants related to practical salinity functions
GSW_ERRVAL = 9e15
GSW_FNAME = 'GSW_Data_v3_0.dat'
_ERRVAL = 1e10
_GSW_NS = (91,45,45)
_DTAS = ((0,0), (0,1), (1,1), (1,0))
_BALTIC = (((12.6,7.,26.),(50.,59.,69.)), ((45.,26.),(50.,69.)))
_PANAMA = ((260.00,272.59,276.50,278.65,280.73,292.00),
    (19.55,13.97,9.6,8.1,9.33,3.4))
_SAARLIM = 100.
_DSAL = 0.087
_LATLIMS = (-86.,90.)
_GSW_GRID = None


### Mass fraction functions
def air_molarmass(airf):
    """Calculate humid air molar mass.
    
    Calculate the average molar mass of humid air for the given dry air
    mass fraction.
    
    :arg float airf: Dry air mass fraction in kg/kg.
    :returns: Molar mass of humid air in kg/mol.
    
    :Example:
    
    >>> air_molarmass(0.5)
    0.022214237490882645
    """
    massmol = airf/_MDRY + (1-airf)/_MWAT
    molmass = 1./massmol
    return molmass

def air_molfractionvap(airf):
    """Calculate water vapour mole fraction.
    
    Calculate the mole fraction of water vapour in humid air for the
    given dry air mass fraction.
    
    :arg float airf: Dry air mass fraction in kg/kg.
    :returns: Mole fraction of water vapour in mol/mol.
    
    :Example:
    
    >>> air_molfractionvap(0.5)
    0.6165391902824494
    """
    ntot = airf/_MDRY + (1-airf)/_MWAT
    nvap = (1-airf)/_MWAT
    vapx = nvap / ntot
    return vapx

def air_molfractiondry(airf):
    """Calculate dry air mole fraction.
    
    Calculate the mole fraction of dry air in humid air for the given
    dry air mass fraction.
    
    :arg float airf: Dry air mass fraction in kg/kg.
    :returns: Mole fraction of dry air in mol/mol.
    
    :Example:
    
    >>> air_molfractiondry(0.5)
    0.3834608097175506
    """
    ntot = airf/_MDRY + (1-airf)/_MWAT
    ndry = airf/_MDRY
    airx = ndry / ntot
    return airx

def air_massfractionvap(airx):
    """Calculate water vapour mass fraction.
    
    Calculate the mass fraction of water vapour in humid air for the
    given dry air mole fraction.
    
    :arg float airx: Dry air mole fraction in mol/mol.
    :returns: Mass fraction of water vapour in kg/kg.
    
    :Example:
    
    >>> air_massfractionvap(0.5)
    0.3834608097175506
    """
    mtot = airx*_MDRY + (1-airx)*_MWAT
    mvap = (1-airx)*_MWAT
    vapf = mvap / mtot
    return vapf

def air_massfractiondry(airx):
    """Calculate dry air mass fraction.
    
    Calculate the mass fraction of dry air in humid air for the given
    dry air mole fraction.
    
    :arg float airx: Dry air mole fraction in mol/mol.
    :returns: Mass fraction of dry air in kg/kg.
    
    :Example:
    
    >>> air_massfractionair(0.5)
    0.6165391902824494
    """
    mtot = airx*_MDRY + (1-airx)*_MWAT
    mdry = airx*_MDRY
    airf = mdry / mtot
    return airf

def sal_molality(salt):
    """Calculate molality of salt in seawater.
    
    Calculate the molality of salt in seawater for the given absolute
    salinity. The molality is the ratio of the amount of solute (salt,
    in moles) to the mass of solvent (pure water, in kg).
    
    :arg float salt: Salinity in kg/kg.
    :returns: Molality of salt in mol/kg.
    
    :Example:
    
    >>> sal_molality(0.035)
    1.15493681893
    """
    molal = salt / ((1-salt) * _MSAL)
    return molal


### Auxiliary functions for the practical salinity unit calculations
def _getind(x,x0):
    """Find the index of a point on a grid.
    
    Find the index k of x0 in a real array x, such that
        x[k] <= x0 < x[k+1].
    To enable usage by the interpolation routine, the edge cases are
        x0 <= x[0]: k = 0
        x0 >= x[-1]: k = len(x) - 2
    so that expressions such as (x[k+1]-x[k]) are always valid.
    
    :arg x: Grid points, monotonically increasing.
    :type x: list[float]
    :arg float x0: Value to find the index of.
    :returns: Index of x0 in the array.
    """
    if x0 <= x[0]:
        k = 0
    elif x0 >= x[-1]:
        k = len(x) - 2
    else:
        k = [(xk > x0) for xk in x].index(True) - 1
    return k

def _xinterp1(x,y,x0):
    """Interpolate linearly for a value.
    
    Linearly interpolate from gridded of (x,y) values to estimate y(x0).
    
    :arg x: Grid points, monotonically increasing.
    :type x: list[float]
    :arg y: Function values at the grid points.
    :type y: list[float]
    :arg float x0: Point to estimate the function at.
    :returns: Linearly interpolated estimate y(x0).
    """
    k = _getind(x,x0)
    r = (x0 - x[k]) / (x[k+1]-x[k])
    y0 = y[k] + r*(y[k+1]-y[k])
    return y0

def _read_saar():
    """Read the SAAR datafile.
    
    Read the data file GSW_FNAME of gridded Salinity Absolute Anomaly
    Ratio (SAAR) reference values. The data are saved to the tuple of
    arrays _GSW_GRID.
    """
    # Check if the array has already been read
    global _GSW_GRID
    if _GSW_GRID is not None:
        return None
    
    # Initialize arrays to hold data
    GSW_NX, GSW_NY, GSW_NZ = _GSW_NS[:]
    GSW_LONS = [0.0 for iX in range(GSW_NX)]
    GSW_LATS = [0.0 for iY in range(GSW_NY)]
    GSW_PS = [0.0 for iZ in range(GSW_NZ)]
    GSW_NDEPTHS = [[0 for iX in range(GSW_NX)] for iY in range(GSW_NY)]
    SAAR_REF = [[[0.0 for iX in range(GSW_NX)] for iY in range(GSW_NY)]
        for iZ in range(GSW_NZ)]
    
    ## Read data file
    with open(GSW_FNAME,'r') as gswdatafile:
        # Coordinate arrays
        for iX in range(GSW_NX):
            GSW_LONS[iX] = float(gswdatafile.readline())
        for iY in range(GSW_NY):
            GSW_LATS[iY] = float(gswdatafile.readline())
        for iZ in range(GSW_NZ):
            GSW_PS[iZ] = float(gswdatafile.readline())
        # Reference depths
        for iX in range(GSW_NX):
            for iY in range(GSW_NY):
                temp = float(gswdatafile.readline())
                if temp < 0:
                    GSW_NDEPTHS[iY][iX] = -1
                else:
                    GSW_NDEPTHS[iY][iX] = int(temp)
        # SAAR values
        for iX in range(GSW_NX):
            for iY in range(GSW_NY):
                for iZ in range(GSW_NZ):
                    SAAR_REF[iZ][iY][iX] = float(gswdatafile.readline())
    
    # Compile variables into a single tuple of arrays and return
    _GSW_GRID = (GSW_LONS, GSW_LATS, GSW_PS, GSW_NDEPTHS, SAAR_REF)
    return None

def _check_baltic(lon0,lat0):
    """Check that a point is within the Baltic sea.
    
    Check that the given point is within the Baltic sea. See _BALTIC for
    the longitude and latitude limits.
    
    :arg float lon0: Reference longitude in degrees East.
    :arg float lat0: Reference latitude in degrees North.
    :returns: True if the point is within the Baltic.
    """
    # Is it within the widest limits?
    XBL, YBL = _BALTIC[0]
    XBR, YBR = _BALTIC[1]
    if (lon0<XBL[-1] or lon0>XBR[0] or lat0<YBL[0] or lat0>YBL[-1]):
        return False
    
    # Is it within the narrow limits?
    xl = _xinterp1(YBL,XBL,lat0)
    xr = _xinterp1(YBR,XBR,lat0)
    if (lon0<xl or lon0>xr):
        return False
    return True

def _gsw_safromspbaltic(spsu,lon0,lat0):
    """Convert practical to absolute salinity for the Baltic sea.
    
    Calculate the absolute salinity from the practical salinity using an
    empirical function derived for the Baltic Sea.
    
    :arg float spsu: Salinity in practical salinity units (psu).
    :arg float lon0: Reference longitude in degrees East.
    :arg float lat0: Reference latitude in degrees North.
    :returns: Absolute salinity in g/kg.
    :raises RuntimeWarning: If the longitude and latitude are outside of
        the Baltic Sea limits, _BALTIC. The error value GSW_ERRVAL is
        returned.
    """
    # Check that the point is actually in the Baltic
    isbaltic = _check_baltic(lon0,lat0)
    if not isbaltic:
        warnmsg = ('Location ({0} E, {1} N) appears to be outside the Baltic '
            'sea').format(lon0,lat0)
        warnmsg = BALTICMSG.format(lon0,lat0)
        warnings.warn(warnmsg,RuntimeWarning)
        return GSW_ERRVAL
    
    # Calculate the salinity
    sgkg = (_SAL1*_KG2G-_DSAL)/(_SAL0*_KG2G)*spsu + _DSAL
    return sgkg

def _gsw_spfromsabaltic(sgkg,lon0,lat0):
    """Convert absolute to practical salinity for the Baltic sea.
    
    Calculate the practical salinity from the absolute salinity using an
    empirical function derived for the Baltic Sea.
    
    :arg float sgkg: Absolute salinity in g/kg.
    :arg float lon0: Reference longitude in degrees East.
    :arg float lat0: Reference latitude in degrees North.
    :returns: Practical salinity in practical salinity units (psu).
    :raises RuntimeWarning: If the longitude and latitude are outside of
        the Baltic Sea limits, _BALTIC. The error value GSW_ERRVAL is
        returned.
    """
    # Check that the point is actually in the Baltic
    isbaltic = _check_baltic(lon0,lat0)
    if not isbaltic:
        warnmsg = ('Location ({0} E, {1} N) appears to be outside the Baltic '
            'sea').format(lon0,lat0)
        warnmsg = BALTICMSG.format(lon0,lat0)
        warnings.warn(warnmsg,RuntimeWarning)
        return GSW_ERRVAL
    
    # Calculate the salinity
    spsu = (_SAL0*_KG2G)/(_SAL1*_KG2G-_DSAL) * (sgkg - _DSAL)
    return spsu

def _gsw_addbarrier(saars_in,lon0,lat0,lonsw,latsw,dlon,dlat):
    """Average data around the Panama barrier.
    
    Averages data on the appropriate side of the Central American
    barrier (Panama) separating the Pacific and Atlantic oceans.
    
    :arg saars_in: SAAR values near the point being calculated. These
        are the gridpoint SAAR values southwest, northwest, southeast,
        and northeast of the central point.
    :type saars_in: list[float]
    :arg float lon0: Reference longitude in degrees East.
    :arg float lat0: Reference latitude in degrees North.
    :arg float lonsw: Grid point longitude directly West of the
        reference.
    :arg float latsw: Grid point latitude directly South of the
        reference.
    :arg float dlon: Grid spacing for longitude in decimal degrees.
    :arg float dlat: Grid spacing for latitude in decimal degrees.
    :returns: SAAR values averaged correctly on each side of the Panama
        barrier.
    :rtype: list[float]
    """
    aboveline = [0 for __ in range(4)]
    PLONS, PLATS = _PANAMA[:]
    
    # Is the reference point above the barrier?
    k0 = _getind(PLONS,lon0)
    r = (lon0 - PLONS[k0]) / (PLONS[k0+1] - PLONS[k0])
    lats_line = PLATS[k0] + r*(PLATS[k0+1] - PLATS[k0])
    aboveline0 = (lats_line < lat0)
    
    # Are the western grid points above the barrier?
    kw = _getind(PLONS,lonsw)
    r = (lonsw - PLONS[kw]) / (PLONS[kw+1] - PLONS[kw])
    lats_line = PLATS[kw] + r*(PLATS[kw+1] - PLATS[kw])
    aboveline[0] = (lats_line < latsw)
    aboveline[3] = (lats_line < latsw + dlat)
    
    # Are the eastern grid points above the barrier?
    ke = _getind(PLONS,lonsw+dlon)
    r = (lonsw+dlon - PLONS[ke]) / (PLONS[ke+1] - PLONS[ke])
    lats_line = PLATS[ke] + r*(PLATS[ke+1] - PLATS[ke])
    aboveline[1] = (lats_line < latsw)
    aboveline[2] = (lats_line < latsw + dlat)
    
    # Add and count the valid SAAR values
    nmean = 0
    avg = 0.
    for k in range(4):
        if (abs(saars_in[k]) <= _SAARLIM and aboveline0 == aboveline[k]):
            nmean += 1
            avg += saars_in[k]
    if nmean == 0:
        avg = 0.
    else:
        avg /= nmean
    
    # Mask the invalid values with the average
    for k in range(4):
        if (abs(saars_in[k]) >= _ERRVAL or aboveline0 != aboveline[k]):
            saars_out[k] = avg
        else:
            saars_out[k] = saars[k]
    return saars_out

def _gsw_addmean(saars_in,lon0,lat0):
    """Fill in invalid numbers on a grid.
    
    Replaces invalid numbers in a data array with the mean of the
    adjacent valid neighbors.
    
    :arg saars_in: SAAR values near the point being calculated. These
        are the gridpoint SAAR values southwest, northwest, southeast,
        and northeast of the central point.
    :type saars_in: list[float]
    :arg float lon0: Reference longitude in degrees East.
    :arg float lat0: Reference latitude in degrees North.
    :returns: SAAR values with invalid entries replaced by the average
        over the neighboring points.
    :rtype: list[float]
    """
    # Add and count the valid SAAR values
    nmean = 0
    avg = 0.
    for k in range(4):
        if (abs(saars_in[k]) <= _SAARLIM):
            nmean += 1
            avg += saars_in[k]
    if nmean == 0:
        avg = 0.
    else:
        avg /= nmean
    
    # Replace invalid values with the average
    saars_out = [0. for __ in range(4)]
    for k in range(4):
        if abs(saars_in[k]) >= _SAARLIM:
            saars_out[k] = avg
        else:
            saars_out[k] = saars_in[k]
    return saars_out

def _gsw_saar(pdbar,lon0,lat0):
    """Calculate SAAR for a location.
    
    Calculate the Salinity Absolute Anomaly Ratio (SAAR) for a location
    by interpolating the reference table in GSW_FNAME. The units here
    are those used by the GSW toolbox, i.e. with pressure as gauge
    pressure in decibars.
    
    :arg float pdbar: Seawater (gauge) pressure in dbar.
    :arg float lon0: Reference longitude in degrees East.
    :arg float lat0: Reference latitude in degrees North.
    :returns: Interpolated SAAR value for this location.
    :raises RuntimeWarning: If the reference point is not over the
        ocean, as determined by the GSW_FNAME gridded data file.
        GSW_ERRVAL is returned.
    """
    # Check input values
    if lat0 < _LATLIMS[0] or lat0 > _LATLIMS[1]:
        warnmsg = ('Latitude {0} is not between {1} S and {2} '
            'N').format(lat0,-_LATLIMS[0],_LATLIMS[1])
        warnings.warn(warnmsg,RuntimeWarning)
        return GSW_ERRVAL
    lon1 = lon0 % 360  # Wrap longitude to between 0 and 360
    
    # Get the GSW data file SAAR values
    GSW_NX, GSW_NY, GSW_NZ = _GSW_NS[:]
    _read_saar()
    GSW_LONS, GSW_LATS, GSW_PS, GSW_NDEPTHS, SAAR_REF = _GSW_GRID
    DLON = GSW_LONS[1] - GSW_LONS[0]
    DLAT = GSW_LATS[1] - GSW_LATS[0]
    
    # Find location on lon/lat/pressure grid
    indx = int(lon1/DLON)
    indy = int((lat0-GSW_LATS[0])/DLAT)
    nzmax = max(GSW_NDEPTHS[indy+dj][indx+di] for (dj,di) in _DTAS) - 1
    if nzmax < 0:
        warnmsg = ('Given point ({0} E, {1} N) appears to be over '
            'land').format(lon0,lat0)
        warnings.warn(warnmsg,RuntimeWarning)
        return GSW_ERRVAL
    pdbar = min(pdbar,GSW_PS[nzmax])
    indz = _getind(GSW_PS,pdbar)
    
    # Interpolate near this point for a salinity upper limit
    r1 = (lon0 - GSW_LONS[indx]) / DLON
    s1 = (lat0 - GSW_LATS[indy]) / DLAT
    t1 = (pdbar - GSW_PS[indz]) / (GSW_PS[indz+1] - GSW_PS[indz])
    saars = [SAAR_REF[indz][indy+dj][indx+di] for (dj,di) in _DTAS]
    # Special treatment of the Panama barrier and out-of-bounds
    PLONS, PLATS = _PANAMA[:]
    if (lon0 > PLONS[0] and lon0 < PLONS[-1]
            and lat0 > PLATS[0] and lat0 < PLATS[-1]):
        saars = _gsw_addbarrier(saars,lon0,lat0,GSW_LONS[indx],GSW_LATS[indy],
            dlon,dlat)
    elif max(saars) >= _ERRVAL:
        saars = _gsw_addmean(saars,lon0,lat0)
    saar_upper = ((1-s1)*((1-r1)*saars[0] + r1*saars[1])
        + s1*((1-r1)*saars[3] + r1*saars[2]))
    
    # Get salinity lower bound
    saars = [SAAR_REF[indz+1][indy+dj][indx+di] for (dj,di) in _DTAS]
    if (lon0 > PLONS[0] and lon0 < PLONS[-1]
            and lat0 > PLATS[0] and lat0 < PLATS[-1]):
        saars = _gsw_addbarrier(saars,lon0,lat0,GSW_LONS[indx],GSW_LATS[indy],
            dlon,dlat)
    elif max(saars) >= _ERRVAL:
        saars = gsw_addmean(saars,lon0,lat0)
    saar_lower = ((1-s1)*((1-r1)*saars[0] + r1*saars[1])
        + s1*((1-r1)*saars[3] + r1*saars[2]))
    if saar_lower >= _ERRVAL:
        saar_lower = saar_upper
    
    # Calculate SAAR value by averaging
    saar = (1-t1)*saar_upper + t1*saar_lower
    if saar >= _ERRVAL:
        saar = GSW_ERRVAL
    return saar


## Public practical salinity unit functions
def gsw_safromsp(spsu,pdbar,lon0,lat0):
    """Convert practical to absolute salinity in GSW units.
    
    Convert from salinity in practical salinity units (psu) to absolute
    salinity in g/kg (the units used by the Gibbs Seawater toolbox).
    
    :arg float spsu: Salinity in practical salinity units (psu).
    :arg float pdbar: Seawater (gauge) pressure in dbar.
    :arg float lon0: Reference longitude in degrees East.
    :arg float lat0: Reference latitude in degrees North.
    :returns: Absolute salinity in g/kg.
    :raises RuntimeWarning: If the reference point is not over the
        ocean, as determined by the GSW_FNAME gridded data file.
        GSW_ERRVAL is returned.
    """
    # Handle Baltic sea locations first
    isbaltic = _check_baltic(lon0,lat0)
    if isbaltic:
        sgkg = gsw_safromspbaltic(spsu,lon0,lat0)
        if sgkg < _ERRVAL:
            return sgkg
    
    # Calculate salinity from SAAR value
    saar = _gsw_saar(pdbar,lon0,lat0)
    if saar >= _ERRVAL:
        return GSW_ERRVAL
    sgkg = _SAL1/_SAL0 * spsu * (1+saar)
    return sgkg

def gsw_spfromsa(sgkg,pdbar,lon0,lat0):
    """Convert absolute to practical salinity in GSW units.
    
    Convert from absolute salinity in g/kg (the units used by the Gibbs
    Seawater toolbox) to practical salinity units (psu).
    
    :arg float sgkg: Absolute salinity in g/kg.
    :arg float pdbar: Seawater (gauge) pressure in dbar.
    :arg float lon0: Reference longitude in degrees East.
    :arg float lat0: Reference latitude in degrees North.
    :returns: Salinity in practical salinity units (psu).
    :raises RuntimeWarning: If the reference point is not over the
        ocean, as determined by the GSW_FNAME gridded data file.
        GSW_ERRVAL is returned.
    """
    # Handle Baltic sea locations first
    isbaltic = _check_baltic(lon0,lat0)
    if isbaltic:
        spsu = gsw_safromspbaltic(sgkg,lon0,lat0)
        if spsu < _ERRVAL:
            return spsu
    
    # Calculate salinity from SAAR value
    saar = _gsw_saar(pdbar,lon0,lat0)
    if saar >= _ERRVAL:
        return GSW_ERRVAL
    spsu = _SAL0/_SAL1 * sgkg / (1+saar)
    return spsu

def sal_asalfrompsal(spsu,lon0,lat0,pres):
    """Convert practical to absolute salinity in SI units.
    
    Convert from salinity in practical salinity units (psu) to absolute
    salinity in kg/kg (the SI standard mass fraction unit).
    
    :arg float spsu: Salinity in practical salinity units (psu).
    :arg float lon0: Reference longitude in degrees East.
    :arg float lat0: Reference latitude in degrees North.
    :arg float pres: Absolute pressure (including atmosphere) in Pa.
    :returns: Absolute salinity in kg/kg.
    :raises RuntimeWarning: If the reference point is not over the
        ocean, as determined by the GSW_FNAME gridded data file.
        GSW_ERRVAL is returned.
    
    :Examples:
    
    >>> patm = 101325.
    >>> sal_asalfrompsal(35.7,201.,-21.,patm+1023e4)
    0.035873322343341174
    >>> sal_asalfrompsal(35.,180.,40.,patm+2e7)
    0.03518867406301395
    >>> sal_asalfrompsal(8.,20.,57.,patm)
    0.008104837714285716
    """
    pdbar = (pres - _PATM)/_DBAR2PA
    sgkg = gsw_safromsp(spsu,pdbar,lon0,lat0)
    if sgkg >= _ERRVAL:
        return GSW_ERRVAL
    salt = sgkg/_KG2G
    return salt

def sal_psalfromasal(salt,lon0,lat0,pres):
    """Convert absolute to practical salinity in SI units.
    
    Convert from absolute salinity in kg/kg (the SI standard mass
    fraction unit) to practical salinity units (psu).
    
    :arg float salt: Absolute salinity in kg/kg.
    :arg float lon0: Reference longitude in degrees East.
    :arg float lat0: Reference latitude in degrees North.
    :arg float pres: Absolute pressure (including atmosphere) in Pa.
    :returns: Salinity in practical salinity units (psu).
    :raises RuntimeWarning: If the reference point is not over the
        ocean, as determined by the GSW_FNAME gridded data file.
        GSW_ERRVAL is returned.
    
    :Examples:
    
    >>> patm = 101325.
    >>> asal = 0.035873322343341172
    >>> sal_psalfromasal(asal,201.,-21.,patm+1023e4)
    35.7
    """
    pdbar = (pres - _PATM)/_DBAR2PA
    sgkg = salt * _KG2G
    spsu = gsw_spfromsa(sgkg,pdbar,lon0,lat0)
    return spsu


