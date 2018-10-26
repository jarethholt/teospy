"""Conversions between SI and non-SI units.

This module provides functions to convert between SI units, which are
used throughout this library, and other common measurement scales.

:Examples:

>>> cnvpressure(4e3,'dbar','mpa')
40.101325
>>> cnvpressure(1023.,'dbar','pa')
1.0331325e7
>>> cnvpressure(1.0331325e7,'pa','torr')
7.74913101406e4
>>> cnvpressure(1.05350196040e5,'kgf','pa')
1.0331325e7
>>> cnvpressure(1.0331325e7,'pa','psi')
1498.42272437
>>> cnvpressure(1023.,'dbar','m')
1011.94563591
>>> cnvtemperature(300.,'k(t48)','degf(t48)')
80.33
>>> cnvtemperature(300.,'k(t48)','degc(t68)')
26.8411099988
>>> cnvtemperature(300.,'k(t48)','k(t90)')
299.984373916
>>> cnvtemperature(299.991108231,'k(t68)','k(t90)')
299.984372149
>>> sal78fromcnd(1.8880911555787682,40.,1e4)
40.
>>> cndfromsal78(40.,40.,1e4)
1.8880911555787682
>>> t0 = 273.15 + 25.5
>>> unitin = 'kg/kg(abs)'
>>> kwargs = {'p0': 101325. + 1023e4, 'lon0': 201., 'lat0': -21.}
>>> cnvsalinity(0.0357,unitin,'cl',**kwargs)
19.6659461767
>>> cnvsalinity(0.0357,unitin,'kn',**kwargs)
35.5270328489
>>> cnvsalinity(0.0357,unitin,'cnd',t0=t0,**kwargs)
1.27556269128
>>> cnvsalinity(0.0357,unitin,'pss',**kwargs)
35.5275150654
>>> cnvsalinity(0.0357,unitin,'kg/kg(ref)',**kwargs)
0.0356950425250
>>> cnvsalinity(8.,'pss',unitin,lon0=20.,lat0=57.)
0.008104837714285714

:Functions:

* :func:`swpres`: Approximate seawater pressure from depth.
* :func:`swdepth`: Approximate seawater depth from pressure.
* :func:`cnvpressure`: Convert between pressure scales.
* :func:`cnvtemperature`: Convert between temperature scales.
* :func:`sal78fromcnd`: Calculate salinity in the PSS-78 scale from
  conductivity, temperature, and pressure.
* :func:`cndfromsal78`: Calculate conductivity from salinity in the
  PSS-78 scale, temperature, and pressure.
* :func:`cnvsalinity`: Convert between salinity scales.

"""

__all__ = ['swpres','swdepth','cnvpressure','cnvtemperature',
    'sal78fromcnd','cndfromsal78','cnvsalinity']

import numpy
import warnings
import constants0
import convert0

_PATM = constants0.SEALEVEL_PRESSURE_SI
_TCELS = constants0.CELSIUS_TEMPERATURE_SI
_SAL1 = constants0.SO_SALINITY_SI
_SAL0 = constants0.SAL0
_GRAV = 9.80665
_MMHG = _PATM / 760.
_PSI = 6894.8
_C_FAHR = (32., 1.8)  # Constants related to Fahrenheit temperatures
_UPS = _SAL1 / 35.  # Reference salinity ratio
_UCL = 1.80655*_UPS  # Cl-normalized salinity ratio
_C35 = 4.2914  # Conductivity of reference seawater
_DBAR2PA = 1e4
_G2KG = 1e-3
_C_KNS = (0.03,1.805)
_CND2MMHO = 10.
_CNDMIN = 5e-4
_SPSUMIN = 0.2
_TSAL0 = 15.
_SALTOL = 1e-12
_ITMAX = 20
_LAT0 = 45.
_LON0 = -20.
_TSAL1 = 298.15

_PUNITS = ('pa', 'hpa', 'mpa', 'kpa', 'dbar', 'mbar', 'bar', 'kbar', 'torr',
    'mmhg', 'kgf', 'atm', 'lbf/in2', 'psi', 'm')
_PCONV2PA = {'pa': 1e0, 'hpa': 1e2, 'mpa': 1e6, 'kpa': 1e3, 'dbar': 1e4,
    'mbar': 1e2, 'bar': 1e5, 'kbar': 1e8, 'torr': _MMHG, 'mmhg': _MMHG,
    'kgf': _GRAV*10, 'atm': _PATM, 'lbf/in2': _PSI, 'psi': _PSI}
_TUNITS = ('degf(t48)', 'degf(t68)', 'degf(t90)', 'degc(t48)', 'degc(t68)',
    'degc(t90)', 'k(t48)', 'k(t68)', 'k(t90)', 'f', 'c', 'k')
_SUNITS = ('cl','kn','psu','pss','one','kg/kg(ref)','g/kg(ref)','kg/kg(abs)',
    'g/kg(abs)','ms/cm','mmho/cm','s/m','cnd')

_T48MIN = 90.  # Minimum temperature for the 48-scale
_T48MAX = 872.49268  # Temperature at which parameterization fails
_T48TOL = 1e-12
_C_4868 = ((273.15,250.846209678803,135.099869965,52.78567590085,
        27.6768548854105,39.1053205376684,65.5613230578069,80.8035868559867,
        70.5242118234052,44.7847589638966,21.2525653556078,7.67976358170846,
        2.1368945938285,0.459843348928069,7.63614629231648e-2,
        9.69328620373121e-3,9.23069154007008e-4,6.38116590952654e-5,
        3.02293237874619e-6,8.7755139130376e-8,1.17702613125477e-9),
    (1.,3.984517e-3,-5.85502e-7,4.35716e-10,-4.35716e-12),
    (0.003984712,-1.174912e-6), (-1.9539e-7,1.9539e-9),
    (4.5e-4,100.,419.58,630.74,-5.87456e-7))
_T68LIMS = (13.8,83.8,903.89,1337.58)
_T68SCLS = (40.,630.,1337.33)
_C_6890 = (
    (-0.005903,0.008174,-0.061924,-0.193388,1.490793,1.252347,-9.835868,
        1.411912,25.277595,-19.183815,-18.437089,27.000895,-8.716324),
    (-0.148759,-0.267408,1.08076,1.269056,-4.089591,-1.871251,7.438081,
        -3.536296),
    (78.687209,-0.47135991,1.0954715e-3,-1.2357884e-6,6.7736583e-10,
        -1.4458081e-13))
_T90LIMS = (13.8,83.8,903.765,1337.33)
_C_9068 = (
    (-0.005903,0.008174,-0.061924,-0.193388,1.490793,1.252347,-9.835868,
        1.411912,25.277595,-19.183815,-18.437089,27.000895,-8.716324),
    (-0.148759,-0.267408,1.08076,1.269056,-4.089591,-1.871251,7.438081,
        -3.536296),
    (78.687209,-0.47135991,1.0954715e-3,-1.2357884e-6,6.7736583e-10,
        -1.4458081e-13))

_C_AFUN = (-3.107e-3,0.4215)
_C_BFUN = (4.464e-4,3.426e-2,1.)
_C_CFUN = (3.989e-15,-6.370e-10,2.070e-5,0.)
_C_RT35 = (1.0031e-9,-6.9698e-7,1.104259e-4,2.00564e-2,0.6766097)
_C_SFUN = (0.0162, (2.7081,-7.0261,14.0941,25.3851,-0.1692,0.0080),
    (-0.0144,0.0636,-0.0375,-0.0066,-0.0056,0.0005))
_C_SWD = (5.92e-3,5.25e-3,4.42e-6)


## Pressure conversion
def swpres(depth,lat=_LAT0):
    """Approximate seawater pressure from depth.
    
    Approximate the seawater pressure at a given depth, based on an
    empirical function of latitude.
    
    :arg float depth: Depth in meters.
    :arg float lat: Latitude in degrees North (default _LAT0).
    :returns: Seawater pressure in decibars.
    :raises ValueError: If the depth is negative.
    :raises ValueError: If the latitude is not between -90 and 90.
    """
    if depth < 0:
        errmsg = 'Depth must be nonnegative.'
        raise ValueError(errmsg)
    if abs(lat) > 90:
        errmsg = 'Latitude must be between -90 and 90.'
        raise ValueError(errmsg)
    
    x = numpy.sin(numpy.abs(lat) * numpy.pi/180)
    c1 = _C_SWD[0] + _C_SWD[1]*x**2
    pdbar = (1-c1 - ((1-c1)**2 - 2*_C_SWD[2]*depth)**.5)/_C_SWD[2]
    return pdbar

def swdepth(pdbar,lat=_LAT0):
    """Approximate seawater depth from pressure.
    
    Approximate the seawater depth at a given pressure, based on an
    empirical function of latitude.
    
    :arg float pdbar: Seawater pressure in decibars.
    :arg float lat: Latitude in degrees North (default _LAT0).
    :returns: Depth in meters.
    :raises ValueError: If the pressure is negative.
    :raises ValueError: If the latitude is not between -90 and 90.
    """
    if pdbar < 0:
        errmsg = 'Pressure must be nonnegative.'
        raise ValueError(errmsg)
    if abs(lat) > 90:
        errmsg = 'Latitude must be between -90 and 90.'
        raise ValueError(errmsg)
    
    x = numpy.sin(numpy.abs(lat) * numpy.pi/180)
    c1 = _C_SWD[0] + _C_SWD[1]*x**2
    depth = ((1-c1)**2 - (1-c1 - pdbar*_C_SWD[2])**2)/(2*_C_SWD[2])
    return depth

def cnvpressure(presin,unitin,unitout,lat=_LAT0):
    """Convert between pressure units.
    
    Convert pressure from one unit system to another.
    
    :arg float presin: Pressure in the input units.
    :arg str unitin: Units of the input pressure; see _PUNITS for
        options. The strings are not case-sensitive.
    :arg str unitout: Units to convert to.
    :arg float lat: Latitude in degrees North (default _LAT0).
    :returns: Pressure in the output units.
    :raises ValueError: If either unitin or unitout are not recognized
        units.
    
    :Examples:
    
    >>> cnvpressure(4e3,'dbar','mpa')
    40.101325
    >>> cnvpressure(1.0331325e7,'pa','dbar')
    1023.
    >>> cnvpressure(1023.,'dbar','pa')
    1.03313250000e7
    >>> cnvpressure(1.0331325e7,'pa','torr')
    7.74913101406e4
    >>> cnvpressure(7.74913101406e4,'torr','pa')
    1.03313250000e7
    >>> cnvpressure(1.0331325e7,'pa','kgf')
    1.05350196040e5
    >>> cnvpressure(1.05350196040e5,'kgf','pa')
    1.03313250000e7
    >>> cnvpressure(1.0331325e7,'pa','atm')
    101.9622501850481
    >>> cnvpressure(1.05350196040e5,'atm','pa')
    10674608613.753
    >>> cnvpressure(1.0331325e7,'pa','lbf/in2')
    1498.42272437
    >>> cnvpressure(1.49842272437e3,'lbf/in2','pa')
    1.03313250000e7
    >>> cnvpressure(1.0331325e7,'pa','psi')
    1498.42272437
    >>> cnvpressure(1.49842272437e3,'psi','pa')
    1.03313250000e7
    >>> cnvpressure(1023.,'dbar','m')
    1011.94563591
    >>> cnvpressure(1011.94563591,'m','dbar')
    1023.
    """
    unitin = unitin.lower()
    unitout = unitout.lower()
    if unitin not in _PUNITS:
        errmsg = 'Input unit {0} must be one of _PUNITS'.format(unitin)
        raise ValueError(errmsg)
    if unitout not in _PUNITS:
        errmsg = ('Output unit {0} must be one of _PUNITS').format(unitout)
        raise ValueError(errmsg)
    
    if unitin == 'm':
        pres = _PATM + _DBAR2PA*swpres(presin,lat=lat)
    else:
        cnv = _PCONV2PA[unitin]
        pres = presin*cnv
        if unitin == 'dbar':
            pres += _PATM
    
    if unitout == 'm':
        presout = swdepth((pres-_PATM)/_DBAR2PA,lat=lat)
    else:
        cnv = _PCONV2PA[unitout]
        if unitout == 'dbar':
            pres -= _PATM
        presout = pres/cnv
    return presout


## Temperature scale functions
def _t68fromt48(t48,tol=_T48TOL,itmax=_ITMAX):
    """Convert IPTS-48 temperature to IPTS-68.
    
    Convert between temperatures based on the IPTS-48 and IPTS-68
    reference systems.
    
    :arg float t48: IPTS-48 temperature in K.
    :arg float tol: Tolerance for halting the iterative calculation
        (default _T48TOL).
    :arg int itmax: Maximum number of allowed iterations (default
        _ITMAX).
    :returns: IPTS-68 temperature in K.
    :raises RuntimeWarning: If t48 is lower than _T48MIN or higher than
        _T48MAX.
    :raises RuntimeWarning: If the iterative calculation reaches the
        maximum number of iterations before reaching the tolerance.
    """
    if t48 < _T48MIN:
        warnmsg = ('No conversion available for temperatures < '
            '{1} K').format(_T48MIN)
        warnings.warn(warnmsg,RuntimeWarning)
        return t48
    if t48 > _T48MAX:
        warnmsg = ('Conversion calculation fails for temperatures > '
            '{0} K').format(_T48MAX)
        warnings.warn(warnmsg,RuntimeWarning)
        return t48
    
    if t48 < _TCELS:
        x = t48 - _TCELS
        larg = 0.
        for coeff in _C_4868[1][::-1]:
            larg = larg*x + coeff
        hg0 = numpy.log(larg)
        t68 = 0.
        for coeff in _C_4868[0][::-1]:
            t68 = t68*hg0 + coeff
        return t68
    
    x = t48 - _TCELS
    hg1 = _C_4868[2][0] + _C_4868[2][1]*x
    hg2 = _C_4868[3][0]*x + _C_4868[3][1]*x**2
    x1 = x
    x2 = x
    for __ in range(itmax):
        fkt = _C_4868[4][0] * x2
        for coeff in _C_4868[4][1:-1]:
            fkt *= (x2/coeff - 1)
        hg3 = _C_4868[4][-1] * (x1 - x - fkt)**2
        x1 = x + fkt + (hg2-hg3)/hg1
        if (abs(x1 - x2) < tol):
            break
        x2 = x1
    else:
        warnmsg = ('Maximum number of iterations {0} reached before step size '
            '{1} is within tolerance {2}').format(itmax,abs(x1-x2),tol)
        warnings.warn(warnmsg,RuntimeWarning)
    t68 = x1 + _TCELS
    return t68

def _t90fromt68(t68):
    """Convert IPTS-68 temperature to ITS-90.
    
    Convert between temperatures based on the IPTS-68 and ITS-90
    reference systems.
    
    :arg float t68: IPTS-68 temperature in K.
    :returns: ITS-90 temperature in K.
    :raises RuntimeWarning: If t68 is lower than _T68LIMS[0].
    """
    deltat = 0.
    if t68 < _T68LIMS[0]:
        warnmsg = ('No conversion available for temperatures < {0} '
            'K').format(_T68LIMS[0])
        warnings.warn(warnmsg,RuntimeWarning)
    elif t68 >= _T68LIMS[3]:
        deltat = -(t68/_T68SCLS[2])**2/4
    else:
        if t68 < _T68LIMS[1]:
            ind = 0
            x = t68/_T68SCLS[0] - 1
        elif t68 < _T68LIMS[2]:
            ind = 1
            x = (t68 - _TCELS) / _T68SCLS[1]
        else:
            ind = 2
            x = t68 - _TCELS
        for coeff in _C_6890[ind][::-1]:
            deltat = x*(deltat + coeff)
    t90 = t68 + deltat
    return t90

def _t68fromt90(t90):
    """Convert ITS-90 temperature to IPTS-68.
    
    Convert between temperatures based on the ITS-90 and IPTS-68
    reference systems.
    
    :arg float t90: ITS-90 temperature in K.
    :returns: IPTS-68 temperature in K.
    :raises RuntimeWarning: If t90 is lower than _T90LIMS[0].
    """
    deltat = 0.
    if t90 < _T90LIMS[0]:
        warnmsg = ('No conversion available for temperatures < {0} '
            'K').format(_T90LIMS[0])
        warnings.warn(warnmsg,RuntimeWarning)
    elif t90 >= _T90LIMS[3]:
        deltat = -(t90/_T68SCLS[2])**2/4
    else:
        if t90 < _T90LIMS[1]:
            ind = 0
            x = t90/_T68SCLS[0] - 1
        elif t90 < _T90LIMS[2]:
            ind = 1
            x = (t90 - _TCELS)/_T68SCLS[1]
        else:
            ind = 2
            x = t90 - _TCELS
        for coeff in _C_9068[ind][::-1]:
            deltat = x*(deltat + coeff)
    t68 = t90 - deltat
    return t68

def cnvtemperature(tempin,unitin,unitout,tol=_T48TOL,itmax=_ITMAX):
    """Convert between temperature units.
    
    Convert temperature from one system of units to another.
    
    :arg float tempin: Temperature in the input units.
    :arg str unitin: Units of the input temperature; see _TUNITS for
        options. The strings are not case-sensitive.
    :arg str unitout: Units to convert to.
    :arg float tol: Tolerance for halting the iterative calculation
        (default _T48TOL). Only used for IPTS-48.
    :arg int itmax: Maximum number of allowed iterations (default
        _ITMAX). Only used for IPTS-48.
    :returns: Temperature in the output units.
    :raises ValueError: If either unitin or unitout are not recognized
        units.
    :raises RuntimeWarning: If trying to convert to the IPTS-48 scale
        from any other scale (IPTS-68 or ITS-90).
    
    :Examples:
    
    >>> cnvtemperature(300.,'k(t48)','degf(t48)')
    80.33
    >>> cnvtemperature(300.,'k(t48)','degf(t68)')
    80.3139979978
    >>> cnvtemperature(300.,'k(t48)','degf(t90)')
    80.3018730496
    >>> cnvtemperature(300.,'k(t48)','degc(t48)')
    26.8500000000
    >>> cnvtemperature(300.,'k(t48)','degc(t68)')
    26.8411099988
    >>> cnvtemperature(300.,'k(t48)','degc(t90)')
    26.8343739165
    >>> cnvtemperature(300.,'k(t48)','k(t68)')
    299.991109999
    >>> cnvtemperature(300.,'k(t48)','k(t90)')
    299.984373916
    >>> cnvtemperature(299.991108231,'k(t68)','degf(t68)')
    80.3139948158
    >>> cnvtemperature(299.991108231,'k(t68)','degf(t90)')
    80.3018698685
    >>> cnvtemperature(299.991108231,'k(t68)','degc(t68)')
    26.8411082310
    >>> cnvtemperature(299.991108231,'k(t68)','degc(t90)')
    26.8343721491
    >>> cnvtemperature(299.991108231,'k(t68)','k(t90)')
    299.984372149
    >>> cnvtemperature(299.984372149,'k(t90)','degf(t90)')
    80.3018698682
    >>> cnvtemperature(299.984373917,'k(t90)','degc(t90)')
    26.8343739170
    """
    unitin = unitin.lower()
    unitout = unitout.lower()
    if unitin not in _TUNITS:
        errmsg = 'Input unit {0} is not in _TUNITS'.format(unitin)
        raise ValueError(errmsg)
    if unitout not in _TUNITS:
        errmsg = 'Output unit {0} is not in in _TUNITS'.format(unitout)
        raise ValueError(errmsg)
    
    if 't48' in unitout:
        if 't48' not in unitin:
            warnmsg = ('Cannot convert to IPTS-48 temperatures unless the '
                'input units are also IPTS-48. Only the standard conversion '
                'will be applied')
            warnings.warn(warnmsg,RuntimeWarning)
        if 'f' in unitin:
            tempin_k = _TCELS + (tempin - _C_FAHR[0])/_C_FAHR[1]
        elif 'c' in unitin:
            tempin_k = _TCELS + tempin
        else:
            tempin_k = tempin
        if 'f' in unitout:
            tempout = (tempin_k - _TCELS)*_C_FAHR[1] + _C_FAHR[0]
        elif 'c' in unitout:
            tempout = tempin_k - _TCELS
        else:
            tempout = tempin_k
        return tempout
    
    if 'f' in unitin:
        tempin_k = _TCELS + (tempin - _C_FAHR[0])/_C_FAHR[1]
    elif 'c' in unitin:
        tempin_k = _TCELS + tempin
    else:
        tempin_k = tempin
    
    if 't48' in unitin:
        tempout_k = _t68fromt48(tempin_k,tol=tol,itmax=itmax)
        if 't90' in unitout:
            tempout_k = _t90fromt68(tempout_k)
    elif 't68' in unitin:
        if 't68' in unitout:
            tempout_k = tempin_k
        else:
            tempout_k = _t90fromt68(tempin_k)
    else:
        if 't68' in unitout:
            tempout_k = _t68fromt90(tempin_k)
        else:
            tempout_k = tempin_k
    
    if 'f' in unitout:
        tempout = (tempout_k - _TCELS)*_C_FAHR[1] + _C_FAHR[0]
    elif 'c' in unitout:
        tempout = tempout_k - _TCELS
    else:
        tempout = tempout_k
    return tempout


## Salinity-conductivity functions
def _abcrtfun(t68,pdbar):
    """Calculate terms in the conductivity-salinity relation.
    
    Calculate terms in the conductivity-salinity relation. The first two
    are temperature-dependent; the third is pressure-dependent. The last
    term is the ratio of the conductivity at the given temperature to
    that at the reference temperature of 15 degrees Celsius.
    
    :arg float t68: Temperature in degrees Celsius (IPTS-68 scale).
    :arg float pdbar: Pressure in decibars.
    :returns: Temperature- and pressure-dependent terms in the
        conductivity-salinity relation (all unitless).
    :rtype: tuple(float*4)
    """
    a = 0.
    for coeff in _C_AFUN:
        a = a*t68 + coeff
    b = 0.
    for coeff in _C_BFUN:
        b = b*t68 + coeff
    c = 0.
    for coeff in _C_CFUN:
        c = c*pdbar + coeff
    rt0 = 0.
    for coeff in _C_RT35:
        rt0 = rt0*t68 + coeff
    return (a, b, c, rt0)

def _salfun(rt,dt):
    """Calculate salinity from conductivity and temperature variables.
    
    Calculate the salinity in the Practical Salinity Scale 1978 (PSS-78)
    from auxiliary variables related to conductivity and temperature.
    
    :arg float rt: Square root of the conductivity ratio, unitless.
    :arg float dt: Temperature difference in degrees Celsius (IPTS-
        69 scale) from the reference temperature (15 degrees).
    :returns: Salinity in PSS-78.
    """
    s1 = 0.
    for coeff in _C_SFUN[1]:
        s1 = s1*rt + coeff
    s2 = 0.
    for coeff in _C_SFUN[2]:
        s2 = s2*rt + coeff
    spsu = s1 + s2*dt/(1 + _C_SFUN[0]*dt)
    return spsu

def _dsalfun(rt,dt):
    """Calculate derivative of salfun with respect to rt.
    
    Calculate the derivative of the function salfun with respect to the
    conductivity variable rt.
    
    :arg float rt: Square root of the conductivity ratio, unitless.
    :arg float dt: Temperature difference in degrees Celsius (IPTS-
        69 scale) from the reference temperature (15 degrees).
    :returns: Derivative of salinity in PSS-78.
    """
    ds1 = 0.
    for (i,coeff) in enumerate(_C_SFUN[1][::-1]):
        if i == 0:
            continue
        ds1 += coeff * i*rt**(i-1)
    ds2 = 0.
    for (i,coeff) in enumerate(_C_SFUN[2][::-1]):
        if i == 0:
            continue
        ds2 += coeff * i*rt**(i-1)
    dspsu = ds1 + ds2*dt/(1 + _C_SFUN[0]*dt)
    return dspsu

def sal78fromcnd(cnd,t68,pdbar):
    """Convert conductivity to salinity.
    
    Convert conductivity to salinity, as measured in the PSS-78 scale.
    
    :arg float cnd: Conductivity ratio, unitless.
    :arg float t68: Temperature in degrees Celsius (IPTS-68 scale).
    :arg float pdbar: Pressure in decibars.
    :returns: Salinity in practical salinity units (PSS-78).
    :raises RuntimeWarning: If the conductivity ratio is less than
        _CNDMIN.
    
    :Examples:
    
    >>> sal78fromcnd(1.8880911555787682,40.,1e4)
    40.
    """
    if cnd <= _CNDMIN:
        warnmsg = ('Zero salinity is assumed for conductivity ratio '
            '<{0}').format(_CNDMIN)
        warnings.warn(warnmsg,RuntimeWarning)
        return 0.
    
    res = cnd
    dt = t68 - _TSAL0
    a, b, c, rt0 = _abcrtfun(t68,pdbar)
    rt = res / (rt0 * (1 + c/(a*res+b)))
    rt = numpy.abs(rt)**.5
    spsu = _salfun(rt,dt)
    return spsu

def cndfromsal78(spsu,t68,pdbar,tol=_SALTOL,itmax=_ITMAX):
    """Convert salinity to conductivity.
    
    Convert salinity to conductivity, as measured in the PSS-78
    reference scale.
    
    :arg float spsu: Salinity in practical salinity units (PSS-78).
    :arg float t68: Temperature in degrees Celsius (IPTS-68 scale).
    :arg float pdbar: Pressure in decibars.
    :arg float tol: Tolerance for halting the iterative calculation
        (default _SALTOL).
    :arg int itmax: Maximum number of allowed iterations (default
        _ITMAX).
    :returns: Conductivity ratio, unitless.
    :raises RuntimeWarning: If the salinity is less than _SPSUMIN.
    
    :Examples:
    
    >>> cndfromsal78(40.,40.,1e4)
    1.8880911555787682
    """
    if spsu <= _SPSUMIN:
        warnmsg = ('Zero conductivity is assumed for salinity < '
            '{0} PSU').format(_SPSUMIN)
        warnings.warn(warnmsg,RuntimeWarning)
        return 0.
    
    dt = t68 - _TSAL0
    rt = (spsu/_SAL0)**.5
    si = _salfun(rt,dt)
    for it in range(itmax):
        rt = rt + (spsu-si) / _dsalfun(rt,dt)
        si = _salfun(rt,dt)
        ds = abs(si - spsu)
        if (it > 0) and (ds < tol):
            break
    else:
        warnmsg = ('Maximum number of iterations {0} reached before step size '
            '{1} is within tolerance {2}').format(itmax,ds,tol)
        warnings.warn(warnmsg,RuntimeWarning)
    
    a, b, c, rt0 = _abcrtfun(t68,pdbar)
    rtt = rt0 * rt**2
    c1 = rtt * (c + b)
    b1 = b - rtt*a
    res = (abs(b1**2 + 4*a*c1))**.5 - b1
    cnd = .5 * res / a
    return cnd

def cnvsalinity(saltin,unitin,unitout,t0=_TSAL1,p0=_PATM,lon0=_LON0,
    lat0=_LAT0,tol=_SALTOL,itmax=_ITMAX):
    """Convert between salinity units.
    
    Convert salinity from one unit system to another.
    
    :arg float saltin: Salinity in the input units.
    :arg str unitin: Units of the input salinity; see _SUNITS for
        options. The strings are not case-sensitive.
    :arg str unitout: Units to convert to.
    :arg float t0: Reference temperature in K (ITS-90 scale; default
        _TSAL1). Only used when converting conductivity.
    :arg float p0: Reference pressure in Pa (default _PATM). Used when
        converting conductivity or practical salinity.
    :arg float lon0: Reference longitude in degrees East. Default
        location (_LON0,_LAT) is in northeastern Atlantic. Only used
        when converting practical salinity.
    :arg float lat0: Reference latitude in degrees North.
    :arg float tol: Tolerance for halting the iterative calculation
        (default _SALTOL). Only used when converting from practical
        salinity to conductivity.
    :arg int itmax: Maximum number of allowed iterations (default
        _ITMAX). Only used when converting from practical salinity to
        conductivity.
    :returns: Salinity in the output units.
    :raises ValueError: If either unitin or unitout are not recognized
        units.
    
    :Examples:
    
    >>> p0 = 101325. + 1023e4
    >>> t0 = 273.15 + 25.5
    >>> unitin = 'kg/kg(abs)'
    >>> lon0, lat0 = 201., -21.
    >>> kwargs = {'p0': p0, 'lon0': lon0, 'lat0': lat0}
    >>> cnvsalinity(0.0357,unitin,'cl',**kwargs)
    19.6659461767
    >>> cnvsalinity(0.0357,unitin,'kn',**kwargs)
    35.5270328489
    >>> cnvsalinity(0.0357,unitin,'cnd',t0=t0,**kwargs)
    1.27556269128
    >>> cnvsalinity(0.0357,unitin,'pss',**kwargs)
    35.5275150654
    >>> cnvsalinity(0.0357,unitin,'kg/kg(ref)',**kwargs)
    0.0356950425250
    >>> cnvsalinity(0.0356951724471,'kg/kg(ref)',unitin,**kwargs)
    0.0357
    >>> cnvsalinity(35.52764437773386,'pss','cl')
    19.6660177563
    >>> cnvsalinity(35.52764437773386,'pss','kn')
    35.5271620502
    >>> cnvsalinity(35.52764437773386,'pss','cnd',t0=t0,p0=p0)
    1.27556680822
    >>> cnvsalinity(35.52764437773386,'pss','kg/kg(ref)')
    0.0356951724471
    >>> cnvsalinity(35.52764437773386,'pss','kg/kg(abs)',**kwargs)
    0.0357
    >>> cnvsalinity(35.,'pss',unitin,p0=101325+2e7,lon0=179.,lat0=40.)
    0.0351890192696528
    >>> cnvsalinity(8.,'pss',unitin,lon0=20.,lat0=57.)
    0.008104837714285714
    >>> cnvsalinity(35.5271620502,'kn','pss')
    35.5276443777
    >>> cnvsalinity(35.5271620502,'kn','kg/kg(abs)',**kwargs)
    0.0357
    >>> cnvsalinity(1.27556680822,'cnd','pss',t0=t0,p0=p0)
    35.5276443779
    >>> cnvsalinity(1.27556680822,'cnd',unitin,t0=t0,**kwargs)
    0.0357
    """
    unitin = unitin.lower()
    unitout = unitout.lower()
    if unitin not in _SUNITS:
        errmsg = 'Input unit {0} must be one of _SUNITS'.format(unitin)
        raise ValueError(errmsg)
    if unitout not in _SUNITS:
        errmsg = ('Output unit {0} must be one of _SUNITS').format(unitout)
        raise ValueError(errmsg)
    tc = t0 - _TCELS
    pb = (p0 - _PATM)/_DBAR2PA
    
    if unitin == 'cl':
        salt = saltin * _UCL
    elif unitin == 'kn':
        salt = (saltin - _C_KNS[0]) / _C_KNS[1] * _UCL
    elif unitin in ('psu','pss','one'):
        salt = saltin * _UPS
    elif unitin == 'kg/kg(ref)':
        salt = saltin
    elif unitin == 'g/kg(ref)':
        salt = saltin * _G2KG
    elif unitin == 'kg/kg(abs)':
        salt = convert0.sal_psalfromasal(saltin,lon0,lat0,p0)*_UPS
    elif unitin == 'g/kg(abs)':
        salt = convert0.sal_psalfromasal(saltin*_G2KG,lon0,lat0,p0)*_UPS
    else:
        # Various forms of conductivity
        t68c = cnvtemperature(t0,'k(t90)','degc(t68)')
        pdb = cnvpressure(p0,'pa','dbar')
        if unitin in ('ms/cm','mmho/cm'):
            cratio = saltin / _CND2MMHO / _C35
        elif unitin == 's/m':
            cratio = saltin / _C35
        else:
            cratio = saltin
        salt = sal78fromcnd(cratio,t68c,pdb) * _UPS
    
    if unitout == 'cl':
        saltout = salt / _UCL
    elif unitout == 'kn':
        saltout = _C_KNS[0] + _C_KNS[1]*salt/_UCL
    elif unitout in ('psu','pss','one'):
        saltout = salt / _UPS
    elif unitout == 'kg/kg(ref)':
        saltout = salt
    elif unitout == 'g/kg(ref)':
        saltout = salt / _G2KG
    elif unitout == 'kg/kg(abs)':
        saltout = convert0.sal_asalfrompsal(salt/_UPS,lon0,lat0,p0)
    elif unitout == 'g/kg(abs)':
        saltout = convert0.sal_asalfrompsal(salt/_UPS,lon0,lat0,p0)/_G2KG
    else:
        # Various forms of conductivity
        t68c = cnvtemperature(t0,'k(t90)','degc(t68)')
        pdb = cnvpressure(p0,'pa','dbar')
        psal = salt / _UPS
        saltout = cndfromsal78(psal,t68c,pdb)
        if unitout in ('ms/cm','mmho/cm'):
            saltout *= _CND2MMHO*_C35
        elif unitout == 's/m':
            saltout *= _C35
    return saltout

