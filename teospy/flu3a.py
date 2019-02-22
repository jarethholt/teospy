"""Fluid water Gibbs free energy.

This module implements the Gibbs free energies of both liquid water and water vapour and their derivatives with respect to temperature and pressure. It also provides several approximations for the equilibrium densities which can be used to initialize the equilibrium function.

:Examples:

>>> eq_tp_liq(300.,1e5)
996.556340389
>>> liq_g(0,0,300.,1e5)
-5265.05045577
>>> liq_g(0,1,300.,1e5)
1.00345555938e-03
>>> liq_g(2,0,300.,1e5)
-13.9354650734
>>> eq_tp_vap(300.,1e3)
7.22603510025e-03
>>> vap_g(0,0,300.,1e3)
-180090.341338
>>> vap_g(0,1,300.,1e3)
138.388478069
>>> vap_g(2,0,300.,1e3)
-6.24707163427

:Functions:

* :func:`eq_tp_liq`: Calculate equilibrium liquid water quantities at
  temperature and pressure.
* :func:`eq_tp_vap`: Calculate equilibrium water vapour quantities at
  temperature and pressure.
* :func:`liq_g`: Liquid water Gibbs free energy with derivatives.
* :func:`vap_g`: Water vapour Gibbs free energy with derivatives.

"""

__all__ = ['eq_tp_liq','eq_tp_vap','liq_g','vap_g']

import warnings
import numpy
from teospy import constants0
from teospy import flu1
from teospy import flu2
from teospy import maths3

_CHKTOL = constants0.CHKTOL
_RWAT = constants0.RWAT
_PATM = constants0.PATM
_TCELS = constants0.TCELS
_DCP = constants0.CP_DENSITY_SI
_TCP = constants0.CP_TEMPERATURE_SI
_PCP = constants0.CP_PRESSURE_SI

_chkflubnds = constants0.chkflubnds
_flu_f = flu1.flu_f
_eq_pressure = flu2.eq_pressure
_newton = maths3.newton

# Constants related to approximation functions
_LIMS = ((623.15, _TCP, 650., 1073.15), (1.6529e7, _PCP, 3.5e7, 1e8))
_C_E80 = ((1.00024, 1e5),
    (19652.21,148.4206,-2.327105,0.01360477,-5.155288e-5),
    (3.239908,0.00143713,1.16092e-4,-5.77905e-7),
    (8.50935e-5,-6.12293e-6,5.2787e-8),
    (999.842594,0.06793952,-0.00909529,1.001685e-4,-1.120083e-6,6.536332e-9))
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
_C_CRIT = ((0.0,-7.60041479494879,118.661872386874),
    (0.0,-17.463827264079,186.040087842884),
    (0.0,0.69701967809328,25.5059905941023),
    (-0.602044738250314,30.8633119943879,14.4873846518829))
_C_IF97L = ((1386., 1.653e7, 461.526, 1.222, 7.1),
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
_C_IF97V = ((540., 1e6, 461.526, 0.5),
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


## Density approximation functions
def _dliq_eos80(temp,pres):
    """Approximate liquid water density using EOS 1980.
    
    Approximate the density of liquid water using the Equation of
    State (1980).
    
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :returns: Liquid water density in kg/m3.
    """
    TSCL, PRED = _C_E80[0]
    t = TSCL * (temp - _TCELS)
    p = (pres - _PATM) / PRED
    
    pcoeffs = list()
    for coeffs in _C_E80[1:4]:
        pcoeff = 0.
        for (k,coeff) in enumerate(coeffs):
            pcoeff += coeff * t**k
        pcoeffs.append(pcoeff)
    denom = 1 - p / (pcoeffs[0] + p*pcoeffs[1] + p**2*pcoeffs[2])
    numer = 0.0
    for (k,coeff) in enumerate(_C_E80[4]):
        numer += coeff * t**k
    dliq = numer / denom
    return dliq

def _dliq_if97(temp,pres):
    """Approximate liquid water density using IF97.
    
    Approximate the density of liquid water using the IAPWS-IF97
    parameterization.
    
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :returns: Liquid water density in kg/m3.
    """
    if temp <= 0 or pres <= 0:
        errmsg = 'Temperature and pressure must be positive'
        raise ValueError(errmsg)
    
    TRED, PRED, RWAT, TAU0, PSI0 = _C_IF97L[0]
    tau = TRED / temp
    psi = pres / PRED
    rt = RWAT * temp
    tt = tau - TAU0
    pp = PSI0 - psi
    
    g1 = 0.
    for (i,j,n) in _C_IF97L[1]:
        if tt == 0:
            if j == 0:
                pwrt = 1.
            else:
                pwrt = 0.
        else:
            pwrt = tt**j
        if pp == 0:
            if i == 1:
                pwrp = 1.
            else:
                pwrp = 0.
        else:
            pwrp = -i * pp**(i - 1)
        g1 += n * pwrp * pwrt
    
    g_p = rt/PRED * g1
    dliq = g_p**(-1.)
    return dliq

def _dvap_if97(temp,pres):
    """Approximate water vapour density using IF97.
    
    Approximate the density of water vapour using the IAPWS-IF97
    parameterization.
    
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :returns: Water vapour density in kg/m3.
    """
    if temp <= 0 or pres <= 0:
        errmsg = 'Temperature and pressure must be positive'
        raise ValueError(errmsg)
    
    TRED, PRED, RWAT, TAU0 = _C_IF97V[0]
    tau = TRED / temp
    psi = pres / PRED
    rt = RWAT * temp
    
    g0 = psi**(-1)
    gr = 0.
    tt = tau - TAU0
    for (i,j,n) in _C_IF97V[2]:
        if tt == 0:
            if j == 0:
                pwrt = 1.
            else:
                pwrt = 0.
        else:
            pwrt = tt**j
        pwrp = i*psi**(i - 1)
        gr += n * pwrp * pwrt
    
    g_p = rt/PRED * (g0 + gr)
    dvap = g_p**(-1)
    return dvap

def _dliq_f03(temp,pres):
    """Approximate liquid water density using Feistel (2003).
    
    Approximate the density of liquid water using the Feistel (2003)
    parameterization.
    
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :returns: Liquid water density in kg/m3.
    """
    TRED, PRED = _C_F03[0]
    y = (temp - _TCELS) / TRED
    z = (pres - _PATM) / PRED
    
    g_p = 0.0
    for (j,k,c) in _C_F03[1:]:
        if k == 0:
            continue
        g_p += c * y**j * k*z**(k-1)/PRED
    dliq = g_p**(-1)
    return dliq

def _dvap_ideal(temp,pres):
    """Approximate water vapour density using ideal gas law.
    
    Approximate the density of water vapour using the ideal gas law.
    
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :returns: Water vapour density in kg/m3.
    """
    dvap = pres / (_RWAT * temp)
    return dvap

def _dliqvap_crit(temp,pres):
    """Approximate fluid water densities near critical point.
    
    Approximate the densities of fluid water near the critical point,
    based on a cubic equation of state. In subcritical conditions, both
    liquid water and water vapour densities are returned. In
    supercritical conditions, one of the returned values may be 0.
    
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :returns: Liquid water and water vapour densities in kg/m3.
    """
    pr = _PCP / pres - 1
    tr = temp / _TCP - 1
    
    coeffs = list()
    for (i,ccs) in enumerate(_C_CRIT):
        coeff = 0.
        for (k,cc) in enumerate(ccs):
            coeff += cc * tr**k
        coeffs.append(coeff)
    a0, a1, a2, a3 = coeffs
    
    if tr < 0:
        # Get the pressure range of 2-phase solutions
        r = a2**2 - 3*a1*a3
        r = r**.5
        d1 = -(a2+r) / (3*a3)
        d2 = -(a2-r) / (3*a3)
        p1 = a0 + d1*(a1 + d1*(a2 + d1*a3))
        p2 = a0 + d2*(a1 + d2*(a2 + d2*a3))
    else:
        p1 = pr
        p2 = pr
    r = a2 / a3
    s = a1 / a3
    t = (a0 - pr) / a3
    
    # Find roots of cubic critical equation
    poly = [1., r, s, t]
    roots = numpy.roots(poly)
    if numpy.all(numpy.imag(roots) == 0):
        # 3 real, positive roots
        roots.sort()
        dvap = numpy.real(roots[0] + 1)*_DCP
        dliq = numpy.real(roots[2] + 1)*_DCP
    else:
        # 1 real root
        ind = numpy.nonzero(numpy.imag(roots) == 0)[0][0]
        x0 = numpy.real(roots[ind])
        dflu = (x0 + 1)*_DCP
        if pr >= p1:
            dvap = dflu
        else:
            dvap = 0.
        if pr <= p2:
            dliq = dflu
        else:
            dliq = 0.
    return dliq, dvap

def _dliq_crit(temp,pres):
    """Approximate liquid water density near critical point.
    
    Approximate the density of liquid water near the critical point,
    based on a cubic equation of state.
    
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :returns: Liquid water density in kg/m3.
    """
    dliq, dvap = _dliqvap_crit(temp,pres)
    if dliq == 0:
        dliq = dvap
    return dliq

def _dvap_crit(temp,pres):
    """Approximate water vapour density near critical point.
    
    Approximate the density of water vapour near the critical point,
    based on a cubic equation of state.
    
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :returns: Water vapour density in kg/m3.
    """
    dliq, dvap = _dliqvap_crit(temp,pres)
    if dvap == 0:
        dvap = dliq
    return dvap

def _dliq_const(temp,pres):
    """Approximate liquid water density as 1e3 kg/m3.
    
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :returns: Liquid water density in kg/m3.
    """
    return 1e3

def _dliq_default(temp,pres):
    """Approximate liquid water density.
    
    Approximate the liquid water density using a mixture of approaches.
    The specific method used is (with inputs in SI units):
    
    * `_dliq_if97`: (temp<=623.15 and pres<=1e8)
        or (623.15<temp<=647.096 and pres<=1.6529e7)
    * `_dliq_const`: (temp<=623.15 and pres>1e8)
        or (623.15<temp<=650. and pres>3.5e7)
    * `_dliq_crit`: (623.15<temp<=647.096 and 1.6529e7<pres<=3.5e7)
        or (647.096<temp<=650. and pres>3.5e7)
    * `_dvap_if97`: 647.096<temp<=1073.15 and pres<=1.6529e7
    * `_dvap_ideal`: All other conditions.
    
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :returns: Liquid water density in kg/m3.
    """
    TLIM0, TLIM1, TLIM2, TLIM3 = _LIMS[0]
    PLIM0, PLIM1, PLIM2, PLIM3 = _LIMS[1]
    c1 = ( (temp <= TLIM0 and pres <= PLIM3)
        or (temp > TLIM0 and temp <= TLIM1 and pres <= PLIM0) )
    c2 = ( (temp <= TLIM0 and pres > PLIM3)
        or (temp > TLIM0 and temp <= TLIM2 and pres > PLIM2) )
    c3 = ( (temp > TLIM0 and temp <= TLIM1
            and pres > PLIM0 and pres <= PLIM2)
        or (temp > TLIM1 and temp <= TLIM2 and pres > PLIM2) )
    c4 = ( (temp > TLIM1 and temp <= TLIM3 and pres <= PLIM0) )
    if c1:
        dliq = _dliq_if97(temp,pres)
    elif c2:
        dliq = _dliq_const(temp,pres)
    elif c3:
        dliq = _dliq_crit(temp,pres)
    elif c4:
        dliq = _dvap_if97(temp,pres)
    else:
        dliq = _dvap_ideal(temp,pres)
    return  dliq

def _dvap_default(temp,pres):
    """Approximate water vapour density.
    
    Approximate the water vapour density using a mixture of approaches.
    The specific method used is (with inputs in SI units):
    
    * `_dvap_if97`: (temp<=623.15 and pres<=2.2064e7)
        or (623.15<temp<=1073.15 and pres<=1.6529e7)
    * `_dliq_if97`: temp<=623.15 and 2.2064e7<pres<=1e8
    * `_dliq_const`: (temp<=623.15 and pres>1e8)
        or (623.15<temp<=650. and pres>3.5e7)
    * `_dvap_crit`: 623.15<temp<=650. and 1.6529e7<pres<=3.5e7
    * `_dvap_ideal`: All other conditions.
    
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :returns: Water vapour density in kg/m3.
    """
    TLIM0, TLIM1, TLIM2, TLIM3 = _LIMS[0]
    PLIM0, PLIM1, PLIM2, PLIM3 = _LIMS[1]
    c1 = ( (temp <= TLIM0 and pres <= PLIM1)
        or (temp > TLIM0 and temp <= TLIM3 and pres <= PLIM0) )
    c2 = ( (temp <= TLIM0 and pres > PLIM1 and pres <= PLIM3) )
    c3 = ( (temp <= TLIM0 and pres > PLIM3)
        or (temp > TLIM0 and temp <= TLIM2 and pres > PLIM2) )
    c4 = (temp > TLIM0 and temp <= TLIM2
        and pres > PLIM0 and pres <= PLIM1)
    if c1:
        dvap = _dvap_if97(temp,pres)
    elif c2:
        dvap = _dliq_if97(temp,pres)
    elif c3:
        dvap = _dliq_const(temp,pres)
    elif c4:
        dvap = _dvap_crit(temp,pres)
    else:
        dvap = _dvap_ideal(temp,pres)
    return dvap

# Dictionaries of approximation methods
_LIQMETHODS = {'eos80': _dliq_eos80, 'if97': _dliq_if97, 'f03': _dliq_f03,
    'crit': _dliq_crit, 'const': _dliq_const, 'default': _dliq_default}
_VAPMETHODS = {'if97': _dvap_if97, 'crit': _dvap_crit, 'ideal': _dvap_ideal,
    'default': _dvap_default}


## Auxiliary equilibrium functions
def _diff_tp(d,temp,pres):
    """Calculate fluid water disequilibrium at TP.
    
    Calculate both sides of the equation
    
        given pressure = pressure of fluid water
    
    and their derivatives with respect to fluid water density. Solving
    this equation gives the equilibrium fluid water density for the
    given temperature and pressure.
    
    :arg float d: Fluid water density in kg/m3.
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :returns: Left-hand side of the equation, right-hand side,
        derivative of LHS, and derivative of RHS.
    :rtype: tuple(float)
    """
    pflu = _eq_pressure(0,0,temp,d)
    lhs = pres
    rhs = pflu
    
    pflu_d = _eq_pressure(0,1,temp,d)
    dlhs = 0.
    drhs = pflu_d
    return lhs, rhs, dlhs, drhs

def eq_tp_liq(temp,pres,dliq=None,chkvals=False,chktol=_CHKTOL,
    dliq0=None,chkbnd=False,mathargs=None):
    """Get primary liquid variables at TP.
    
    Get the value of the equilibrium liquid water density for the given
    temperature and pressure.
    
    If the calculation has already been done, the results can be passed
    to avoid unnecessary repeat calculations. If enough values are
    passed, they will be checked for consistency if chkvals is True.
    
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg dliq: Liquid water density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dliq: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg dliq0: Initial guess for the liquid water density in kg/m3. If
        None (default) then `_dliq_default` is used. A string specifier
        of the method (e.g. 'crit') can also be passed. See _LIQMETHODS
        for valid specifiers.
    :type dliq0: float or str or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Liquid water density in kg/m3.
    :raises RuntimeWarning: If a string is passed for `dliq0` that does
        not match an available method. The default is used instead.
    :raises UserWarning: If a string is passed for `dliq0` that
        specifies a function intended for water vapour.
    :raises RuntimeWarning: If the value of dliq is more consistent with
        water vapour in the subcritical region.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> eq_tp_liq(300.,1e5)
    996.556340389
    """
    if dliq is None:
        if dliq0 is None or isinstance(dliq0,float):
            dliq0fun = _dliq_default
        elif isinstance(dliq0,str):
            if dliq0 in _LIQMETHODS.keys():
                dliq0fun = _LIQMETHODS[dliq0]
            elif dliq0 in _VAPMETHODS.keys():
                warnmsg = ('Method {0} is in _VAPMETHODS and is intended for '
                    'water vapour calculations').format(dliq0)
                warnings.warn(warnmsg,UserWarning)
                dliq0fun = _VAPMETHODS[dliq0]
            else:
                warnmsg = ('Method {0} is not recognized; the default '
                    'approximation will be used instead').format(dliq0)
                warnings.warn(warnmsg,RuntimeWarning)
                dliq0fun = _LIQMETHODS['default']
        fargs = (temp,pres)
        if mathargs is None:
            mathargs = dict()
        dliq = _newton(_diff_tp,dliq0,dliq0fun,fargs=fargs,**mathargs)
    
    # Avoid accidental vapour density
    _chkflubnds(temp,dliq,chkbnd=chkbnd)
    if (temp < _TCP and pres < _PCP and dliq < _DCP):
        warnmsg = ('Vapour density {0} was found during iteration. Try again '
            'with different settings').format(dliq)
        warnings.warn(warnmsg,RuntimeWarning)
    if not chkvals:
        return dliq
    
    lhs, rhs, __, __ = _diff_tp(dliq,temp,pres)
    errs = [abs(lhs/rhs-1)]
    if max(errs) > chktol:
        warnmsg = ('Given value {0} and solution {1} disagree to more than the '
            'tolerance {2}').format(lhs,rhs,chktol)
        warnings.warn(warnmsg,RuntimeWarning)
    return dliq

def eq_tp_vap(temp,pres,dvap=None,chkvals=False,chktol=_CHKTOL,
    dvap0=None,chkbnd=False,mathargs=None):
    """Get primary vapour variables at TP.
    
    Get the value of the equilibrium water vapour density for the given
    temperature and pressure.
    
    If the calculation has already been done, the results can be passed
    to avoid unnecessary repeat calculations. If enough values are
    passed, they will be checked for consistency if chkvals is True.
    
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg dvap: Water vapour density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dvap: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg dvap0: Initial guess for the water vapour density in kg/m3. If
        None (default) then `_dvap_default` is used. A string specifier
        of the method (e.g. 'crit') can also be passed. See _VAPMETHODS
        for valid specifiers.
    :type dvap0: float or str or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Water vapour density in kg/m3.
    :raises RuntimeWarning: If a string is passed for `dvap0` that does
        not match an available method. The default is used instead.
    :raises UserWarning: If a string is passed for `dvap0` that
        specifies a function intended for liquid water.
    :raises RuntimeWarning: If the value of dvap is more consistent with
        liquid water in the subcritical region.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> eq_tp_vap(300.,1e3)
    7.22603510025e-03
    """
    if dvap is None:
        if dvap0 is None or isinstance(dvap0,float):
            dvap0fun = _dvap_default
        elif isinstance(dvap0,str):
            if dvap0 in _VAPMETHODS.keys():
                dvap0fun = _VAPMETHODS[dvap0]
            elif dvap0 in _LIQMETHODS.keys():
                warnmsg = ('Method {0} is in _LIQMETHODS and is intended for '
                    'liquid water calculations').format(dvap0)
                warnings.warn(warnmsg,UserWarning)
                dvap0fun = _LIQMETHODS[dvap0]
            else:
                warnmsg = ('Method {0} is not recognized; the default '
                    'approximation will be used instead').format(dvap0)
                warnings.warn(warnmsg,RuntimeWarning)
                dvap0fun = _VAPMETHODS['default']
        fargs = (temp,pres)
        if mathargs is None:
            mathargs = dict()
        dvap = _newton(_diff_tp,dvap0,dvap0fun,fargs=fargs,**mathargs)
    
    # Avoid accidental liquid density
    _chkflubnds(temp,dvap,chkbnd=chkbnd)
    if (temp < _TCP and pres < _PCP and dvap > _DCP):
        warnmsg = ('Liquid density {0} was found during iteration. Try again '
            'with different settings').format(dvap)
        warnings.warn(warnmsg,RuntimeWarning)
    if not chkvals:
        return dvap
    
    lhs, rhs, __, __ = _diff_tp(dvap,temp,pres)
    errs = [abs(lhs/rhs-1)]
    if max(errs) > chktol:
        warnmsg = ('Given value {0} and solution {1} disagree to more than the '
            'tolerance {2}').format(lhs,rhs,chktol)
        warnings.warn(warnmsg,RuntimeWarning)
    return dvap


## Gibbs functions
def _flu_g(drvt,drvp,temp,dflu):
    """Calculate fluid water Gibbs energy from TD.
    
    Calculate the specific Gibbs free energy of fluid water or its
    derivatives with respect to temperature and pressure from the
    temperature and density. This function is used to simplify liq_g and
    vap_g which use temperature and pressure.
    
    :arg int drvt: Number of temperature derivatives.
    :arg int drvp: Number of pressure derivatives.
    :arg float temp: Temperature in K.
    :arg float dflu: Fluid water density in kg/m3.
    :returns: Gibbs free energy in units of
        (J/kg) / K^drvt / Pa^drvp.
    """
    if (drvt,drvp) == (0,0):
        f = _flu_f(0,0,temp,dflu)
        f_d = _flu_f(0,1,temp,dflu)
        g = f + dflu*f_d
        return g
    elif (drvt,drvp) == (1,0):
        f_t = _flu_f(1,0,temp,dflu)
        g_t = f_t
        return g_t
    elif (drvt,drvp) == (0,1):
        g_p = dflu**(-1)
        return g_p
    
    # Second-order derivatives require this denominator
    p_d = _eq_pressure(0,1,temp,dflu)
    if (drvt,drvp) == (2,0):
        f_tt = _flu_f(2,0,temp,dflu)
        f_td = _flu_f(1,1,temp,dflu)
        p_t = _eq_pressure(1,0,temp,dflu)
        d_t = -p_t / p_d
        g_tt = f_tt + f_td*d_t
        return g_tt
    elif (drvt,drvp) == (1,1):
        f_td = _flu_f(1,1,temp,dflu)
        d_p = p_d**(-1)
        g_tp = f_td*d_p
        return g_tp
    elif (drvt,drvp) == (0,2):
        d_p = p_d**(-1)
        g_pp = -d_p / dflu**2
        return g_pp
    errmsg = 'Derivatives {0} not recognized'.format((drvt,drvp))
    raise ValueError(errmsg)

def liq_g(drvt,drvp,temp,pres,dliq=None,chkvals=False,chktol=_CHKTOL,
    dliq0=None,chkbnd=False,mathargs=None):
    """Calculate liquid water Gibbs free energy with derivatives.
    
    Calculate the specific Gibbs free energy of liquid water or its
    derivatives with respect to temperature and pressure.
    
    :arg int drvt: Number of temperature derivatives.
    :arg int drvp: Number of pressure derivatives.
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg dliq: Liquid water density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dliq: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg dliq0: Initial guess for the liquid water density in kg/m3. If
        None (default) then `_dliq_default` is used. A string specifier
        of the method (e.g. 'crit') can also be passed. See _LIQMETHODS
        for valid specifiers.
    :type dliq0: float or str or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Gibbs free energy in units of
        (J/kg) / K^drvt / Pa^drvp.
    :raises RuntimeWarning: If a string is passed for `dliq0` that does
        not match an available method. The default is used instead.
    :raises UserWarning: If a string is passed for `dliq0` that
        specifies a function intended for water vapour.
    :raises RuntimeWarning: If the value of dliq is more consistent with
        water vapour in the subcritical region.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> liq_g(0,0,300.,1e5)
    -5265.05045577
    >>> liq_g(0,1,300.,1e5)
    1.00345555938e-03
    >>> liq_g(2,0,300.,1e5)
    -13.9354650734
    """
    dliq = eq_tp_liq(temp,pres,dliq=dliq,chkvals=chkvals,chktol=chktol,
        dliq0=dliq0,chkbnd=chkbnd,mathargs=mathargs)
    gliq = _flu_g(drvt,drvp,temp,dliq)
    return gliq

def vap_g(drvt,drvp,temp,pres,dvap=None,chkvals=False,chktol=_CHKTOL,
    dvap0=None,chkbnd=False,mathargs=None):
    """Calculate water vapour Gibbs free energy with derivatives.
    
    Calculate the specific Gibbs free energy of water vapour or its
    derivatives with respect to temperature and pressure.
    
    :arg int drvt: Number of temperature derivatives.
    :arg int drvp: Number of pressure derivatives.
    :arg float temp: Temperature in K.
    :arg float pres: Pressure in Pa.
    :arg dvap: Water vapour density in kg/m3. If unknown, pass None
        (default) and it will be calculated.
    :type dvap: float or None
    :arg bool chkvals: If True (default False) and all values are given,
        this function will calculate the disequilibrium and raise a
        warning if the results are not within a given tolerance.
    :arg float chktol: Tolerance to use when checking values (default
        _CHKTOL).
    :arg dvap0: Initial guess for the water vapour density in kg/m3. If
        None (default) then `_dvap_default` is used. A string specifier
        of the method (e.g. 'crit') can also be passed. See _VAPMETHODS
        for valid specifiers.
    :type dvap0: float or str or None
    :arg bool chkbnd: If True then warnings are raised when the given
        values are valid but outside the recommended bounds (default
        False).
    :arg mathargs: Keyword arguments to the root-finder
        :func:`_newton <maths3.newton>` (e.g. maxiter, rtol). If None
        (default) then no arguments are passed and default parameters
        will be used.
    :returns: Gibbs free energy in units of
        (J/kg) / K^drvt / Pa^drvp.
    :raises RuntimeWarning: If a string is passed for `dvap0` that does
        not match an available method. The default is used instead.
    :raises UserWarning: If a string is passed for `dvap0` that
        specifies a function intended for liquid water.
    :raises RuntimeWarning: If the value of dvap is more consistent with
        liquid water in the subcritical region.
    :raises RuntimeWarning: If the relative disequilibrium is more than
        chktol, if chkvals is True and all values are given.
    
    :Examples:
    
    >>> vap_g(0,0,300.,1e3)
    -180090.341338
    >>> vap_g(0,1,300.,1e3)
    138.388478069
    >>> vap_g(2,0,300.,1e3)
    -6.24707163427
    """
    dvap = eq_tp_vap(temp,pres,dvap=dvap,chkvals=chkvals,chktol=chktol,
        dvap0=dvap0,chkbnd=chkbnd,mathargs=mathargs)
    gvap = _flu_g(drvt,drvp,temp,dvap)
    return gvap

