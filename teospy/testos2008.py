#!/usr/bin/env python
"""Test library performance against OS2008 tables.

This module provides tests of the accuracy of the fluid water, ice, and
seawater thermodynamic functions. The reference values are taken from
tables in OS2008 (Feistel et al., Ocean Sciences, 2008) and were
generated from a quadruple precision (float128) implementation of this
library.

This module can be called from the command line as::

    python testos2008.py [tab1 tab2 ...]

which will run the tests associated with the tables ``tab1``, ``tab2``,
etc. If no additional arguments are given, all tests are run. The
available tables to test are 2, 3, and A1 through A8. For example, the
command::

    python testos2008.py 2 a4 a5

will test tables 2, A4, and A5.
"""

__all__ = ['gentable2','gentable3','gentablea1','gentablea2','gentablea3',
    'gentablea4','gentablea5','gentablea6','gentablea7','gentablea8']

import warnings
from tester import Tester
import constants0
_TTP = constants0.TTP
_PTPI = constants0.PTPI
_PTPE = constants0.PTPE
_TCP = constants0.TCP
_DCP = constants0.DCP
_TCELS = constants0.TCELS
_PATM = constants0.PATM
_SAL1 = constants0.SAL1
_DERS2 = ((0,0),(1,0),(0,1),(2,0),(1,1),(0,2))
_DERS3 = ((0,0,0),(1,0,0),(0,1,0),(0,0,1),(1,0,1),(0,2,0),(0,1,1),(0,0,2))


def gentable2():
    """Generate Testers for table 2.
    
    Generate Tester instances for OS2008 table 2: liquid water and water
    vapour at the triple point temperature.
    """
    import flu1
    import flu2
    import flu3b
    import liqvap4
    temp = _TTP
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore','Step sizes are smaller',
            RuntimeWarning)
        warnings.filterwarnings('ignore','Maximum number of iterations',
            RuntimeWarning)
        __, pres, dvap, dliq = liqvap4.eq_tp(temp=temp)
    
    funs = liqvap4.pressure
    fargs = tuple()
    fkwargs = {'temp': temp}
    refs = 611.6547710078944
    fnames = 'ptriple'
    argfmt = '    '
    header = 'OS2008 table 2, TP pressure'
    chktol = 1e-12
    test_pt = Tester(funs,fargs,refs,fnames,argfmt,header=header,
        fkwargs=fkwargs,chktol=chktol)
    
    funs = flu2.pressure
    fargs = [(temp,dvap), (temp,dliq)]
    refs = [611.6547710079,611.6548]
    fnames = 'pressure'
    header = 'OS2008 table 2, fluid pressure'
    chktol = 1e-7
    test_pvl = Tester(funs,fargs,refs,fnames,argfmt,header=header,
        chktol=chktol)
    
    funs = [flu3b.vap_gibbsenergy,flu3b.liq_gibbsenergy]
    fargs = (temp,pres)
    refs = [0.6117817,0.6117817]
    fnames = ['gvap','gliq']
    header = 'OS2008 table 2, Gibbs energy'
    test_gvl = Tester(funs,fargs,refs,fnames,argfmt,header=header,
        chktol=chktol)
    return (test_pt,test_pvl,test_gvl)

def gentable3():
    """Generate Testers for table 3.
    
    Generate Tester instances for OS2008 table 3: liquid, vapour, and
    ice at the triple point.
    """
    import ice1
    import ice2
    import flu3a
    import flu3b
    temp = _TTP
    pres = _PTPI
    
    funs = flu3a.liq_g
    args1 = (temp,pres)
    fargs = [(der+args1) for der in _DERS2]
    refs = [0.611781686436,9.14297012087e-11,1.00020752302525e-3,
        -15.448497277676,-6.797921518e-8,-5.09062317110e-13]
    refs_alt = [0.6117817,5e-11,None,None,None,None]
    fnames = 'liq_g'
    argfmt = '({0:1d},{1:1d},{2:6.2f},{3:8.4f})'
    header = 'OS2008 table 3, liquid Gibbs energy'
    eqfun = flu3a.eq_tp_liq
    eqargs = args1
    eqkeys = ['dliq']
    test_lg = Tester(funs,fargs,refs,fnames,argfmt,header=header,eqfun=eqfun,
        eqargs=eqargs,eqkeys=eqkeys,refs_alt=refs_alt)
    
    funs = [flu3b.liq_enthalpy,flu3b.liq_helmholtzenergy,
        flu3b.liq_internalenergy,flu3b.liq_entropy,flu3b.liq_density,
        flu3b.liq_cp,flu3b.liq_expansion,flu3b.liq_kappa_t,flu3b.liq_kappa_s]
    fargs = args1
    refs = [0.611781661461,-2.23942234602e-10,-2.51988794168e-08,
        -9.14297012087e-11,999.7925200316,4219.91151637,-6.7965110856e-5,
        5.08956696877e-10,5.08657624753e-10]
    refs_alt = [0.6117817,5e-9,2e-9,-5e-11,None,None,None,None,None]
    fnames = ['enthalpy','helmholtzenergy','internalenergy','entropy','density',
        'cp','expansion','kappa_t','kappa_s']
    argfmt = '({0:6.2f},{1:8.4f})'
    header = 'OS2008 table 3, liquid properties'
    test_lp = Tester(funs,fargs,refs,fnames,argfmt,header=header,eqfun=eqfun,
        eqargs=eqargs,eqkeys=eqkeys,refs_alt=refs_alt)
    
    funs = flu3a.vap_g
    fargs = [(der+args1) for der in _DERS2]
    refs = [0.61178170,-9155.49340929856,205.991224917108,-6.89834540986158,
        0.75819292780093,-0.336992433127456]
    fnames = 'vap_g'
    argfmt = '({0:1d},{1:1d},{2:6.2f},{3:8.4f})'
    header = 'OS2008 table 3, vapour Gibbs energy'
    eqfun = flu3a.eq_tp_vap
    eqargs = args1
    eqkeys = ['dvap']
    test_vg = Tester(funs,fargs,refs,fnames,argfmt,header=header,eqfun=eqfun,
        eqargs=eqargs,eqkeys=eqkeys)
    
    funs = [flu3b.vap_enthalpy,flu3b.vap_helmholtzenergy,
        flu3b.vap_internalenergy,flu3b.vap_entropy,flu3b.vap_density,
        flu3b.vap_cp,flu3b.vap_expansion,flu3b.vap_kappa_t,flu3b.vap_kappa_s]
    fargs = args1
    refs = [2500915.19146570,-125994.90372461,2374919.67595939,9155.49340929856,
        4.8545757247786e-3,1884.35203215779,3.68070498200120e-3,
        1.63595528529462e-3,1.2314112503315e-3]
    fnames = ['enthalpy','helmholtzenergy','internalenergy','entropy','density',
        'cp','expansion','kappa_t','kappa_s']
    argfmt = '({0:6.2f},{1:8.4f})'
    header = 'OS2008 table 3, vapour properties'
    test_vp = Tester(funs,fargs,refs,fnames,argfmt,header=header,eqfun=eqfun,
        eqargs=eqargs,eqkeys=eqkeys)
    
    funs = ice1.ice_g
    fargs = [(der+args1) for der in _DERS2]
    refs = [0.611781703,1220.69433939648,1.09085812736669e-3,-7.67602985875191,
        1.74387964700076e-7,-1.28495941571693e-13]
    fnames = 'ice_g'
    argfmt = '({0:1d},{1:1d},{2:6.2f},{3:8.4f})'
    header = 'OS2008 table 3, ice Gibbs energy'
    test_ig = Tester(funs,fargs,refs,fnames,argfmt,header=header)
    
    funs = [ice2.enthalpy,ice2.helmholtzenergy,ice2.internalenergy,ice2.entropy,
        ice2.density,ice2.cp,ice2.expansion,ice2.kappa_t,ice2.kappa_s]
    fargs = args1
    refs = [-333444.25396784,-5.5446875e-2,-333444.92119642,-1220.69433939648,
        916.709492199488,2096.78431621667,1.59863102565909e-4,
        1.17793449347882e-10,1.14161597778763e-10,]
    fnames = ['enthalpy','helmholtzenergy','internalenergy','entropy','density',
        'cp','expansion','kappa_t','kappa_s']
    argfmt = '({0:6.2f},{1:8.4f})'
    header = 'OS2008 table 3, ice properties'
    test_ip = Tester(funs,fargs,refs,fnames,argfmt,header=header)
    return (test_lg,test_lp,test_vg,test_vp,test_ig,test_ip)

def gentablea1():
    """Generate Testers for table A1.
    
    Generate Tester instances for OS2008 table A1: Ideal and residual
    components of the fluid water Helmholtz potential.
    """
    import flu1
    temp = 500.
    dliq = 838.025
    tau = _TCP/temp
    dta = dliq/_DCP
    funs = [flu1._phi0,flu1._phi0_t,flu1._phi0_d,flu1._phi0_tt,flu1._phi0_td,
        flu1._phi0_dd,flu1._phir,flu1._phir_t,flu1._phir_d,flu1._phir_tt,
        flu1._phir_td,flu1._phir_dd]
    fargs = (tau,dta)
    refs = [2.0479773347960,9.0461110617524,0.38423674711375,-1.93249185013052,
        0.,-0.147637877832556,-3.4269320568156,-5.8140343523842,
        -0.36436665036388,-2.23440736884336,-1.1217691467031,0.85606370097461]
    fnames = ['phi0','phi0_t','phi0_d','phi0_tt','phi0_td','phi0_dd','phir',
        'phir_t','phir_d','phir_tt','phir_td','phir_dd']
    argfmt = '({0:7g},{1:7g})'
    header = 'OS2008 table A1, fluid Helmholtz potential'
    test_hp = Tester(funs,fargs,refs,fnames,argfmt,header=header)
    return (test_hp,)

def gentablea2():
    """Generate Testers for table A2.
    
    Generate Tester instances for OS2008 table A2: Fluid water
    properties for several single-phase temperatures and densities.
    """
    import flu2
    temps = (300.,)*3 + (500.,)*4 + (647.,) + (900.,)*3
    denss = (996.5560,1005.308,1188.202,0.435,4.532,838.025,1084.564,358.,0.241,
        52.615,870.769)
    funs = [flu2.pressure,flu2.cv,flu2.soundspeed,flu2.entropy]
    fargs = list(zip(temps,denss))
    refs = [
        [0.992418352e5,2.00022515e7,7.00004704e8,0.999679423e5,0.999938125e6,
            1.00003858e7,7.00000405e8,2.20384756e7,1.00062559e5,2.00000690e7,
            7.00000006e8],
        [4130.18112,4067.98347,3461.35580,1508.17541,1669.91025,3221.06219,
            3074.37693,6183.15728,1758.90657,1935.10526,2664.22350],
        [1501.51914,1534.92501,2443.57992,548.314253,535.739001,1271.28441,
            2412.00877,252.145078,724.027147,698.445674,2019.33608],
        [393.062643,387.405401,132.609616,7944.88271,6825.02725,2566.90919,
            2032.37509,4320.92307,9166.53194,6590.70225,4172.23802]
    ]
    fnames = ['pressure','cv','soundspeed','entropy']
    argfmt = '({0:3g},{1: 8.3f})'
    header = 'OS2008 table A2, fluid water properties'
    test_fp = Tester(funs,fargs,refs,fnames,argfmt,header=header)
    return (test_fp,)

def gentablea3():
    """Generate Testers for table A3.
    
    Generate Tester instances for OS2008 table A3: Fluid water
    properties for several two-phase temperatures.
    """
    import flu2
    import liqvap4
    temps = (275.,450.,625.)
    dliqs, dvaps = list(), list()
    for temp in temps:
        __, __, dvap, dliq = liqvap4.eq_tp(temp=temp)
        dliqs.append(dliq)
        dvaps.append(dvap)
    densfun = lambda temp, dens: dens
    
    funs = [flu2.pressure,densfun,flu2.enthalpy,flu2.entropy]
    fargs = list(zip(temps,dliqs))
    refs = [
        [698.45117,9.32203563628e5,1.6908269318578e7],
        [999.887406120,890.34124976167,567.0903851464],
        [7.759722016e3,0.74916158501217e6,1.6862697594697e6],
        [28.3094669595,2108.65844688447,3801.9468301114]
    ]
    fnames = ['pressure','density','enthalpy','entropy']
    argfmt = '({0:3g},{1:7.4f})'
    header = 'OS2008 table A3 liquid properties'
    test_liq = Tester(funs,fargs,refs,fnames,argfmt,header=header)
    
    fargs = list(zip(temps,dvaps))
    refs = [
        [6.98451166701e2,9.322035636282e5,1.6908269318578e7],
        [5.506649185041e-3,4.812003601257,118.2902804512],
        [2.5042899500405e6,2.77441077988962e6,2.5507162456235e6],
        [9.1066012052322e3,6.6092122132788e3,5.1850612079574e3]
    ]
    argfmt = '({0:3g},{1:9.3e})'
    header = 'OS2008 table A3 vapour properties'
    test_vap = Tester(funs,fargs,refs,fnames,argfmt,header=header)
    return (test_liq,test_vap)

def gentablea4():
    """Generate Testers for table A4.
    
    Generate Tester instances for OS2008 table A4: Ice properties at the
    triple point, at atmospheric pressure, and at low temperature.
    """
    import ice1
    import ice2
    temps = (_TTP,273.152519,100.)
    press = (_PTPE,_PATM,1e8)
    args1 = list(zip(temps,press))
    
    funs = ice1.ice_g
    # This order places change in derivatives before change in variables
    fargs = [(der+args) for args in args1 for der in _DERS2]
    refs = [0.611784135,1220.69433940,1.09085812737e-3,-7.67602985875,
        1.74387964700e-7,-1.28495941571e-13,101.34274069,1220.76932550,
        1.09084388214e-3,-7.67598233365,1.74362219972e-7,-1.28485364928e-13,
        -222296.513088,2611.95122589,1.06193389260e-3,-8.66333195517,
        2.74505162488e-8,-9.41807981761e-14]
    fnames = 'ice_g'
    argfmt = '({0:1d},{1:1d},{2:8.4f},{3:8g})'
    header = 'OS2008 table A4 ice Gibbs energy'
    test_der = Tester(funs,fargs,refs,fnames,argfmt,header=header)
    
    funs = [ice2.enthalpy,ice2.helmholtzenergy,ice2.internalenergy,ice2.entropy,
        ice2.cp,ice2.density,ice2.expansion,ice2.pcoefficient,ice2.kappa_t,
        ice2.kappa_s]
    fargs = args1
    refs = [
        [-333444.253966,-333354.873637,-483491.635676],
        [-0.055446875,-9.18701567,-328489.902347],
        [-333444.921197,-333465.403393,-589685.024936],
        [-1220.69433940,-1220.76932550,-2611.95122589],
        [2096.78431622,2096.71391024,866.333195517],
        [916.709492200,916.721463419,941.678203297],
        [1.59863102566e-4,1.59841589458e-4,2.58495528207e-5],
        [1.35714764659e6,1.35705899321e6,2.91466166994e5],
        [1.17793449348e-10,1.17785291765e-10,8.86880048115e-11],
        [1.14161597779e-10,1.14154442556e-10,8.86060982687e-11]
    ]
    fnames = ['enthalpy','helmholtzenergy','internalenergy','entropy','cp',
        'density','expansion','pcoefficient','kappa_t','kappa_s']
    argfmt = '({0:8.4f},{1:8g})'
    header = 'OS2008 table A4 ice properties'
    test_prop = Tester(funs,fargs,refs,fnames,argfmt,header=header)
    return (test_der,test_prop)

def gentablea5():
    """Generate Testers for table A5.
    
    Generate Tester instances for OS2008 table A5: Seawater component
    properties under standard conditions.
    """
    import sal2
    import flu3a
    import flu3b
    import sea3a
    salt = _SAL1
    temp = _TCELS
    pres = _PATM
    args1 = (salt,temp,pres)
    sal_g = sal2.sal_g
    
    funs = flu3a.liq_g
    fargs = [(der + args1[1:]) for der in _DERS2]
    refs = [101.342742,0.147643376,1.00015693912169e-3,-15.4473542320,
        -6.777003179e-8,-5.08928894643e-13]
    fnames = 'liq_g'
    argfmt = '({0:1d},{1:1d},{2:6.2f},{3:6g})'
    header = 'OS2008 table A5 liquid Gibbs energy'
    eqfun = flu3a.eq_tp_liq
    eqargs = args1[1:]
    eqkeys = ['dliq']
    test_lg = Tester(funs,fargs,refs,fnames,argfmt,header=header,
        eqfun=eqfun,eqargs=eqargs,eqkeys=eqkeys)
    
    funs = [flu3b.liq_enthalpy,flu3b.liq_helmholtzenergy,
        flu3b.liq_internalenergy,flu3b.liq_entropy,flu3b.liq_density,
        flu3b.liq_cp,flu3b.liq_soundspeed,flu3b.liq_gibbsenergy]
    fargs = args1[1:]
    refs = [61.0139535,1.83989364064e-3,-40.3269484,-0.147643376,
        999.84308550433,4219.44480846,1402.38253109,101.342742]
    refs_alt = [None,1.83987540918528e-3,None,None,None,None,None,None]
    fnames = ['enthalpy','helmholtzenergy','internalenergy','entropy','density',
        'cp','soundspeed','liqpotential']
    argfmt = '({0:6.2f},{1:6g})'
    header = 'OS2008 table A5 liquid properties'
    test_lp = Tester(funs,fargs,refs,fnames,argfmt,header=header,
        eqfun=eqfun,eqargs=eqargs,eqkeys=eqkeys,refs_alt=refs_alt)
    
    funs = sal_g
    fargs = [(der + args1) for der in _DERS3]
    refs = [-101.34274172939,63997.4067312299,-1.47643376346e-1,
        -2.74957224268433e-5,-7.59615411515309e-4,0.85286115117592,
        1.1928678741396e-7,5.8153517233288e-14]
    fnames = 'sal_g'
    argfmt = '({0:1d},{1:1d},{2:1d},{3:8.6f},{4:6.2f},{5:6g})'
    header = 'OS2008 table A5 salt Gibbs energy'
    test_sg = Tester(funs,fargs,refs,fnames,argfmt,header=header)
    
    funs = [
        lambda s,t,p: sal_g(0,0,0,s,t,p) - t*sal_g(0,1,0,s,t,p),
        lambda s,t,p: sal_g(0,0,0,s,t,p) - p*sal_g(0,0,1,s,t,p),
        lambda s,t,p: (sal_g(0,0,0,s,t,p) - t*sal_g(0,1,0,s,t,p)
            - p*sal_g(0,0,1,s,t,p)),
        lambda s,t,p: -sal_g(0,1,0,s,t,p),
        lambda s,t,p: -t*sal_g(0,2,0,s,t,p),
        sal2.liqpot
    ]
    fargs = args1
    refs = [-61.0139534804,-98.556737654491,-58.2279494055,0.147643376346,
        -232.95902344370,-2351.81410932936]
    fnames = ['enthalpy','helmholtzenergy','internalenergy','entropy','cp',
        'liqpot']
    argfmt = '({0:8.6f},{1:6.2f},{2:6g})'
    header = 'OS2008 table A5 salt properties'
    test_sp = Tester(funs,fargs,refs,fnames,argfmt,header=header)
    
    funs = sea3a.sea_g
    fargs = [(der + args1) for der in _DERS3]
    refs = [-1.990440523514e-8,63997.4067312299,3.748848453888e-11,
        9.7266121669485e-4,-7.59615411515309e-4,-14.59449308080,5.1516755627e-8,
        -4.507753774102e-13]
    fnames = 'sea_g'
    argfmt = '({0:1d},{1:1d},{2:1d},{3:8.6f},{4:6.2f},{5:6g})'
    header = 'OS2008 table A5 seawater Gibbs energy'
    test_bg = Tester(funs,fargs,refs,fnames,argfmt,header=header,
        eqfun=eqfun,eqargs=eqargs,eqkeys=eqkeys)
    
    funs = [sea3a.enthalpy,sea3a.helmholtzenergy,sea3a.internalenergy,
        sea3a.entropy,sea3a.density,sea3a.cp,sea3a.soundspeed,sea3a.liqpot]
    fargs = args1
    refs = [-3.0144384786934e-8,-98.55489778,-98.5548978,3.748848453888e-11,
        1028.10719995401,3986.48578502,1449.00246362,-2250.471368]
    fnames = ['enthalpy','helmholtzenergy','internalenergy','entropy','density',
        'cp','soundspeed','liqpotential']
    argfmt = '({0:8.6f},{1:6.2f},{2:6g})'
    header = 'OS2008 table A5 seawater properties'
    test_bp = Tester(funs,fargs,refs,fnames,argfmt,header=header,
        eqfun=eqfun,eqargs=eqargs,eqkeys=eqkeys)
    return (test_lg,test_lp,test_sg,test_sp,test_bg,test_bp)

def gentablea6():
    """Generate Testers for table A6.
    
    Generate Tester instances for OS2008 table A6: Seawater component
    properties at high salinity and temperature.
    """
    import sal2
    import flu3a
    import flu3b
    import sea3a
    salt = 0.1
    temp = 353.
    pres = _PATM
    args1 = (salt,temp,pres)
    sal_g = sal2.sal_g
    
    funs = flu3a.liq_g
    fargs = [(der + args1[1:]) for der in _DERS2]
    refs = [-44611.4968996,-1073.7599318875,1.02892955635611e-3,
        -11.888500004755,6.59051552339e-7,-4.746728193611e-13]
    fnames = 'liq_g'
    argfmt = '({0:1d},{1:1d},{2:3g},{3:6g})'
    header = 'OS2008 table A6 liquid Gibbs energy'
    eqfun = flu3a.eq_tp_liq
    eqargs = args1[1:]
    eqkeys = ['dliq']
    test_lg = Tester(funs,fargs,refs,fnames,argfmt,header=header,
        eqfun=eqfun,eqargs=eqargs,eqkeys=eqkeys)
    
    funs = [flu3b.liq_enthalpy,flu3b.liq_helmholtzenergy,
        flu3b.liq_internalenergy,flu3b.liq_entropy,flu3b.liq_density,
        flu3b.liq_cp,flu3b.liq_soundspeed,flu3b.liq_gibbsenergy]
    fargs = args1[1:]
    refs = [334425.7590567,-44715.7531869,334321.5027694,1073.7599318875,
        971.88383191308,4196.6405016784,1554.4629665347,-44611.4968996]
    fnames = ['enthalpy','helmholtzenergy','internalenergy','entropy','density',
        'cp','soundspeed','liqpotential']
    argfmt = '({0:3g},{1:6g})'
    header = 'OS2008 table A6 liquid properties'
    test_lp = Tester(funs,fargs,refs,fnames,argfmt,header=header,
        eqfun=eqfun,eqargs=eqargs,eqkeys=eqkeys)
    
    funs = sal_g
    fargs = [(der + args1) for der in _DERS3]
    refs = [15087.174003705,251957.275851413,156.230907404291,
        -5.7922728577126e-5,-3.0595780244234e-4,1.27922649315507,
        8.0306159575153e-7,2.1308615424374e-13]
    fnames = 'sal_g'
    argfmt = '({0:1d},{1:1d},{2:1d},{3:3.1f},{4:3g},{5:6g})'
    header = 'OS2008 table A6 salt Gibbs energy'
    test_sg = Tester(funs,fargs,refs,fnames,argfmt,header=header)
    
    funs = [
        lambda s,t,p: sal_g(0,0,0,s,t,p) - t*sal_g(0,1,0,s,t,p),
        lambda s,t,p: sal_g(0,0,0,s,t,p) - p*sal_g(0,0,1,s,t,p),
        lambda s,t,p: (sal_g(0,0,0,s,t,p) - t*sal_g(0,1,0,s,t,p)
            - p*sal_g(0,0,1,s,t,p)),
        lambda s,t,p: -sal_g(0,1,0,s,t,p),
        lambda s,t,p: -t*sal_g(0,2,0,s,t,p),
        sal2.liqpot
    ]
    fargs = args1
    refs = [-40062.33631001,15093.043024178,-40056.467289536,-156.230907404291,
        -451.566952083741,-10108.5535814360]
    fnames = ['enthalpy','helmholtzenergy','internalenergy','entropy','cp',
        'liqpot']
    argfmt = '({0:3.1f},{1:3g},{2:6g})'
    header = 'OS2008 table A6 salt properties'
    test_sp = Tester(funs,fargs,refs,fnames,argfmt,header=header)
    
    funs = sea3a.sea_g
    fargs = [(der + args1) for der in _DERS3]
    refs = [-29524.3228959,251957.275851413,-917.5290244832,9.7100682777898e-4,
        -3.0595780244234e-4,-10.609273511600,1.462113148091e-6,
        -2.61586665117e-13]
    fnames = 'sea_g'
    argfmt = '({0:1d},{1:1d},{2:1d},{3:3.1f},{4:3g},{5:6g})'
    header = 'OS2008 table A6 seawater Gibbs energy'
    test_bg = Tester(funs,fargs,refs,fnames,argfmt,header=header,
        eqfun=eqfun,eqargs=eqargs,eqkeys=eqkeys)
    
    funs = [sea3a.enthalpy,sea3a.helmholtzenergy,sea3a.internalenergy,
        sea3a.entropy,sea3a.density,sea3a.cp,sea3a.soundspeed,sea3a.liqpot]
    fargs = args1
    refs = [294363.422747,-29622.7101627,294265.0354799,917.5290244832,
        1029.85887574790,3745.073549595,3961.2783529,-54720.0504810]
    fnames = ['enthalpy','helmholtzenergy','internalenergy','entropy','density',
        'cp','soundspeed','liqpotential']
    argfmt = '({0:3.1f},{1:3g},{2:6g})'
    header = 'OS2008 table A6 seawater properties'
    test_bp = Tester(funs,fargs,refs,fnames,argfmt,header=header,
        eqfun=eqfun,eqargs=eqargs,eqkeys=eqkeys)
    return (test_lg,test_lp,test_sg,test_sp,test_bg,test_bp)

def gentablea7():
    """Generate Testers for table A7.
    
    Generate Tester instances for OS2008 table A7: Seawater component
    properties at high pressure.
    """
    import sal2
    import flu3a
    import flu3b
    import sea3a
    salt = _SAL1
    temp = _TCELS
    pres = 1e8
    args1 = (salt,temp,pres)
    sal_g = sal2.sal_g
    
    funs = flu3a.liq_g
    fargs = [(der + args1[1:]) for der in _DERS2]
    refs = [97730.38621954,8.5146650206,9.5668332915351e-4,-14.29698733876,
        1.99079570803e-7,-3.715308894234e-13]
    fnames = 'liq_g'
    argfmt = '({0:1d},{1:1d},{2:6.2f},{3:5g})'
    header = 'OS2008 table A7 liquid Gibbs energy'
    eqfun = flu3a.eq_tp_liq
    eqargs = args1[1:]
    eqkeys = ['dliq']
    test_lg = Tester(funs,fargs,refs,fnames,argfmt,header=header,
        eqfun=eqfun,eqargs=eqargs,eqkeys=eqkeys)
    
    funs = [flu3b.liq_enthalpy,flu3b.liq_helmholtzenergy,
        flu3b.liq_internalenergy,flu3b.liq_entropy,flu3b.liq_density,
        flu3b.liq_cp,flu3b.liq_soundspeed,flu3b.liq_gibbsenergy]
    fargs = args1[1:]
    refs = [95404.6054692,2062.05330419,-263.7274462,-8.5146650206,
        1045.27796139692,3905.222091582,1575.422398486,97730.3862195]
    fnames = ['enthalpy','helmholtzenergy','internalenergy','entropy','density',
        'cp','soundspeed','liqpotential']
    argfmt = '({0:6.2f},{1:5g})'
    header = 'OS2008 table A7 liquid properties'
    test_lp = Tester(funs,fargs,refs,fnames,argfmt,header=header,
        eqfun=eqfun,eqargs=eqargs,eqkeys=eqkeys)
    
    funs = sal_g
    fargs = [(der + args1) for der in _DERS3]
    refs = [-2600.93050730637,-5458.6158064880,7.5404568488117,
        -2.2912384179113e-5,-6.4075761854575e-4,0.488076973942251,
        4.6628441224121e-8,3.57345735845327e-14]
    fnames = 'sal_g'
    argfmt = '({0:1d},{1:1d},{2:1d},{3:8.6f},{4:3g},{5:5g})'
    header = 'OS2008 table A7 salt Gibbs energy'
    test_sg = Tester(funs,fargs,refs,fnames,argfmt,header=header)
    
    funs = [
        lambda s,t,p: sal_g(0,0,0,s,t,p) - t*sal_g(0,1,0,s,t,p),
        lambda s,t,p: sal_g(0,0,0,s,t,p) - p*sal_g(0,0,1,s,t,p),
        lambda s,t,p: (sal_g(0,0,0,s,t,p) - t*sal_g(0,1,0,s,t,p)
            - p*sal_g(0,0,1,s,t,p)),
        lambda s,t,p: -sal_g(0,1,0,s,t,p),
        lambda s,t,p: -t*sal_g(0,2,0,s,t,p),
        sal2.liqpot
    ]
    fargs = args1
    refs = [-4660.6062955593,-309.69208939506,-2369.3678776480,-7.5404568488117,
        -133.318225432326,-2408.9780641266]
    fnames = ['enthalpy','helmholtzenergy','internalenergy','entropy','cp',
        'liqpot']
    argfmt = '({0:8.6f},{1:6.2f},{2:6g})'
    header = 'OS2008 table A7 salt properties'
    test_sp = Tester(funs,fargs,refs,fnames,argfmt,header=header)
    
    funs = sea3a.sea_g
    fargs = [(der + args1) for der in _DERS3]
    refs = [95129.45571223,-5458.6158064880,16.0551218694,9.3377094497440e-4,
        -6.40757618545748e-4,-13.80891036482,2.45708012027e-7,
        -3.35796315839e-13]
    fnames = 'sea_g'
    argfmt = '({0:1d},{1:1d},{2:1d},{3:8.6f},{4:3g},{5:5g})'
    header = 'OS2008 table A7 seawater Gibbs energy'
    test_bg = Tester(funs,fargs,refs,fnames,argfmt,header=header,
        eqfun=eqfun,eqargs=eqargs,eqkeys=eqkeys)
    
    funs = [sea3a.enthalpy,sea3a.helmholtzenergy,sea3a.internalenergy,
        sea3a.entropy,sea3a.density,sea3a.cp,sea3a.soundspeed,sea3a.liqpot]
    fargs = args1
    refs = [90743.9991736,1752.36121479,-2633.0953238,-16.055121869,
        1070.9264465574,3771.90386615,1621.989976499,95321.40815541]
    fnames = ['enthalpy','helmholtzenergy','internalenergy','entropy','density',
        'cp','soundspeed','liqpotential']
    argfmt = '({0:8.6f},{1:6.2f},{2:5g})'
    header = 'OS2008 table A7 seawater properties'
    test_bp = Tester(funs,fargs,refs,fnames,argfmt,header=header,
        eqfun=eqfun,eqargs=eqargs,eqkeys=eqkeys)
    return (test_lg,test_lp,test_sg,test_sp,test_bg,test_bp)

def gentablea8():
    """Generate Testers for table A8.
    
    Generate Tester instances for OS2008 table A8: Properties of liquid
    water, ice, and seawater under standard conditions.
    """
    import ice1
    import ice2
    import flu3a
    import flu3b
    import sea3a
    salt = _SAL1
    temp = _TCELS
    pres = _PATM
    args1 = (salt,temp,pres)
    
    funs = flu3a.liq_g
    fargs = [(der+args1[1:]) for der in _DERS2]
    refs = [101.3427417,0.1476433763,1.00015693912169e-3,-15.44735423198,
        -6.777003179e-8,-5.089288946435e-13]
    fnames = 'liq_g'
    argfmt = '({0:1d},{1:1d},{2:6.2f},{3:6g})'
    header = 'OS2008 table A8 liquid Gibbs energy'
    eqfun = flu3a.eq_tp_liq
    eqargs = args1[1:]
    eqkeys = ['dliq']
    test_lg = Tester(funs,fargs,refs,fnames,argfmt,header=header,
        eqfun=eqfun,eqargs=eqargs,eqkeys=eqkeys)
    
    funs = [flu3b.liq_enthalpy,flu3b.liq_helmholtzenergy,
        flu3b.liq_internalenergy,flu3b.liq_entropy,flu3b.liq_density,
        flu3b.liq_cp,flu3b.liq_expansion,flu3b.liq_kappa_t,flu3b.liq_kappa_s,
        flu3b.liq_gibbsenergy]
    fargs = args1[1:]
    refs = [61.0139535,1.83989364064e-3,-40.3269484,-0.147643376,
        999.84308550433,4219.44480846,-6.7759397686e-5,5.08849036323e-10,
        5.08551764928e-10,101.3427417]
    fnames = ['enthalpy','helmholtzenergy','internalenergy','entropy','density',
        'cp','expansion','kappa_t','kappa_s','gibbsenergy']
    argfmt = '({0:6.2f},{1:6g})'
    header = 'OS2008 table A8 liquid properties'
    test_lp = Tester(funs,fargs,refs,fnames,argfmt,header=header,
        eqfun=eqfun,eqargs=eqargs,eqkeys=eqkeys)
    
    funs = ice1.ice_g
    fargs = [(der+args1[1:]) for der in _DERS2]
    refs = [98.267598403,1220.78866129995,1.09084344292644e-3,-7.67598511156675,
        1.74360824960850e-7,-1.28484824639762e-13]
    fnames = 'ice_g'
    argfmt = '({0:1d},{1:1d},{2:6.2f},{3:6g})'
    header = 'OS2008 table A8 ice Gibbs energy'
    test_ig = Tester(funs,fargs,refs,fnames,argfmt,header=header)
    
    funs = [ice2.enthalpy,ice2.helmholtzenergy,ice2.internalenergy,ice2.entropy,
        ice2.density,ice2.cp,ice2.expansion,ice2.kappa_t,ice2.kappa_s,
        ice2.chempot]
    fargs = args1[1:]
    refs = [-333360.155235679,-12.26211345,-333470.6849475,-1220.78866129995,
        916.721832527382,2096.69533322446,1.59840374979096e-4,
        1.17784843895722e-10,1.14154052637228e-10,98.267598403]
    fnames = ['enthalpy','helmholtzenergy','internalenergy','entropy','density',
        'cp','expansion','kappa_t','kappa_s','gibbsenergy']
    argfmt = '({0:6.2f},{1:6g})'
    header = 'OS2008 table A8 ice properties'
    test_ip = Tester(funs,fargs,refs,fnames,argfmt,header=header)
    
    funs = sea3a.sea_g
    fargs = [(der+args1) for der in _DERS3]
    refs = [4e-9,63997.4067312299,-6e-11,9.7266121669485e-4,
        -7.59615411515309e-4,-14.59449308080,5.15167556274e-8,
        -4.507753774102e-13]
    fnames = 'sea_g'
    argfmt = '({0:1d},{1:1d},{2:1d},{3:8.6f},{4:6.2f},{5:6g})'
    header = 'OS2008 table A8 seawater Gibbs energy'
    test_bg = Tester(funs,fargs,refs,fnames,argfmt,header=header,
        eqfun=eqfun,eqargs=eqargs,eqkeys=eqkeys)
    
    funs = [sea3a.enthalpy,sea3a.helmholtzenergy,sea3a.internalenergy,
        sea3a.entropy,sea3a.density,sea3a.cp,sea3a.expansion_t,sea3a.kappa_t,
        sea3a.kappa_s,sea3a.liqpot]
    fargs = args1
    refs = [2e-8,-98.55489778,-98.5548978,6e-11,1028.10719995401,3986.48578502,
        5.2964747379e-5,4.634454110774e-10,4.63258452069e-10,-2250.471368]
    fnames = ['enthalpy','helmholtzenergy','internalenergy','entropy','density',
        'cp','expansion','kappa_t','kappa_s','liqpot']
    argfmt = '({0:8.6f},{1:6.2f},{2:6g})'
    header = 'OS2008 table A8 seawater properties'
    test_bp = Tester(funs,fargs,refs,fnames,argfmt,header=header,
        eqfun=eqfun,eqargs=eqargs,eqkeys=eqkeys)
    return (test_lg,test_lp,test_ig,test_ip,test_bg,test_bp)


## Dictionary relating modules to functions
_GENDICT = {'2': gentable2, '3': gentable3, 'A1': gentablea1, 'A2': gentablea2,
    'A3': gentablea3, 'A4': gentablea4, 'A5': gentablea5, 'A6': gentablea6,
    'A7': gentablea7, 'A8': gentablea8}


## See if all values fall within the given tolerances
if __name__ == "__main__":
    # Figure out which dictionaries to include
    import sys
    testlist = list()
    if len(sys.argv) == 1:
        for (tableno,genfun) in _GENDICT.items():
            testlist += list(genfun())
    else:
        for arg in sys.argv[1:]:
            for (tableno,genfun) in _GENDICT.items():
                if arg.upper() == tableno:
                    testlist += list(genfun())
    
    # Run tests
    for test in testlist:
        test.run()
        test.printresults()

