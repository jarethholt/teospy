"""Test accuracy of the level 4 mixture modules.

This module provides tests of the accuracy of the level 4 mixture modules. The prefixes of these modules list the substances in equilibrium; `liqvap4` is for pure liquid water and pure water vapour, whereas `iceair4b` is for ice and humid air.

This module can also be called from the command line as

    python testmix.py arg1 arg2 ...

The arguments can be module names, parts of module names, or substance prefixes. This will run all available tests for modules that contain the arguments in their name. For example, 'iceair4' will test the three modules `iceair4a`, `iceair4b`, `iceair4c`; 'sea' will test all modules for which seawater is a component. If no arguments are given, all available tests are run.

The functions provided by this module generate the tests for the mixture
module of the same name. Each function returns a tuple of
:class:`~tester.Tester` instances which include the functions checked,
values of the arguments, and tables of reference values. Use the `run`
method of a Tester to run the test, and `printresults` to print a
summary.

:Available modules to test:

* :mod:`liqvap4`
* :mod:`iceliq4`
* :mod:`icevap4`
* :mod:`iceair4a`

"""

__all__ = ['genliqvap4','geniceliq4','genicevap4','geniceair4a']
from tester import Tester
_DERS3 = ((0,0,0),(1,0,0),(0,1,0),(0,0,1),(2,0,0),(1,1,0),(1,0,1),
    (0,2,0),(0,1,1),(0,0,2))


## Generating Tester instances
def genliqvap4():
    """Generate liqvap4 Testers.
    """
    import liqvap4
    funs = [liqvap4.pressure,liqvap4.densityvap,liqvap4.densityliq,
        liqvap4.chempot,liqvap4.entropyliq,liqvap4.entropyvap,
        liqvap4.enthalpyliq,liqvap4.enthalpyvap,liqvap4.volumeevap,
        liqvap4.entropyevap,liqvap4.enthalpyevap]
    fargs = tuple()
    fkwargs = {'temp': 300.}
    refs = [3536.80675227,2.55896736829e-2,996.513027468,-5361.84908682,
        393.089029801,8517.38650061,112564.859854,2549854.10109,39.0772595686,
        8124.29747080,2437289.24124]
    fnames = ['pressure','densityvap','densityliq','chempot','entropyliq',
        'entropyvap','enthalpyliq','enthalpyvap','volumeevap','entropyevap',
        'enthalpyevap']
    argfmt = '({0:s}={1:3g})'
    header = 'Liquid-vapour equilibrium at temp'
    eqfun = liqvap4.eq_tp
    eqargs = fargs
    eqkwargs = fkwargs
    eqkeys = ['temp','pres','dvap','dliq']
    test_t = Tester(funs,fargs,refs,fnames,argfmt,header=header,
        fkwargs=fkwargs,eqfun=eqfun,eqargs=eqargs,eqkwargs=eqkwargs,
        eqkeys=eqkeys)
    
    funs = [liqvap4.temperature,liqvap4.densityvap,liqvap4.densityliq,
        liqvap4.chempot,liqvap4.entropyliq,liqvap4.entropyvap,
        liqvap4.enthalpyliq,liqvap4.enthalpyvap,liqvap4.volumeevap,
        liqvap4.entropyevap,liqvap4.enthalpyevap]
    fargs = tuple()
    fkwargs = {'pres': 1e4}
    refs = [318.956328924,6.81657223094e-2,989.833275365,-15259.1024273,
        649.195605196,8148.82019425,191805.944559,2583858.67179,14.6691196141,
        7499.62458905,2392052.72723]
    fnames = ['temperature','densityvap','densityliq','chempot','entropyliq',
        'entropyvap','enthalpyliq','enthalpyvap','volumeevap','entropyevap',
        'enthalpyevap']
    argfmt = '({0:s}={1:5g})'
    header = 'Liquid-vapour equilibrium at pres'
    eqfun = liqvap4.eq_tp
    eqargs = fargs
    eqkwargs = fkwargs
    eqkeys = ['temp','pres','dvap','dliq']
    test_p = Tester(funs,fargs,refs,fnames,argfmt,header=header,
        fkwargs=fkwargs,eqfun=eqfun,eqargs=eqargs,eqkwargs=eqkwargs,
        eqkeys=eqkeys)
    return (test_t, test_p)

def geniceliq4():
    """Generate iceliq4 Testers.
    """
    import iceliq4
    funs = [iceliq4.pressure,iceliq4.densityliq,iceliq4.chempot,
        iceliq4.densityice,iceliq4.enthalpyice,iceliq4.enthalpyliq,
        iceliq4.enthalpymelt,iceliq4.entropyice,iceliq4.entropyliq,
        iceliq4.entropymelt,iceliq4.volumemelt]
    fargs = tuple()
    fkwargs = {'temp': 270.}
    refs = [39313338.8825,1019.05568894,38870.0605192,921.359428514,
        -299055.938629,26110.7481094,325166.686739,-1251.57777462,
        -47.2567126291,1204.32106199,-1.04052121182e-4]
    fnames = ['pressure','densityliq','chempot','densityice','enthalpyice',
        'enthalpyliq','enthalpymelt','entropyice','entropyliq','entropymelt',
        'volumemelt']
    argfmt = '({0:s}={1:3g})'
    header = 'Ice-liquid equilibrium at temp'
    eqfun = iceliq4.eq_tp
    eqargs = fargs
    eqkwargs = fkwargs
    eqkeys = ['temp','pres','dliq']
    test_t = Tester(funs,fargs,refs,fnames,argfmt,header=header,
        fkwargs=fkwargs,eqfun=eqfun,eqargs=eqargs,eqkwargs=eqkwargs,
        eqkeys=eqkeys)
    
    funs = [iceliq4.temperature,iceliq4.densityliq,iceliq4.chempot,
        iceliq4.densityice,iceliq4.enthalpyice,iceliq4.enthalpyliq,
        iceliq4.enthalpymelt,iceliq4.entropyice,iceliq4.entropyliq,
        iceliq4.entropymelt,iceliq4.volumemelt]
    fargs = tuple()
    fkwargs = {'pres': 1e7}
    refs = [272.401648869,1004.79353660,9972.88171381,917.896690831,
        -324602.983812,6945.92700483,331548.910817,-1228.24464139,
        -11.1121012723,1217.13254011,-9.42178903288e-5]
    fnames = ['temperature','densityliq','chempot','densityice','enthalpyice',
        'enthalpyliq','enthalpymelt','entropyice','entropyliq','entropymelt',
        'volumemelt']
    argfmt = '({0:s}={1:5g})'
    header = 'Ice-liquid equilibrium at pres'
    eqfun = iceliq4.eq_tp
    eqargs = fargs
    eqkwargs = fkwargs
    eqkeys = ['temp','pres','dliq']
    test_p = Tester(funs,fargs,refs,fnames,argfmt,header=header,
        fkwargs=fkwargs,eqfun=eqfun,eqargs=eqargs,eqkwargs=eqkwargs,
        eqkeys=eqkeys)
    return (test_t,test_p)

def genicevap4():
    """Generate icevap4 Testers.
    """
    import icevap4
    funs = [icevap4.pressure,icevap4.densityvap,icevap4.chempot,
        icevap4.densityice,icevap4.enthalpyice,icevap4.enthalpyvap,
        icevap4.entropyice,icevap4.entropyvap,icevap4.volumesubl,
        icevap4.entropysubl,icevap4.enthalpysubl]
    fargs = tuple()
    fkwargs = {'temp': 270.}
    refs = [470.059067981,3.77406140772e-3,-3895.26747392,917.170465733,
        -340033.434649,2495132.21977,-1244.95617472,9255.65736018,264.965451558,
        10500.6135349,2835165.65442]
    fnames = ['pressure','densityvap','chempot','densityice','enthalpyice',
        'enthalpyvap','entropyice','entropyvap','volumesubl','entropysubl',
        'enthalpysubl']
    argfmt = '({0:s}={1:3g})'
    header = 'Ice-vapour equilibrium at temp'
    eqfun = icevap4.eq_tp
    eqargs = fargs
    eqkwargs = fkwargs
    eqkeys = ['temp','pres','dvap']
    test_t = Tester(funs,fargs,refs,fnames,argfmt,header=header,
        fkwargs=fkwargs,eqfun=eqfun,eqargs=eqargs,eqkwargs=eqkwargs,
        eqkeys=eqkeys)
    
    funs = [icevap4.temperature,icevap4.densityvap,icevap4.chempot,
        icevap4.densityice,icevap4.enthalpyice,icevap4.enthalpyvap,
        icevap4.entropyice,icevap4.entropyvap,icevap4.volumesubl,
        icevap4.entropysubl,icevap4.enthalpysubl]
    fargs = tuple()
    fkwargs = {'pres': 100.}
    refs = [252.817910215,8.57185487853e-4,-26421.2820403,919.600269745,
        -374576.247867,2463525.19629,-1377.09771247,9848.77406912,1166.60755699,
        11225.8717816,2838101.44416]
    fnames = ['temperature','densityvap','chempot','densityice','enthalpyice',
        'enthalpyvap','entropyice','entropyvap','volumesubl','entropysubl',
        'enthalpysubl']
    argfmt = '({0:s}={1:3g})'
    header = 'Ice-vapour equilibrium at pres'
    eqfun = icevap4.eq_tp
    eqargs = fargs
    eqkwargs = fkwargs
    eqkeys = ['temp','pres','dvap']
    test_p = Tester(funs,fargs,refs,fnames,argfmt,header=header,
        fkwargs=fkwargs,eqfun=eqfun,eqargs=eqargs,eqkwargs=eqkwargs,
        eqkeys=eqkeys)
    return (test_t, test_p)

def geniceair4a():
    """Generate iceair4a Testers.
    """
    import iceair4a
    funs = [iceair4a.massfractionair,iceair4a.enthalpysubl,iceair4a.densityair,
        iceair4a.densityvap,iceair4a.densityice]
    fargs = tuple()
    fkwargs = {'temp': 270., 'pres': 1e5}
    refs = [0.997058854720,2833397.471581049,1.28880078014,3.79055033080e-3,
        917.181167192]
    refs_alt = [None,2833359.27614,None,None,None]
    fnames = ['massfractionair','enthalpysubl','densityair','densityvap',
        'densityice']
    argfmt = '({0:s}={1:3g},{2:s}={3:6g})'
    header = 'Icy air at temp and pres'
    eqfun = iceair4a.eq_atpe
    eqargs = fargs
    eqkwargs = fkwargs
    eqkeys = ['airf','temp','pres','dhum']
    test_tp = Tester(funs,fargs,refs,fnames,argfmt,header=header,
        fkwargs=fkwargs,eqfun=eqfun,eqargs=eqargs,eqkwargs=eqkwargs,
        eqkeys=eqkeys,refs_alt=refs_alt)
    
    funs = [iceair4a.pressure,iceair4a.enthalpysubl,iceair4a.densityair,
        iceair4a.densityvap,iceair4a.densityice]
    fkwargs = {'airf': 0.997, 'temp': 270.}
    refs = [98034.4511233,2833421.4055958856,1.26340801028,3.79022403085e-3,
        917.180955861]
    refs_alt = [None,2833386.54980,None,None,None]
    fnames = ['pressure','enthalpysubl','densityair','densityvap','densityice']
    argfmt = '({0:s}={1:5.3f},{2:s}={3:3g})'
    header = 'Icy air at airf and temp'
    eqkwargs = fkwargs
    test_at = Tester(funs,fargs,refs,fnames,argfmt,header=header,
        fkwargs=fkwargs,eqfun=eqfun,eqargs=eqargs,eqkwargs=eqkwargs,
        eqkeys=eqkeys,refs_alt=refs_alt)
    
    funs = [iceair4a.temperature,iceair4a.enthalpysubl,iceair4a.densityair,
        iceair4a.densityvap,iceair4a.densityice]
    fkwargs = {'airf': 0.997, 'pres': 1e5}
    refs = [270.23481612550177,2833334.552629217,1.28763121402,3.86289364206e-3,
        917.147060527]
    refs_alt = [270.232024746,2833296.51317,None,None,None]
    fnames = ['temperature','enthalpysubl','densityair','densityvap',
        'densityice']
    argfmt = '({0:s}={1:5.3f},{2:s}={3:6g})'
    header = 'Icy air at airf and temp'
    eqkwargs = fkwargs
    test_ap = Tester(funs,fargs,refs,fnames,argfmt,header=header,
        fkwargs=fkwargs,eqfun=eqfun,eqargs=eqargs,eqkwargs=eqkwargs,
        eqkeys=eqkeys,refs_alt=refs_alt)
    
    funs = [iceair4a.enthalpysubl,iceair4a.densityair,iceair4a.densityvap,
        iceair4a.pressure,iceair4a.temperature,iceair4a.densityice]
    fkwargs = {'airf': 0.997, 'entr': 100.}
    refs = [2834605.615243472,0.949325026119,2.84797507836e-3,72721.4579415,
        266.514349350,917.681749114]
    refs_alt = [2834612.42351,None,None,None,None,None]
    fnames = ['enthalpysubl','densityair','densityvap','pressure','temperature',
        'densityice']
    argfmt = '({0:s}={1:5.3f},{2:s}={3:3g})'
    header = 'Icy air at airf and entr'
    eqkwargs = fkwargs
    test_ae = Tester(funs,fargs,refs,fnames,argfmt,header=header,
        fkwargs=fkwargs,eqfun=eqfun,eqargs=eqargs,eqkwargs=eqkwargs,
        eqkeys=eqkeys,refs_alt=refs_alt)
    
    funs = iceair4a.frostpoint
    fargs = (0.997,1e5)
    refs = 270.234816126
    fnames = 'frostpoint'
    argfmt = '({0:5.3f},{1:6g})'
    header = 'Frost point'
    test_fp = Tester(funs,fargs,refs,fnames,argfmt,header=header)
    
    funs = iceair4a.sublimationpressure
    fargs = (270.,1e5)
    refs = 472.041823975
    fnames = 'sublimationpressure'
    argfmt = '({0:3g},{1:6g})'
    header = 'Sublimation pressure'
    test_sp = Tester(funs,fargs,refs,fnames,argfmt,header=header)
    
    funs = iceair4a.condensationpressure
    fargs = (0.997,270.)
    refs = 98034.4511233
    fnames = 'condensationpressure'
    argfmt = '({0:5.3f},{1:3g})'
    header = 'Condensation pressure'
    test_cp = Tester(funs,fargs,refs,fnames,argfmt,header=header)
    
    funs = [iceair4a.ict,iceair4a.icl]
    fargs = (0.997,300.,1e5)
    refs = [265.224998411,64988.3931838]
    fnames = ['ict','icl']
    argfmt = '({0:5.3f},{1:3g},{2:6g})'
    header = 'ICL functions'
    eqfun = iceair4a.eq_icl
    eqargs = fargs
    eqkeys = ['dhum','ticl','picl','dhicl']
    test_icl = Tester(funs,fargs,refs,fnames,argfmt,header=header,
        eqfun=eqfun,eqargs=eqargs,eqkeys=eqkeys)
    
    funs = [iceair4a.rhfromairf_cct,iceair4a.rhfromairf_wmo]
    fargs = (0.998,270.,1e5)
    refs = [0.680395740553,0.679365943331]
    fnames = ['rhfromairf_cct','rhfromairf_wmo']
    argfmt = '({0:5.3f},{1:3g},{2:6g})'
    header = 'RH from airf'
    test_rh1 = Tester(funs,fargs,refs,fnames,argfmt,header=header)
    
    funs = [iceair4a.airffromrh_cct,iceair4a.airffromrh_wmo]
    fargs = (0.8,270.,1e5)
    refs = [0.997647924743,0.997645698908]
    fnames = ['airffromrh_cct','airffromrh_wmo']
    argfmt = '({0:3.1f},{1:3g},{2:6g})'
    header = 'airf from RH'
    test_rh2 = Tester(funs,fargs,refs,fnames,argfmt,header=header)
    return (test_tp,test_at,test_ap,test_ae,test_fp,test_sp,test_cp,test_icl,
        test_rh1,test_rh2)

def geniceair4b():
    """Generate iceair4b Testers.
    """
    CHK_IA4B_1 = {'modname': 'ice_air_4b',
        'type': 'der',
        'args': (0.5,270.,1e5),
        'geteqvals': ice_air_4b.geteqvals_wtp,
        'eqkws': ('airf','dhum'),
        'funs': ice_air_4b.ice_air_g,
        'names': 'ice_air_g',
        'ders': DERS3,
        'refs': (-2595.57166634,2382.35592988,610.264515318,0.389645501224,
            0.,-1269.41767949,0.777110408175,-7.00810930740,1.60095965101e-3,
            -3.91178603885e-6)}

    CHK_IA4B_1_ALT = {'modname': 'ice_air_4b',
        'type': 'der',
        'args': (0.5,270.,1e5),
        'geteqvals': ice_air_4b.geteqvals_wtp,
        'eqkws': ('airf','dhum'),
        'funs': ice_air_4b.ice_air_g,
        'names': 'ice_air_g',
        'ders': DERS3,
        'refs': (-2595.5716663375033,2382.35592988098,610.2645167188605,
            0.38964550105835377,0.0,-1269.4176766911055,0.777110407844233,
            -7.011965016178234,0.0016014153056509374,-3.911839886579897e-6)}

    CHK_IA4B_2 = {'modname': 'ice_air_4b',
        'type': 'fun',
        'args': (0.5,270.,1e5),
        'geteqvals': ice_air_4b.geteqvals_wtp,
        'eqkws': ('airf','dhum'),
        'funs': (ice_air_4b.ice_air_solidfraction,ice_air_4b.ice_air_vapourfraction,
            ice_air_4b.ice_air_g_cp,ice_air_4b.ice_air_g_density,
            ice_air_4b.ice_air_g_enthalpy,ice_air_4b.ice_air_g_entropy,
            ice_air_4b.ice_air_g_expansion,ice_air_4b.ice_air_g_kappa_t,
            ice_air_4b.ice_air_g_lapserate),
        'names': ('solidfraction','vapourfraction','cp','density','enthalpy',
            'entropy','expansion','kappa_t','lapserate'),
        'refs': (0.498525089434,1.47491056602e-3,1892.18951300,2.5664353800,
            -167366.990802,-610.264515318,4.10875949031e-3,1.00393460891e-5,
            2.28443875629e-4)}

    CHK_IA4B_2_ALT = {'modname': 'ice_air_4b',
        'type': 'fun',
        'args': (0.5,270.,1e5),
        'geteqvals': ice_air_4b.geteqvals_wtp,
        'eqkws': ('airf','dhum'),
        'funs': (ice_air_4b.ice_air_solidfraction,ice_air_4b.ice_air_vapourfraction,
            ice_air_4b.ice_air_g_cp,ice_air_4b.ice_air_g_density,
            ice_air_4b.ice_air_g_enthalpy,ice_air_4b.ice_air_g_entropy,
            ice_air_4b.ice_air_g_expansion,ice_air_4b.ice_air_g_kappa_t,
            ice_air_4b.ice_air_g_lapserate),
        'names': ('solidfraction','vapourfraction','cp','density','enthalpy',
            'entropy','expansion','kappa_t','lapserate'),
        'refs': (0.498525089567494,0.0014749105660236683,1893.2305543681232,
            2.566435381093336,-167366.99118042982,-610.2645167188605,
            0.004109928900246964,1.003948429009079e-5,0.00022838324589545076)}
    return None

def geniceair4c():
    """Generate iceair4c Testers.
    """
    CHK_IA4C_1 = {'modname': 'ice_air_4c',
        'type': 'der',
        'args': (0.5,1e5),
        'kwargs': {'entr': -600.},
        'geteqvals': ice_air_4c.geteqvals_wep,
        'eqkws': ('entr','airf','temp','dhum'),
        'funs': ice_air_4c.ice_air_h,
        'names': 'ice_air_h',
        'ders': DERS3,
        'refs': (-164588.106002,543.016638396,271.449994437,0.391981878510,
            224958.525864,-177.457078495,0.781782661019,0.139985868894,
            2.26912930199e-4,-3.56976697603e-6)}

    CHK_IA4C_1_ALT = {'modname': 'ice_air_4c',
        'type': 'der',
        'args': (0.5,1e5),
        'kwargs': {'entr': -600.},
        'geteqvals': ice_air_4c.geteqvals_wep,
        'eqkws': ('entr','airf','temp','dhum'),
        'funs': ice_air_4c.ice_air_h,
        'names': 'ice_air_h',
        'ders': DERS3,
        'refs': (-164588.10601652123,543.0166476377772,271.44999442990206,
            0.3919818784985212,224806.0619232994,-177.336808300561,
            0.4942223281964142,0.13989099452736128,0.00022684010853205218,
            -3.5698228708362018e-6)}

    CHK_IA4C_2 = {'modname': 'ice_air_4c',
        'type': 'fun',
        'args': (0.9,1e5),
        'kwargs': {'entr': -100.},
        'geteqvals': ice_air_4c.geteqvals_wep,
        'eqkws': ('entr','airf','temp','dhum'),
        'funs': (ice_air_4c.ice_air_h_temperature,ice_air_4c.ice_air_h_lapserate,
            ice_air_4c.ice_air_h_cp,ice_air_4c.ice_air_h_kappa_s,
            ice_air_4c.ice_air_h_density),
        'names': ('temperature','lapserate','cp','kappa_s','density'),
        'refs': (270.383680119,4.42457786755e-4,1766.52051488,8.23031581047e-6,
            1.42531895993)}

    CHK_IA4C_2_ALT = {'modname': 'ice_air_4c',
        'type': 'fun',
        'args': (0.9,1e5),
        'kwargs': {'entr': -100.},
        'geteqvals': ice_air_4c.geteqvals_wep,
        'eqkws': ('entr','airf','temp','dhum'),
        'funs': (ice_air_4c.ice_air_h_temperature,ice_air_4c.ice_air_h_lapserate,
            ice_air_4c.ice_air_h_cp,ice_air_4c.ice_air_h_kappa_s,
            ice_air_4c.ice_air_h_density),
        'names': ('temperature','lapserate','cp','kappa_s','density'),
        'refs': (270.3836801325995,0.00044209256818978005,1768.514396747884,
            8.23141751514431e-6,1.4253189598569878)}

    CHK_IA4C_3 = {'modname': 'ice_air_4c',
        'type': 'fun',
        'args': (0.9,230.,5e4,1e5),
        'geteqvals': ice_air_4c.geteqvals_pot,
        'eqkws': ('airf','dhum','apot','tpot','dpot'),
        'funs': (ice_air_4c.ice_air_pottemp,ice_air_4c.ice_air_potdensity,
            ice_air_4c.ice_air_potenthalpy),
        'names': ('pottemp','potdensity','potenthalpy'),
        'refs': (266.105208871,1.45048110422,-35781.2564451)}
    return None

def genliqair4a():
    """Generate liqair4a Testers.
    """
    CHK_LA4A_1 = {'modname': 'liq_air_4a',
        'type': 'fun',
        'args': (),
        'kwargs': {'temp': 300., 'pres': 1e5},
        'geteqvals': liq_air_4a.geteqvals_atp,
        'eqkws': ('airf','temp','pres','entr','dhum','dliq'),
        'funs': (liq_air_4a.liq_air_enthalpy_evap,liq_air_4a.liq_air_entropy_air,
            liq_air_4a.liq_air_density_air,liq_air_4a.liq_air_density_vap,
            liq_air_4a.liq_air_density_liq),
        'names': ('enthalpy_evap','entropy_air','density_air','density_vap',
            'density_liq'),
        'refs': (2434585.53919,296.711483507,1.14614215827,2.56669393257e-2,
            996.556340389)}

    CHK_LA4A_1_ALT = {'modname': 'liq_air_4a',
        'type': 'fun',
        'args': (),
        'kwargs': {'temp': 300., 'pres': 1e5},
        'geteqvals': liq_air_4a.geteqvals_atp,
        'eqkws': ('airf','temp','pres','entr','dhum','dliq'),
        'funs': (liq_air_4a.liq_air_enthalpy_evap,liq_air_4a.liq_air_entropy_air,
            liq_air_4a.liq_air_density_air,liq_air_4a.liq_air_density_vap,
            liq_air_4a.liq_air_density_liq),
        'names': ('enthalpy_evap','entropy_air','density_air','density_vap',
            'density_liq'),
        'refs': (2434606.2895444683,296.711483507,1.14614215827,2.56669393257e-2,
            996.556340389)}

    CHK_LA4A_2 = {'modname': 'liq_air_4a',
        'type': 'fun',
        'args': (),
        'kwargs': {'airf': 0.99, 'pres': 1e5},
        'geteqvals': liq_air_4a.geteqvals_atp,
        'eqkws': ('airf','temp','pres','entr','dhum','dliq'),
        'funs': (liq_air_4a.liq_air_enthalpy_evap,liq_air_4a.liq_air_entropy_air,
            liq_air_4a.liq_air_temperature,liq_air_4a.liq_air_density_air,
            liq_air_4a.liq_air_density_vap,liq_air_4a.liq_air_density_liq),
        'names': ('enthalpy_evap','entropy_air','temperature','density_air',
            'density_vap','density_liq'),
        'refs': (2465656.38630,145.863545194,287.078299795,1.20675806022,
            0.0120675806022,999.256685197)}

    CHK_LA4A_2_ALT = {'modname': 'liq_air_4a',
        'type': 'fun',
        'args': (),
        'kwargs': {'airf': 0.99, 'pres': 1e5},
        'geteqvals': liq_air_4a.geteqvals_atp,
        'eqkws': ('airf','temp','pres','entr','dhum','dliq'),
        'funs': (liq_air_4a.liq_air_enthalpy_evap,liq_air_4a.liq_air_entropy_air,
            liq_air_4a.liq_air_temperature,liq_air_4a.liq_air_density_air,
            liq_air_4a.liq_air_density_vap,liq_air_4a.liq_air_density_liq),
        'names': ('enthalpy_evap','entropy_air','temperature','density_air',
            'density_vap','density_liq'),
        'refs': (2465683.4351508557,145.863545194,287.078299795,1.20675806022,
            0.0120675806022,999.256685197)}

    CHK_LA4A_3 = {'modname': 'liq_air_4a',
        'type': 'fun',
        'args': (),
        'kwargs': {'airf': 0.99, 'temp': 300.},
        'geteqvals': liq_air_4a.geteqvals_atp,
        'eqkws': ('airf','temp','pres','entr','dhum','dliq'),
        'funs': (liq_air_4a.liq_air_enthalpy_evap,liq_air_4a.liq_air_entropy_air,
            liq_air_4a.liq_air_pressure,liq_air_4a.liq_air_density_air,
            liq_air_4a.liq_air_density_vap,liq_air_4a.liq_air_density_liq),
        'names': ('enthalpy_evap','entropy_air','pressure','density_air',
            'density_vap','density_liq'),
        'refs': (2433111.29696,-41.9991507402,223057.741750,2.57657653270,
            2.57657653270e-2,996.611581662)}

    CHK_LA4A_3_ALT = {'modname': 'liq_air_4a',
        'type': 'fun',
        'args': (),
        'kwargs': {'airf': 0.99, 'temp': 300.},
        'geteqvals': liq_air_4a.geteqvals_atp,
        'eqkws': ('airf','temp','pres','entr','dhum','dliq'),
        'funs': (liq_air_4a.liq_air_enthalpy_evap,liq_air_4a.liq_air_entropy_air,
            liq_air_4a.liq_air_pressure,liq_air_4a.liq_air_density_air,
            liq_air_4a.liq_air_density_vap,liq_air_4a.liq_air_density_liq),
        'names': ('enthalpy_evap','entropy_air','pressure','density_air',
            'density_vap','density_liq'),
        'refs': (2433303.9209508635,-41.9991507402,223057.741750,2.57657653270,
            2.57657653270e-2,996.611581662)}

    CHK_LA4A_4 = {'modname': 'liq_air_4a',
        'type': 'fun',
        'args': (),
        'kwargs': {'airf': 0.99, 'entr': 100.},
        'geteqvals': liq_air_4a.geteqvals_atp,
        'eqkws': ('airf','temp','pres','entr','dhum','dliq'),
        'funs': (liq_air_4a.liq_air_enthalpy_evap,liq_air_4a.liq_air_temperature,
            liq_air_4a.liq_air_pressure,liq_air_4a.liq_air_density_air,
            liq_air_4a.liq_air_density_vap,liq_air_4a.liq_air_density_liq),
        'names': ('enthalpy_evap','temperature','pressure','density_air',
            'density_vap','density_liq'),
        'refs': (2458121.74961,290.107386673,121546.373652,1.45154665083,
            1.45154665083e-2,998.794738784)}

    CHK_LA4A_4_ALT = {'modname': 'liq_air_4a',
        'type': 'fun',
        'args': (),
        'kwargs': {'airf': 0.99, 'entr': 100.},
        'geteqvals': liq_air_4a.geteqvals_atp,
        'eqkws': ('airf','temp','pres','entr','dhum','dliq'),
        'funs': (liq_air_4a.liq_air_enthalpy_evap,liq_air_4a.liq_air_temperature,
            liq_air_4a.liq_air_pressure,liq_air_4a.liq_air_density_air,
            liq_air_4a.liq_air_density_vap,liq_air_4a.liq_air_density_liq),
        'names': ('enthalpy_evap','temperature','pressure','density_air',
            'density_vap','density_liq'),
        'refs': (2458179.898368215,290.107386673,121546.373652,1.45154665083,
            1.45154665083e-2,998.794738784)}

    CHK_LA4A_5 = {'modname': 'liq_air_4a',
        'type': 'fun',
        'args': (0.9,300.),
        'funs': (liq_air_4a.liq_air_condensationpressure,),
        'names': ('condensationpressure',),
        'refs': (23381.2332935,)}

    CHK_LA4A_6 = {'modname': 'liq_air_4a',
        'type': 'fun',
        'args': (300.,1e5),
        'funs': (liq_air_4a.liq_air_massfraction_air,),
        'names': ('massfraction_air',),
        'refs': (0.977605797727,)}

    CHK_LA4A_7 = {'modname': 'liq_air_4a',
        'type': 'fun',
        'args': (0.99,1e5),
        'funs': (liq_air_4a.liq_air_dewpoint,),
        'names': ('dewpoint',),
        'refs': (287.078299795,)}

    CHK_LA4A_8 = {'modname': 'liq_air_4a',
        'type': 'fun',
        'args': (0.99,300.,1e5),
        'geteqvals': liq_air_4a.geteqvals_icl,
        'eqkws': ('dhum','dliq','ticl','dhicl','dlicl'),
        'funs': (liq_air_4a.liq_air_ict,liq_air_4a.liq_air_icl),
        'names': ('ict','icl'),
        'refs': (284.200207629,82723.6047631)}

    CHK_LA4A_9 = {'modname': 'liq_air_4a',
        'type': 'fun',
        'args': (0.99,300.,1e5),
        'funs': (liq_air_4a.liq_air_rhfromairf_cct,
            liq_air_4a.liq_air_rhfromairf_wmo),
        'names': ('rhfromairf_cct','rhfromairf_wmo'),
        'refs': (0.449887886959,0.440953686019)}

    CHK_LA4A_10 = {'modname': 'liq_air_4a',
        'type': 'fun',
        'args': (0.8,300.,1e5),
        'funs': (liq_air_4a.liq_air_airffromrh_cct,
            liq_air_4a.liq_air_airffromrh_wmo),
        'names': ('airffromrh_cct','airffromrh_wmo'),
        'refs': (0.982133277948,0.982004037135)}
    return None

def genliqair4b():
    """Generate liqair4b Testers.
    """
    CHK_LA4B_1 = {'modname': 'liq_air_4b',
        'type': 'der',
        'args': (0.5,300.,1e5),
        'funs': liq_air_4b.liq_air_g,
        'names': 'liq_air_g',
        'geteqvals': liq_air_4b.geteqvals_wtp,
        'eqkws': ('airf','dhum','dliq'),
        'ders': DERS3,
        'refs': (-5396.77820137,-263.455491203,-343.783393872,0.446729465555,
            0.,98.5580798842,0.891452019991,-14.0995955397,2.43183979422e-3,
            -4.62360294023e-6)}

    CHK_LA4B_1_ALT = {'modname': 'liq_air_4b',
        'type': 'der',
        'args': (0.5,300.,1e5),
        'funs': liq_air_4b.liq_air_g,
        'names': 'liq_air_g',
        'geteqvals': liq_air_4b.geteqvals_wtp,
        'eqkws': ('airf','dhum','dliq'),
        'ders': DERS3,
        'refs': (-5396.77820137,-263.455491203,-343.783393872,0.446729465555,
            0.,98.5580798842,0.891452019991,-14.226522368319394,
            0.0024533597286680994,-4.627251558752521e-6)}

    CHK_LA4B_2 = {'modname': 'liq_air_4b',
        'type': 'fun',
        'args': (0.5,300.,1e5),
        'geteqvals': liq_air_4b.geteqvals_wtp,
        'eqkws': ('airf','dhum','dliq'),
        'funs': (liq_air_4b.liq_air_g_cp,liq_air_4b.liq_air_g_density,
            liq_air_4b.liq_air_g_enthalpy,liq_air_4b.liq_air_g_entropy,
            liq_air_4b.liq_air_g_expansion,liq_air_4b.liq_air_g_kappa_t,
            liq_air_4b.liq_air_g_lapserate,liq_air_4b.liq_air_liquidfraction,
            liq_air_4b.liq_air_vapourfraction),
        'names': ('cp','density','enthalpy','entropy','expansion','kappa_t',
            'lapserate','liquidfraction','vapourfraction'),
        'refs': (4229.87866191,2.23849125053,97738.2399604,343.783393872,
            5.44365210207e-3,1.03498947276e-5,1.72475854884e-4,0.488546404734,
            1.14535952655e-2)}

    CHK_LA4B_2_ALT = {'modname': 'liq_air_4b',
        'type': 'fun',
        'args': (0.5,300.,1e5),
        'geteqvals': liq_air_4b.geteqvals_wtp,
        'eqkws': ('airf','dhum','dliq'),
        'funs': (liq_air_4b.liq_air_g_cp,liq_air_4b.liq_air_g_density,
            liq_air_4b.liq_air_g_enthalpy,liq_air_4b.liq_air_g_entropy,
            liq_air_4b.liq_air_g_expansion,liq_air_4b.liq_air_g_kappa_t,
            liq_air_4b.liq_air_g_lapserate,liq_air_4b.liq_air_liquidfraction,
            liq_air_4b.liq_air_vapourfraction),
        'names': ('cp','density','enthalpy','entropy','expansion','kappa_t',
            'lapserate','liquidfraction','vapourfraction'),
        'refs': (4267.956710495818,2.23849125053,97738.2399604,343.783393872,
            0.005491824287032767,1.0358062128281218e-5,0.00017244971505695663,
            0.488546404734,1.14535952655e-2)}
    return None

def genliqair4c():
    """Generate liqair4c Testers.
    """
    CHK_LA4C_1 = {'modname': 'liq_air_4c',
        'type': 'der',
        'args': (0.5,1e5),
        'kwargs': {'entr': 100.},
        'funs': liq_air_4c.liq_air_h,
        'names': 'liq_air_h',
        'geteqvals': liq_air_4c.geteqvals_wep,
        'eqkws': ('entr','airf','temp','dhum','dliq'),
        'ders': DERS3,
        'refs': (26898.5215492,-1681.79366113,280.393544899,0.406872930019,
            35.7689708915,1.78599925196,0.811745643965,8.91776656830e-2,
            1.55067379031e-4,-3.83770118470e-6)}

    CHK_LA4C_1_ALT = {'modname': 'liq_air_4c',
        'type': 'der',
        'args': (0.5,1e5),
        'kwargs': {'entr': 100.},
        'funs': liq_air_4c.liq_air_h,
        'names': 'liq_air_h',
        'geteqvals': liq_air_4c.geteqvals_wep,
        'eqkws': ('entr','airf','temp','dhum','dliq'),
        'ders': DERS3,
        'refs': (26898.5215492,-1681.79366113,280.393544899,0.406872930019,
            35.72888824975039,1.7839978645440038,0.814851029626437,
            0.08907773335823847,0.00015505664070260158,-3.837702338594866e-6)}

    CHK_LA4C_2 = {'modname': 'liq_air_4c',
        'type': 'fun',
        'args': (0.5,1e5),
        'kwargs': {'entr': 100.},
        'geteqvals': liq_air_4c.geteqvals_wep,
        'eqkws': ('entr','airf','temp','dhum','dliq'),
        'funs': (liq_air_4c.liq_air_h_temperature,liq_air_4c.liq_air_h_lapserate,
            liq_air_4c.liq_air_h_cp,liq_air_4c.liq_air_h_kappa_s,
            liq_air_4c.liq_air_h_density),
        'names': ('temperature','lapserate','cp','kappa_s','density'),
        'refs': (280.393544899,1.55067379031e-4,3144.21265404,9.43218607469e-6,
            2.45776980040)}

    CHK_LA4C_2_ALT = {'modname': 'liq_air_4c',
        'type': 'fun',
        'args': (0.5,1e5),
        'kwargs': {'entr': 100.},
        'geteqvals': liq_air_4c.geteqvals_wep,
        'eqkws': ('entr','airf','temp','dhum','dliq'),
        'funs': (liq_air_4c.liq_air_h_temperature,liq_air_4c.liq_air_h_lapserate,
            liq_air_4c.liq_air_h_cp,liq_air_4c.liq_air_h_kappa_s,
            liq_air_4c.liq_air_h_density),
        'names': ('temperature','lapserate','cp','kappa_s','density'),
        'refs': (280.393544899,0.00015505664070260158,3147.7400055847493,
            9.432188910707534e-6,2.45776980040)}

    CHK_LA4C_3 = {'modname': 'liq_air_4c',
        'type': 'fun',
        'args': (0.5,300.,1e4,1e5),
        'geteqvals': liq_air_4c.geteqvals_pot,
        'eqkws': ('airf','dhum','dliq','apot','tpot','dhpot','dlpot'),
        'funs': (liq_air_4c.liq_air_pottemp,liq_air_4c.liq_air_potdensity,
            liq_air_4c.liq_air_potenthalpy),
        'names': ('pottemp','potdensity','potenthalpy'),
        'refs': (348.222379217,1.22550664945,655155.797982)}
    return None

def genliqiceair4():
    """Generate liqiceair4 Testers.
    """
    CHK_LIA4_1 = {'modname': 'liq_ice_air_4',
        'type': 'fun',
        'args': (),
        'kwargs': {'airf': 0.99},
        'geteqvals': liq_ice_air_4.geteqvals_atp,
        'eqkws': ('airf','temp','pres','dhum','dliq'),
        'funs': (liq_ice_air_4.liq_ice_air_pressure,
            liq_ice_air_4.liq_ice_air_temperature),
        'names': ('pressure','temperature'),
        'refs': (38338.9622424,273.157198087)}

    CHK_LIA4_2 = {'modname': 'liq_ice_air_4',
        'type': 'fun',
        'args': (),
        'kwargs': {'temp': 273.155},
        'geteqvals': liq_ice_air_4.geteqvals_atp,
        'eqkws': ('airf','temp','pres','dhum','dliq'),
        'funs': (liq_ice_air_4.liq_ice_air_pressure,
            liq_ice_air_4.liq_ice_air_airfraction),
        'names': ('pressure','airfraction'),
        'refs': (67931.6010764,0.994366063923)}

    CHK_LIA4_3 = {'modname': 'liq_ice_air_4',
        'type': 'fun',
        'args': (),
        'kwargs': {'pres': 1e4},
        'geteqvals': liq_ice_air_4.geteqvals_atp,
        'eqkws': ('airf','temp','pres','dhum','dliq'),
        'funs': (liq_ice_air_4.liq_ice_air_airfraction,
            liq_ice_air_4.liq_ice_air_temperature),
        'names': ('airfraction','temperature'),
        'refs': (0.961024307544,273.159302793)}

    CHK_LIA4_4 = {'modname': 'liq_ice_air_4',
        'type': 'fun',
        'args': (0.1,),
        'kwargs': {'wliq': 0.2, 'wice': 0.3},
        'geteqvals': liq_ice_air_4.geteqvals_wa,
        'eqkws': ('entr','wetf','wliq','wice','airf','temp','pres','dhum','dliq'),
        'funs': (liq_ice_air_4.liq_ice_air_entropy,
            liq_ice_air_4.liq_ice_air_enthalpy,liq_ice_air_4.liq_ice_air_density,
            liq_ice_air_4.liq_ice_air_pressure,
            liq_ice_air_4.liq_ice_air_temperature),
        'names': ('entropy','enthalpy','density','pressure','temperature'),
        'refs': (3496.16306903,900361.135280,474.974398769,706.817425301,
            273.159992933)}

    CHK_LIA4_4_ALT = {'modname': 'liq_ice_air_4',
        'type': 'fun',
        'args': (0.1,),
        'kwargs': {'wliq': 0.2, 'wice': 0.3},
        'geteqvals': liq_ice_air_4.geteqvals_wa,
        'eqkws': ('entr','wetf','wliq','wice','airf','temp','pres','dhum','dliq'),
        'funs': (liq_ice_air_4.liq_ice_air_entropy,
            liq_ice_air_4.liq_ice_air_enthalpy,liq_ice_air_4.liq_ice_air_density,
            liq_ice_air_4.liq_ice_air_pressure,
            liq_ice_air_4.liq_ice_air_temperature),
        'names': ('entropy','enthalpy','density','pressure','temperature'),
        'refs': (3496.16306903,900361.135280,0.012136403756794166,706.817425301,
            273.159992933)}

    CHK_LIA4_5 = {'modname': 'liq_ice_air_4',
        'type': 'fun',
        'args': (0.99,),
        'kwargs': {'entr': 0.0, 'wetf': 0.5},
        'geteqvals': liq_ice_air_4.geteqvals_wa,
        'eqkws': ('entr','wetf','wliq','wice','airf','temp','pres','dhum','dliq'),
        'funs': (liq_ice_air_4.liq_ice_air_enthalpy,
            liq_ice_air_4.liq_ice_air_density,liq_ice_air_4.liq_ice_air_pressure,
            liq_ice_air_4.liq_ice_air_airfraction,
            liq_ice_air_4.liq_ice_air_solidfraction,
            liq_ice_air_4.liq_ice_air_liquidfraction,
            liq_ice_air_4.liq_ice_air_vapourfraction,
            liq_ice_air_4.liq_ice_air_temperature),
        'names': ('enthalpy','density','pressure','airfraction','solidfraction',
            'liquidfraction','vapourfraction','temperature'),
        'refs': (7356.12943724,7.74757979404,112016.075795,0.996583352944,
            3.30296152581e-3,3.30296152581e-3,3.39407694837e-3,
            273.151724970)}

    CHK_LIA4_5_ALT = {'modname': 'liq_ice_air_4',
        'type': 'fun',
        'args': (0.99,),
        'kwargs': {'entr': 0.0, 'wetf': 0.5},
        'geteqvals': liq_ice_air_4.geteqvals_wa,
        'eqkws': ('entr','wetf','wliq','wice','airf','temp','pres','dhum','dliq'),
        'funs': (liq_ice_air_4.liq_ice_air_enthalpy,
            liq_ice_air_4.liq_ice_air_density,liq_ice_air_4.liq_ice_air_pressure,
            liq_ice_air_4.liq_ice_air_airfraction,
            liq_ice_air_4.liq_ice_air_solidfraction,
            liq_ice_air_4.liq_ice_air_liquidfraction,
            liq_ice_air_4.liq_ice_air_vapourfraction,
            liq_ice_air_4.liq_ice_air_temperature),
        'names': ('enthalpy','density','pressure','airfraction','solidfraction',
            'liquidfraction','vapourfraction','temperature'),
        'refs': (7356.12943724,1.436115286795335,112016.075795,0.996583352944,
            3.30296152581e-3,3.30296152581e-3,3.39407694837e-3,
            273.151724970)}

    CHK_LIA4_6 = {'modname': 'liq_ice_air_4',
        'type': 'fun',
        'args': (.99,100.),
        'funs': (liq_ice_air_4.liq_ice_air_ifl,liq_ice_air_4.liq_ice_air_iml),
        'names': ('ifl','iml'),
        'refs': (83234.7314360,81605.5557728)}
    return None

def genseaair4():
    """Generate seaair4 Testers.
    """
    CHK_SA4_1 = {'modname': 'sea_air_4',
        'type': 'fun',
        'args': (),
        'kwargs': {'salt': 0.035, 'temp': 300., 'pres': 1e5},
        'geteqvals': sea_air_4.geteqvals_satp,
        'eqkws': ('salt','airf','temp','pres','dhum','dliq'),
        'funs': (sea_air_4.sea_air_massfraction_air,
            sea_air_4.sea_air_vapourpressure,sea_air_4.sea_air_entropy_air,
            sea_air_4.sea_air_enthalpy_evap,sea_air_4.sea_air_density_air,
            sea_air_4.sea_air_density_vap),
        'names': ('massfraction_air','vapourpressure','entropy_air',
            'enthalpy_evap','density_air','density_vap'),
        'refs': (0.978029483888,3485.92986681,293.150672957,2434549.89770,
            1.14642944448,2.51876465812e-2)}

    CHK_SA4_1_ALT = {'modname': 'sea_air_4',
        'type': 'fun',
        'args': (),
        'kwargs': {'salt': 0.035, 'temp': 300., 'pres': 1e5},
        'geteqvals': sea_air_4.geteqvals_satp,
        'eqkws': ('salt','airf','temp','pres','dhum','dliq'),
        'funs': (sea_air_4.sea_air_massfraction_air,
            sea_air_4.sea_air_vapourpressure,sea_air_4.sea_air_entropy_air,
            sea_air_4.sea_air_enthalpy_evap,sea_air_4.sea_air_density_air,
            sea_air_4.sea_air_density_vap),
        'names': ('massfraction_air','vapourpressure','entropy_air',
            'enthalpy_evap','density_air','density_vap'),
        'refs': (0.978029483888,3485.92986681,293.150672957,2434570.552291743,
            1.14642944448,2.51876465812e-2)}

    CHK_SA4_2 = {'modname': 'sea_air_4',
        'type': 'fun',
        'args': (),
        'kwargs': {'salt': 0.035, 'airf': 0.99, 'pres': 1e5},
        'geteqvals': sea_air_4.geteqvals_satp,
        'eqkws': ('salt','airf','temp','pres','dhum','dliq'),
        'funs': (sea_air_4.sea_air_enthalpy_evap,
            sea_air_4.sea_air_density_air,sea_air_4.sea_air_density_vap,
            sea_air_4.sea_air_condense_temp),
        'names': ('enthalpy_evap','density_air','density_vap',
            'condense_temp'),
        'refs': (2464738.79136,1.20553988598,1.20553988598e-2,287.367456468)}

    CHK_SA4_2_ALT = {'modname': 'sea_air_4',
        'type': 'fun',
        'args': (),
        'kwargs': {'salt': 0.035, 'airf': 0.99, 'pres': 1e5},
        'geteqvals': sea_air_4.geteqvals_satp,
        'eqkws': ('salt','airf','temp','pres','dhum','dliq'),
        'funs': (sea_air_4.sea_air_enthalpy_evap,
            sea_air_4.sea_air_density_air,sea_air_4.sea_air_density_vap,
            sea_air_4.sea_air_condense_temp),
        'names': ('enthalpy_evap','density_air','density_vap',
            'condense_temp'),
        'refs': (2464765.637211577,1.20553988598,1.20553988598e-2,287.367456468)}

    CHK_SA4_3 = {'modname': 'sea_air_4',
        'type': 'fun',
        'args': (0.9,0.035,300.,1e5),
        'funs': (sea_air_4.sea_air_chempot_evap,),
        'names': ('chempot_evap',),
        'refs': (1.45584069071,)}
    return None

def genseaice4():
    """Generate seaice4 Testers.
    """
    CHK_SI4_1 = {'modname': 'sea_ice_4',
        'type': 'fun',
        'args': (),
        'kwargs': {'temp': 270., 'pres': 1e5},
        'geteqvals': sea_ice_4.geteqvals_stp,
        'eqkws': ('salt','temp','pres','dliq'),
        'funs': (sea_ice_4.sea_ice_brinesalinity,sea_ice_4.sea_ice_density_ice,
            sea_ice_4.sea_ice_density_sea,sea_ice_4.sea_ice_enthalpy_melt,
            sea_ice_4.sea_ice_volume_melt,sea_ice_4.sea_ice_enthalpy_ice,
            sea_ice_4.sea_ice_enthalpy_sea,sea_ice_4.sea_ice_entropy_ice,
            sea_ice_4.sea_ice_entropy_sea,sea_ice_4.sea_ice_salinity),
        'names': ('brinesalinity','density_ice','density_sea','enthalpy_melt',
            'volume_melt','enthalpy_ice','enthalpy_sea','entropy_ice','entropy_sea',
            'salinity'),
        'refs': (0.0560264150322,917.181167192,1045.16805918,328249.119579,
            -9.181869179e-5,-339929.555499,-12742.8664892,-1244.97335506,
            -53.1667911144,0.0560264150322)}

    CHK_SI4_2 = {'modname': 'sea_ice_4',
        'type': 'fun',
        'args': (),
        'kwargs': {'salt': 0.035, 'temp': 270.},
        'geteqvals': sea_ice_4.geteqvals_stp,
        'eqkws': ('salt','temp','pres','dliq'),
        'funs': (sea_ice_4.sea_ice_meltingpressure,sea_ice_4.sea_ice_density_ice,
            sea_ice_4.sea_ice_density_sea,sea_ice_4.sea_ice_enthalpy_melt,
            sea_ice_4.sea_ice_volume_melt,sea_ice_4.sea_ice_enthalpy_ice,
            sea_ice_4.sea_ice_enthalpy_sea,sea_ice_4.sea_ice_entropy_ice,
            sea_ice_4.sea_ice_entropy_sea,sea_ice_4.sea_ice_pressure),
        'names': ('meltingpressure','density_ice','density_sea','enthalpy_melt',
            'volume_melt','enthalpy_ice','enthalpy_sea','entropy_ice','entropy_sea',
            'pressure'),
        'refs': (16132047.4385,918.898527655,1035.73670169,326829.393605,
            -9.67135426848e-5,-323205.968289,2832.94910407,-1247.71314646,
            -46.7361169560,16132047.4385)}

    CHK_SI4_3 = {'modname': 'sea_ice_4',
        'type': 'fun',
        'args': (),
        'kwargs': {'salt': 0.035, 'pres': 1e5},
        'geteqvals': sea_ice_4.geteqvals_stp,
        'eqkws': ('salt','temp','pres','dliq'),
        'funs': (sea_ice_4.sea_ice_freezingtemperature,sea_ice_4.sea_ice_dtfdp,
            sea_ice_4.sea_ice_dtfds,sea_ice_4.sea_ice_density_ice,
            sea_ice_4.sea_ice_density_sea,sea_ice_4.sea_ice_enthalpy_melt,
            sea_ice_4.sea_ice_volume_melt,sea_ice_4.sea_ice_enthalpy_ice,
            sea_ice_4.sea_ice_enthalpy_sea,sea_ice_4.sea_ice_entropy_ice,
            sea_ice_4.sea_ice_entropy_sea,sea_ice_4.sea_ice_temperature),
        'names': ('freezingtemperature','dtfdp','dtfds','density_ice','density_sea',
            'enthalpy_melt','volume_melt','enthalpy_ice','enthalpy_sea',
            'entropy_ice','entropy_sea','temperature'),
        'refs': (271.240373585,7.48210942879e-8,-56.8751336296,917.000739687,
            1028.05199645,329942.976285,-9.10140854473e-5,-337351.999358,
            -7613.19337919,-1235.44872812,-27.9264598103,271.240373585)}

    CHK_SI4_4 = {'modname': 'sea_ice_4',
        'type': 'fun',
        'args': (0.035,270.,1e5),
        'kwargs': {},
        'geteqvals': sea_ice_4.geteqvals_seaice,
        'eqkws': ('salt','dliq'),
        'funs': (sea_ice_4.sea_ice_brinefraction_seaice,sea_ice_4.sea_ice_cp_seaice,
            sea_ice_4.sea_ice_density_seaice,sea_ice_4.sea_ice_enthalpy_seaice,
            sea_ice_4.sea_ice_entropy_seaice,sea_ice_4.sea_ice_expansion_seaice,
            sea_ice_4.sea_ice_kappa_t_seaice),
        'names': ('brinefraction','cp','density','enthalpy','entropy','expansion',
            'kappa_t'),
        'refs': (0.6247053283,62868.9015126,993.156434117,-135534.287503,
            -500.445444181,-1.64731328738e-2,1.56513441348e-9)}

    CHK_SI4_5 = {'modname': 'sea_ice_4',
        'type': 'der',
        'args': (0.035,270.,1e5),
        'kwargs': {},
        'geteqvals': sea_ice_4.geteqvals_seaice,
        'eqkws': ('salt','dliq'),
        'funs': sea_ice_4.sea_ice_g,
        'names': 'sea_ice_g',
        'ders': DERS3,
        'refs': (-414.017574547,96363.7730495,500.445444181,1.00689072300e-3,0.,
            -1144.02883419,-8.62856321467e-4,-232.847783380,-1.65866446694e-5,
            -1.57591932118e-12)}

    CHK_SI4_5_ALT = {'modname': 'sea_ice_4',
        'type': 'der',
        'args': (0.035,270.,1e5),
        'kwargs': {},
        'geteqvals': sea_ice_4.geteqvals_seaice,
        'eqkws': ('salt','dliq'),
        'funs': sea_ice_4.sea_ice_g,
        'names': 'sea_ice_g',
        'ders': DERS3,
        'refs': (-414.017574547,96363.7730495,500.445444181,1.00689072300e-3,0.,
            -21272.226025171047,-0.002383040378214491,-232.847783380,
            -1.65866446694e-5,-1.57591932118e-12)}
    return None

def genseaicevap4():
    """Generate seaicevap4 Testers.
    """
    CHK_SIV4_1 = {'modname': 'sea_ice_vap_4',
        'type': 'fun',
        'args': (),
        'kwargs': {'salt': 0.035},
        'geteqvals': sea_ice_vap_4.geteqvals_stp,
        'eqkws': ('salt','temp','pres','dliq','dvap'),
        'funs': (sea_ice_vap_4.sea_ice_vap_density_vap,
            sea_ice_vap_4.sea_ice_vap_temperature,
            sea_ice_vap_4.sea_ice_vap_pressure),
        'names': ('density_vap','temperature','pressure'),
        'refs': (4.17156419318e-3,271.247815057,521.950349225)}

    CHK_SIV4_2 = {'modname': 'sea_ice_vap_4',
        'type': 'fun',
        'args': (),
        'kwargs': {'temp': 270.},
        'geteqvals': sea_ice_vap_4.geteqvals_stp,
        'eqkws': ('salt','temp','pres','dliq','dvap'),
        'funs': (sea_ice_vap_4.sea_ice_vap_salinity,
            sea_ice_vap_4.sea_ice_vap_density_vap,
            sea_ice_vap_4.sea_ice_vap_pressure),
        'names': ('salinity','density_vap','pressure'),
        'refs': (5.61489288506e-2,3.77406140772e-3,470.059067981)}

    CHK_SIV4_3 = {'modname': 'sea_ice_vap_4',
        'type': 'fun',
        'args': (),
        'kwargs': {'pres': 500.},
        'geteqvals': sea_ice_vap_4.geteqvals_stp,
        'eqkws': ('salt','temp','pres','dliq','dvap'),
        'funs': (sea_ice_vap_4.sea_ice_vap_salinity,
            sea_ice_vap_4.sea_ice_vap_density_vap,
            sea_ice_vap_4.sea_ice_vap_temperature),
        'names': ('salinity','density_vap','temperature'),
        'refs': (4.38955878828e-2,4.00364833230e-3,270.734430917)}
    return None

def gensealiq4():
    """Generate sealiq4 Testers.
    """
    CHK_SL4 = {'modname': 'sea_liq_4',
        'type': 'fun',
        'args': (0.035,300.,1e5),
        'funs': (sea_liq_4.sea_liq_osmoticpressure,),
        'names': ('osmoticpressure',),
        'refs': (2594603.20968,)}
    return None

def genseavap4():
    """Generate seavap4 Testers.
    """
    CHK_SV4_1 = {'modname': 'sea_vap_4',
        'type': 'fun',
        'args': (0.035,640.),
        'funs': (sea_vap_4.sea_vap_boilingtemperature,),
        'names': ('boilingtemperature',),
        'refs': (274.042416829,)}

    CHK_SV4_2 = {'modname': 'sea_vap_4',
        'type': 'fun',
        'args': (0.035,274.),
        'funs': (sea_vap_4.sea_vap_vapourpressure,),
        'names': ('vapourpressure',),
        'refs': (638.044692615,)}

    CHK_SV4_3 = {'modname': 'sea_vap_4',
        'type': 'fun',
        'args': (274.,640.),
        'funs': (sea_vap_4.sea_vap_brinesalinity,),
        'names': ('brinesalinity',),
        'refs': (2.94396298294e-2,)}

    CHK_SV4_4 = {'modname': 'sea_vap_4',
        'type': 'fun',
        'args': (),
        'kwargs': {'salt': 0.035, 'pres': 640.},
        'geteqvals': sea_vap_4.geteqvals_stp,
        'eqkws': ('salt','temp','pres','dlsea','dvap'),
        'funs': (sea_vap_4.sea_vap_density_sea,sea_vap_4.sea_vap_density_vap,
            sea_vap_4.sea_vap_enthalpy_evap,sea_vap_4.sea_vap_enthalpy_sea,
            sea_vap_4.sea_vap_enthalpy_vap,sea_vap_4.sea_vap_entropy_sea,
            sea_vap_4.sea_vap_entropy_vap,sea_vap_4.sea_vap_volume_evap),
        'names': ('density_sea','density_vap','enthalpy_evap','enthalpy_sea',
            'enthalpy_vap','entropy_sea','entropy_vap','volume_evap'),
        'refs': (1027.87349556,5.06324890264e-3,2498295.32187,3465.11896144,
            2502546.89358,13.0616891215,9140.56256065,197.500648110)}

    CHK_SV4_4_ALT = {'modname': 'sea_vap_4',
        'type': 'fun',
        'args': (),
        'kwargs': {'salt': 0.035, 'pres': 640.},
        'geteqvals': sea_vap_4.geteqvals_stp,
        'eqkws': ('salt','temp','pres','dlsea','dvap'),
        'funs': (sea_vap_4.sea_vap_density_sea,sea_vap_4.sea_vap_density_vap,
            sea_vap_4.sea_vap_enthalpy_evap,sea_vap_4.sea_vap_enthalpy_sea,
            sea_vap_4.sea_vap_enthalpy_vap,sea_vap_4.sea_vap_entropy_sea,
            sea_vap_4.sea_vap_entropy_vap,sea_vap_4.sea_vap_volume_evap),
        'names': ('density_sea','density_vap','enthalpy_evap','enthalpy_sea',
            'enthalpy_vap','entropy_sea','entropy_vap','volume_evap'),
        'refs': (1027.87349556,5.06324890264e-3,2498295.32187,3465.122066144071,
            2502546.89358,13.061700450797833,9140.56256065,197.500648110)}

    CHK_SV4_5 = {'modname': 'sea_vap_4',
        'type': 'fun',
        'args': (),
        'kwargs': {'salt': 0.035, 'temp': 274.},
        'geteqvals': sea_vap_4.geteqvals_stp,
        'eqkws': ('salt','temp','pres','dlsea','dvap'),
        'funs': (sea_vap_4.sea_vap_density_sea,sea_vap_4.sea_vap_density_vap,
            sea_vap_4.sea_vap_enthalpy_evap,sea_vap_4.sea_vap_enthalpy_sea,
            sea_vap_4.sea_vap_enthalpy_vap,sea_vap_4.sea_vap_entropy_sea,
            sea_vap_4.sea_vap_entropy_vap,sea_vap_4.sea_vap_pressure,
            sea_vap_4.sea_vap_volume_evap),
        'names': ('density_sea','density_vap','enthalpy_evap','enthalpy_sea',
            'enthalpy_vap','entropy_sea','entropy_vap','pressure','volume_evap'),
        'refs': (1027.87626132,5.04855547811e-3,2498395.40101,3295.96629299,
            2502469.07187,12.4443983378,9141.68990452,638.044692615,198.075461154)}

    CHK_SV4_6 = {'modname': 'sea_vap_4',
        'type': 'fun',
        'args': (),
        'kwargs': {'temp': 274., 'pres': 640.},
        'geteqvals': sea_vap_4.geteqvals_stp,
        'eqkws': ('salt','temp','pres','dlsea','dvap'),
        'funs': (sea_vap_4.sea_vap_density_sea,sea_vap_4.sea_vap_density_vap,
            sea_vap_4.sea_vap_enthalpy_evap,sea_vap_4.sea_vap_enthalpy_sea,
            sea_vap_4.sea_vap_enthalpy_vap,sea_vap_4.sea_vap_entropy_sea,
            sea_vap_4.sea_vap_entropy_vap,sea_vap_4.sea_vap_salinity,
            sea_vap_4.sea_vap_volume_evap),
        'names': ('density_sea','density_vap','enthalpy_evap','enthalpy_sea',
            'enthalpy_vap','entropy_sea','entropy_vap','salinity','volume_evap'),
        'refs': (1023.42713047,5.06403699513e-3,2498551.19875,3405.93353730,
            2502466.96633,14.0256815112,9140.27087793,2.94396298294e-2,
            197.469911653)}

    CHK_SV4_7 = {'modname': 'sea_vap_4',
        'type': 'fun',
        'args': (0.035,274.,610.),
        'geteqvals': sea_vap_4.geteqvals_seavap,
        'eqkws': ('salt','dlsea','dvap'),
        'funs': (sea_vap_4.sea_vap_cp_seavap,sea_vap_4.sea_vap_density_seavap,
            sea_vap_4.sea_vap_enthalpy_seavap,sea_vap_4.sea_vap_entropy_seavap,
            sea_vap_4.sea_vap_expansion_seavap,sea_vap_4.sea_vap_kappa_t_seavap,
            sea_vap_4.sea_vap_brinefraction_seavap),
        'names': ('cp','density','enthalpy','entropy','expansion','kappa_t',
            'brinefraction'),
        'refs': (756270.431593,7.27092786882e-3,1661118.41089,6072.50817709,
            0.458863421347,1.19990585451e-2,0.336191581475)}

    CHK_SV4_7_ALT = {'modname': 'sea_vap_4',
        'type': 'fun',
        'args': (0.035,274.,610.),
        'geteqvals': sea_vap_4.geteqvals_seavap,
        'eqkws': ('salt','dlsea','dvap'),
        'funs': (sea_vap_4.sea_vap_cp_seavap,sea_vap_4.sea_vap_density_seavap,
            sea_vap_4.sea_vap_enthalpy_seavap,sea_vap_4.sea_vap_entropy_seavap,
            sea_vap_4.sea_vap_expansion_seavap,sea_vap_4.sea_vap_kappa_t_seavap,
            sea_vap_4.sea_vap_brinefraction_seavap),
        'names': ('cp','density','enthalpy','entropy','expansion','kappa_t',
            'brinefraction'),
        'refs': (756270.431593,7.27092786882e-3,1661118.41089,6072.50817709,
            0.458863421347,1.19990585451e-2,0.3361915894821463)}

    CHK_SV4_8 = {'modname': 'sea_vap_4',
        'type': 'der',
        'args': (0.035,274.,610.),
        'geteqvals': sea_vap_4.geteqvals_seavap,
        'eqkws': ('salt','dlsea','dvap'),
        'funs': sea_vap_4.sea_vap_g,
        'names': 'sea_vap_g',
        'ders': DERS3,
        'refs': (-2748.82963245,151028.257424,-6072.50817709,137.534028399,0.,
            14965.0677011,-321.591932572,-2760.11106421,63.1093348229,
            -1.65027885871)}

    CHK_SV4_8_ALT = {'modname': 'sea_vap_4',
        'type': 'der',
        'args': (0.035,274.,610.),
        'geteqvals': sea_vap_4.geteqvals_seavap,
        'eqkws': ('salt','dlsea','dvap'),
        'funs': sea_vap_4.sea_vap_g,
        'names': 'sea_vap_g',
        'ders': DERS3,
        'refs': (-2748.82963245,151028.257424,-6072.50817709,137.534028399,0.,
            88286.38618253275,-1990.1384855543138,-2760.11106421,63.1093348229,
            -1.65027885871)}

    CHK_SV4_9 = {'modname': 'sea_vap_4',
        'type': 'fun',
        'args': (0.0035,274.,640.),
        'funs': (sea_vap_4.sea_vap_brinefraction_seavap,),
        'names': ('brinefraction',),
        'refs': (0.118887364425,)}
    return None


## Dictionary relating modules to functions
_GENDICT = {'liqvap4': genliqvap4, 'iceliq4': geniceliq4, 'icevap4': genicevap4,
    'iceair4a': geniceair4a}


## See if all values fall within the given tolerances
if __name__ == "__main__":
    # Figure out which dictionaries to include
    import sys
    if len(sys.argv) == 1:
        testlist = list()
        for (modname,genfun) in _GENDICT.items():
            testlist += list(genfun())
    else:
        modlist = list()
        testlist = list()
        for arg in sys.argv[1:]:
            for (modname,genfun) in _GENDICT.items():
                if arg in modname and modname not in modlist:
                    modlist.append(modname)
                    testlist += list(genfun())
    
    # Run tests
    for test in testlist:
        test.run()
        test.printresults()

