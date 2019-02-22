"""Test accuracy of the sea modules.

This module provides tests of the accuracy of the seawater thermodynamic
functions. This module can be called from the command line as::

    python testsea.py mod1 mod2 ...

which will run all tests whose modules match one of ``mod1``, ``mod2``,
etc. For example, one can use ``sea3`` to check all of the level 3
modules (``sea3a``,``sea3b``,``sea3c``,``sea3d``). If no arguments are
given, all available tests are run.

The functions provided by this module generate the tests for the ``sea``
module of the same name. Each function returns a tuple of
:class:`~teospy.tests.tester.Tester` instances which include the
functions checked, values of the arguments, and tables of reference
values. Use the ``run`` method of a Tester to run the test, and
``printresults`` to print a summary.

:Available modules to test:

* :mod:`~teospy.sea3a`
* :mod:`~teospy.sea3b`
* :mod:`~teospy.sea3c`
* :mod:`~teospy.sea3d`
* :mod:`~teospy.sea5`

"""

__all__ = ['gensea3a','gensea3b','gensea3c','gensea3d','gensea5']

from teospy.tests.tester import Tester, _DERS3


## Generating Tester instances
def gensea3a():
    """Generate sea3a Testers.
    """
    from teospy import sea3a
    funs = sea3a.sea_g
    args1 = (0.035,300.,1e5)
    fargs = [(der+args1) for der in _DERS3]
    refs = [-5114.99619857,78593.7757371,-374.452240540,9.77858615182e-4,
        2247550.41118,789.934255688,-7.16682401265e-4,-13.3358324655,
        3.04605539768e-7,-4.10945807960e-13]
    fnames = 'sea_g'
    argfmt = '({0:1d},{1:1d},{2:1d},{3:5.3f},{4:3g},{5:5g})'
    header = 'Seawater Gibbs function'
    eqfun = sea3a._eq_tp_liq
    eqargs = args1[1:]
    eqkeys = ['dliq']
    testsea3a_1 = Tester(funs,fargs,refs,fnames,argfmt,header=header,
        eqfun=eqfun,eqargs=eqargs,eqkeys=eqkeys)
    
    funs = [sea3a.liqpot,sea3a.salpot,sea3a.contraction_t,sea3a.cp,
        sea3a.density,sea3a.enthalpy,sea3a.entropy,sea3a.expansion_t,
        sea3a.gibbsenergy,sea3a.internalenergy,sea3a.kappa_s,sea3a.kappa_t,
        sea3a.lapserate,sea3a.osmcoeff,sea3a.soundspeed]
    fargs = args1
    refs = [-7865.77834937,78593.7757371,0.732910044599,4000.74973964,
        1022.64272613,107220.675963,374.452240540,3.11502639583e-4,
        -5114.99619857,107122.890102,4.13135667732e-10,4.20250741344e-10,
        2.28411342567e-8,0.902777495349,1538.47940766]
    fnames = ['liqpot','salpot','contraction_t','cp','density','enthalpy',
        'entropy','expansion_t','gibbsenergy','internalenergy','kappa_s',
        'kappa_t','lapserate','osmcoeff','soundspeed']
    argfmt = '({0:5.3f},{1:3g},{2:5g})'
    header = 'Seawater properties'
    eqfun = sea3a._eq_tp_liq
    eqargs = args1[1:]
    eqkeys = ['dliq']
    testsea3a_2 = Tester(funs,fargs,refs,fnames,argfmt,header=header,
        eqfun=eqfun,eqargs=eqargs,eqkeys=eqkeys)
    
    funs = sea3a.temp_maxdensity
    fargs = (.01,1e5)
    refs = 274.950121503
    fnames = 'temp_maxdensity'
    argfmt = '({0:4g},{1:5g})'
    header = 'Maximum seawater density'
    testsea3a_3 = Tester(funs,fargs,refs,fnames,argfmt,header=header)
    return (testsea3a_1, testsea3a_2, testsea3a_3)

def gensea3b():
    """Generate sea3b Testers.
    """
    from teospy import sea3b
    funs = sea3b.sea_h
    args1 = (0.035,500.,1e5)
    fargs = [(der+args1) for der in _DERS3]
    refs = [145481.970750,86860.7980486,309.557955853,9.81092930969e-4,
        2393730.16716,72.5298236488,-6.84629317367e-4,7.72873234085e-2,
        2.86305358702e-8,-3.96880481108e-13]
    fnames = 'sea_h'
    argfmt = '({0:1d},{1:1d},{2:1d},{3:5.3f},{4:3g},{5:5g})'
    header = 'Seawater enthalpy function'
    eqfun = sea3b.eq_sep
    eqargs = args1
    eqkeys = ['temp','dliq']
    testsea3b_1 = Tester(funs,fargs,refs,fnames,argfmt,header=header,
        eqfun=eqfun,eqargs=eqargs,eqkeys=eqkeys)
    
    funs = [sea3b.expansion_t,sea3b.contraction_t]
    fargs = args1
    refs = [3.77581809091e-4,0.725209049049]
    fnames = ['expansion_t','contraction_t']
    argfmt = '({0:5.3f},{1:3g},{2:5g})'
    header = 'Seawater isothermal coefficients'
    eqfun = sea3b.eq_sep
    eqargs = args1
    eqkeys = ['temp','dliq']
    testsea3b_2 = Tester(funs,fargs,refs,fnames,argfmt,header=header,
        eqfun=eqfun,eqargs=eqargs,eqkeys=eqkeys)
    
    funs = [sea3b.pottemp,sea3b.potdensity,sea3b.potenthalpy]
    fargs = (0.035,1e7,1e5)
    fkwargs = {'temp': 300.}
    refs = [299.771869405,1022.71520130,106307.996083]
    fnames = ['pottemp','potdensity','potenthalpy']
    argfmt = '({0:5.3f},{1:5g},{2:5g},{3:s}={4:3g})'
    header = 'Seawater potential functions'
    eqfun = sea3b.eq_pot
    eqargs = fargs
    eqkwargs = fkwargs
    eqkeys = ['entr','temp','dliq','tpot','dlpot']
    testsea3b_3 = Tester(funs,fargs,refs,fnames,argfmt,header=header,
        fkwargs=fkwargs,eqfun=eqfun,eqargs=eqargs,eqkwargs=eqkwargs,
        eqkeys=eqkeys)
    
    funs = [sea3b.expansion_theta,sea3b.expansion_h,sea3b.contraction_theta,
        sea3b.contraction_h]
    fargs = (0.035,1e7,1e5)
    fkwargs = {'entr': 500.}
    refs = [3.84755380181e-4,9.60618615640e-8,0.717342103505,0.697779873590]
    fnames = ['expansion_theta','expansion_h','contraction_theta',
        'contraction_h']
    argfmt = '({0:5.3f},{1:5g},{2:5g},{3:s}={4:3g})'
    header = 'Seawater isentropic/isenthalpic coefficients'
    eqfun = sea3b.eq_pot
    eqargs = fargs
    eqkwargs = fkwargs
    eqkeys = ['entr','temp','dliq','tpot','dlpot']
    testsea3b_4 = Tester(funs,fargs,refs,fnames,argfmt,header=header,
        fkwargs=fkwargs,eqfun=eqfun,eqargs=eqargs,eqkwargs=eqkwargs,
        eqkeys=eqkeys)
    return (testsea3b_1, testsea3b_2, testsea3b_3, testsea3b_4)

def gensea3c():
    """Generate sea3c Testers.
    """
    from teospy import sea3c
    funs = sea3c.entropy
    fargs = (0.035,1e5)
    fkwargs = {'enth': 1e5}
    refs = 350.310622663
    fnames = 'entropy'
    argfmt = '({0:5.3f},{1:5g},{2:s}={3:5g})'
    header = 'Seawater entropy at {0:3g} Pa'.format(fargs[-1])
    eqfun = sea3c.eq_shp
    eqargs = fargs
    eqkwargs = fkwargs
    eqkeys = ['enth','temp','dliq']
    testsea3c_1 = Tester(funs,fargs,refs,fnames,argfmt,fkwargs=fkwargs,
        header=header,eqfun=eqfun,eqargs=eqargs,eqkwargs=eqkwargs,eqkeys=eqkeys)
    
    funs = [sea3c.temperature,sea3c.density,sea3c.contraction_t,
        sea3c.expansion_t]
    fargs = (0.035,1e7)
    fkwargs = {'enth': 1e5}
    refs = [295.985682129,1028.10986556,0.728755239644,2.89480851145e-4]
    fnames = ['temperature','density','contraction_t','expansion_t']
    argfmt = '({0:5.3f},{1:5g},{2:s}={3:5g})'
    header = 'Seawater functions at {0:3g} Pa'.format(fargs[-1])
    eqfun = sea3c.eq_shp
    eqargs = fargs
    eqkwargs = fkwargs
    eqkeys = ['enth','temp','dliq']
    testsea3c_2 = Tester(funs,fargs,refs,fnames,argfmt,fkwargs=fkwargs,
        header=header,eqfun=eqfun,eqargs=eqargs,eqkwargs=eqkwargs,eqkeys=eqkeys)
    
    funs = [sea3c.temperature,sea3c.pottemp,sea3c.density,sea3c.potdensity,
        sea3c.contraction_t,sea3c.contraction_theta,sea3c.contraction_h,
        sea3c.expansion_t,sea3c.expansion_theta,sea3c.expansion_h]
    fargs = (0.035,1e7,1e5)
    fkwargs = {'hpot': 1e5}
    refs = [298.413424848,298.194955182,1027.36529798,1023.20527737,
        0.726349317428,0.726099732703,0.714531922616,3.07242256461e-4,
        3.09134848554e-4,7.72876772245e-8]
    fnames = ['temperature','pottemp','density','potdensity','contraction_t',
        'contraction_theta','contracion_h','expansion_t','expansion_theta',
        'expansion_h']
    argfmt = '({0:5.3f},{1:5g},ppot={2:5g},{3:s}={4:5g})'
    header = 'Seawater functions with potential enthalpy'
    eqfun = sea3c.eq_pot
    eqargs = fargs
    eqkwargs = fkwargs
    eqkeys = ['enth','temp','dliq','hpot','tpot','dlpot']
    testsea3c_3 = Tester(funs,fargs,refs,fnames,argfmt,fkwargs=fkwargs,
        header=header,eqfun=eqfun,eqargs=eqargs,eqkwargs=eqkwargs,eqkeys=eqkeys)
    
    funs = [sea3c.contraction_t,sea3c.contraction_theta,sea3c.contraction_h,
        sea3c.expansion_t,sea3c.expansion_theta,sea3c.expansion_h]
    fargs = (0.035,1e7,1e5)
    fkwargs = {'temp': 300.}
    refs = [0.724913833446,0.724667977117,0.712069013013,3.18513471410e-4,
        3.20454167783e-4,8.01009066333e-8]
    fnames = ['contraction_t','contraction_theta','contraction_h',
        'expansion_t','expansionta','expansion_h']
    argfmt = '({0:5.3f},{1:5g},ppot={2:5g},{3:s}={4:3g})'
    header = 'Seawater potential functions at in-situ temperature'
    eqfun = sea3c.eq_pot
    eqargs = fargs
    eqkwargs = fkwargs
    eqkeys = ['enth','temp','dliq','hpot','tpot','dlpot']
    testsea3c_4 = Tester(funs,fargs,refs,fnames,argfmt,fkwargs=fkwargs,
        header=header,eqfun=eqfun,eqargs=eqargs,eqkwargs=eqkwargs,eqkeys=eqkeys)
    
    funs = [sea3c.contraction_t,sea3c.contraction_theta,sea3c.contraction_h,
        sea3c.expansion_t,sea3c.expansion_theta,sea3c.expansion_h]
    fargs = (0.035,1e7,1e5)
    fkwargs = {'tpot': 300.}
    refs = [0.724714253918,0.724468894946,0.711718411190,3.20122324740e-4,
        3.2206971083e-4,8.05023387611e-8]
    fnames = ['contraction_t','contraction_theta','contraction_h',
        'expansion_t','expansionta','expansion_h']
    argfmt = '({0:5.3f},{1:5g},ppot={2:5g},{3:s}={4:3g})'
    header = 'Seawater potential functions at potential temperature'
    eqfun = sea3c.eq_pot
    eqargs = fargs
    eqkwargs = fkwargs
    eqkeys = ['enth','temp','dliq','hpot','tpot','dlpot']
    testsea3c_5 = Tester(funs,fargs,refs,fnames,argfmt,fkwargs=fkwargs,
        header=header,eqfun=eqfun,eqargs=eqargs,eqkwargs=eqkwargs,eqkeys=eqkeys)
    return (testsea3c_1, testsea3c_2, testsea3c_3, testsea3c_4, testsea3c_5)

def gensea3d():
    """Generate sea3d Testers.
    """
    from teospy import sea3d
    funs = sea3d.salinity
    fargs = (273.15,101325.,1028.)
    refs = 3.50315257709e-2
    fnames = 'salinity'
    argfmt = '({0:6.2f},{1:6g},{2:4g})'
    header = 'Seawater salinity function'
    testsea3d = Tester(funs,fargs,refs,fnames,argfmt,header=header)
    return (testsea3d,)

def gensea5():
    """Generate sea5 Testers.
    """
    from teospy import sea5
    funs = sea5.tconfromtpot
    fargs = (0.035,300.)
    refs = 300.010069445
    fnames = 'tconfromtpot'
    argfmt = '({0:5.3f},{1:3g})'
    header = 'Conservative from potential temp'
    test_tc = Tester(funs,fargs,refs,fnames,argfmt,header=header)
    
    funs = sea5.tpotfromtcon
    fargs = (0.035,300.010069445)
    refs = 300.
    fnames = 'tpotfromtcon'
    argfmt = '({0:5.3f},{1:13.9f})'
    header = 'Potential from conservative temp'
    test_tp = Tester(funs,fargs,refs,fnames,argfmt,header=header)
    
    funs = [sea5.expansion_tcon,sea5.expansion_tpot,sea5.expansion_t,
        sea5.contraction_tcon,sea5.contraction_tpot,sea5.contraction_t]
    fargs = (0.035,300.,1e8)
    refs = [3.91772847589e-4,3.92515634064e-4,3.73608885178e-4,0.649596383654,
        0.663973579411,0.666238827368]
    fnames = ['expansion_tcon','expansion_tpot','expansion_t',
        'contraction_tcon','contraction_tpot','contraction_t']
    argfmt = '({0:5.3f},{1:3g},{2:3g})'
    header = 'Expansion/contraction coefficients'
    eqfun = sea5._eq_pot
    eqargs = (0.035,1e8,sea5._PATM)
    eqkwargs = {'temp': 300.}
    eqkeys = ['enth','temp','dliq','hpot','tpot','dlpot']
    keepkeys = ['dliq','tpot','dlpot']
    test_ec = Tester(funs,fargs,refs,fnames,argfmt,header=header,eqfun=eqfun,
        eqargs=eqargs,eqkwargs=eqkwargs,eqkeys=eqkeys,keepkeys=keepkeys)
    
    funs = [sea5.cabb_tcon,sea5.cabb_tpot,sea5.thrmb_tcon,sea5.thrmb_tpot]
    fargs = (0.035,300.,1e5)
    refs = [8.87852779114e-6,8.33874458992e-6,1.48109265005e-12,
        1.45941002500e-12]
    refs_alt = [8.61252567438267e-6,8.33874537690444e-6,1.48109271668362e-12,
        1.45941010702991e-12]
    fnames = ['cabb_tcon','cabb_tpot','thrmb_tcon','thrmb_tpot']
    argfmt = '({0:5.3f},{1:3g},{2:6g})'
    header = 'Cabbeling/thermobaric coefficients'
    eqfun = sea5._eq_pot
    eqargs = (0.035,1e5,sea5._PATM)
    eqkwargs = {'temp': 300.}
    chktol = 1e-7  # Because finite differences are inaccurate
    test_ct = Tester(funs,fargs,refs,fnames,argfmt,header=header,eqfun=eqfun,
        eqargs=eqargs,eqkwargs=eqkwargs,eqkeys=eqkeys,keepkeys=keepkeys,
        refs_alt=refs_alt,chktol=chktol)
    return (test_tc, test_tp, test_ec, test_ct)


## Dictionary relating modules to functions
_GENDICT = {'sea3a': gensea3a, 'sea3b': gensea3b, 'sea3c': gensea3c,
    'sea3d': gensea3d, 'sea5': gensea5}


## See if all values fall within the given tolerances
if __name__ == "__main__":
    # Figure out which dictionaries to include
    import sys
    testlist = list()
    if len(sys.argv) == 1:
        for (modname,genfun) in _GENDICT.items():
            testlist += list(genfun())
    else:
        for arg in sys.argv[1:]:
            for (modname,genfun) in _GENDICT.items():
                if arg in modname:
                    testlist += list(genfun())
    
    # Run tests
    for test in testlist:
        test.run()
        test.printresults()

