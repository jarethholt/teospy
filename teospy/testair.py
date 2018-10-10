"""Test accuracy of the air modules.

This module provides tests of the accuracy of the dry and humid air
thermodynamic functions. This module can be called from the command
line as

    python testair.py mod1 mod2 ...

which will run all tests whose modules match one of mod1, mod2, etc. For
example, one can use 'air3' to check all of the level 3 modules (air3a,
air3b, air3c). If no arguments are given, all available tests are run.

The functions provided by this module generate the tests for the `air`
module of the same name. Each function returns a tuple of
:class:`~tester.Tester` instances which include the functions checked,
values of the arguments, and tables of reference values. Use the `run`
method of a Tester to run the test, and `printresults` to print a
summary.

:Available modules to test:

* :mod:`convert0`
* :mod:`air1`
* :mod:`air2`
* :mod:`air3a`
* :mod:`air3b`
* :mod:`air3c`
* :mod:`air5`

"""

__all__ = ['gencnv0','genair1','genair2','genair3a','genair3b','genair3c',
    'genair5']
from tester import Tester
_DERS2 = ((0,0),(1,0),(0,1),(2,0),(1,1),(0,2))
_DERS3 = ((0,0,0),(1,0,0),(0,1,0),(0,0,1),(2,0,0),(1,1,0),(1,0,1),
    (0,2,0),(0,1,1),(0,0,2))


## Generating Tester instances
def gencnv0():
    """Generate convert0 Testers.
    """
    import convert0
    funs = [convert0.air_molarmass,convert0.air_molfractiondry,
        convert0.air_molfractionvap,convert0.air_massfractiondry,
        convert0.air_massfractionvap]
    fargs = (0.5,)
    refs = [0.222142374909e-1,0.383460809718,0.616539190282,0.616539190282,
        0.383460809718]
    fnames = ['molarmass','molfracdry','molfracvap','massfracdry','massfracvap']
    argfmt = '({0:3.1f})'
    header = 'convert0 air functions'
    testcnv0 = Tester(funs,fargs,refs,fnames,argfmt,header=header)
    return (testcnv0,)

def genair1():
    """Generate air1 Testers.
    """
    import air1
    funs = air1.dry_f
    args1 = (300.,1e-3)
    fargs = [(der+args1) for der in _DERS2]
    refs = [-696239.965190,-2124.55145456,86114714.9596,-2.39242390806,
        287049.624545,-86114738036.1]
    fnames = 'dry_f'
    argfmt = '({0:1d},{1:1d},{2:3g},{3:5g})'
    header = 'air1 dry_f derivatives'
    testair1 = Tester(funs,fargs,refs,fnames,argfmt,header=header)
    return (testair1,)

def genair2():
    """Generate air2 Testers.
    """
    import air2
    funs = air2._air_f_mix
    args1 = (0.9,300.,1.)
    fargs = [(der+args1) for der in _DERS3]
    refs = [-25.8379179182,233.827370434,0.164195952060,-26.2357498619,
        500.273928155,-1.53932744055,241.520643317,-0.687329742959e-3,
        0.172192606103,-0.795663887493]
    fnames = 'f_mix'
    argfmt = '({0:1d},{1:1d},{2:1d},{3:3.1f},{4:3g},{5:1g})'
    header = 'air2 f_mix derivatives'
    testair2_1 = Tester(funs,fargs,refs,fnames,argfmt,header=header)
    
    funs = air2.air_f
    args1 = (0.9,300.,1.)
    fargs = [(der+args1) for der in _DERS3]
    refs = [-95019.5943231,-205645.554995,-940.175394023,91175.3848662,
        1447768.46379,7443.09771950,-48847.9096826,-2.96482218054,
        312.063110700,-91421.4440689]
    fnames = 'air_f'
    argfmt = '({0:1d},{1:1d},{2:1d},{3:3.1f},{4:3g},{5:1g})'
    header = 'air2 air_f derivatives'
    testair2_2 = Tester(funs,fargs,refs,fnames,argfmt,header=header)
    
    funs = [air2.cp,air2.cv,air2.enthalpy,air2.entropy,air2.expansion,
        air2.gibbsenergy,air2.internalenergy,air2.kappa_s,air2.kappa_t,
        air2.lapserate,air2.pressure,air2.soundspeed]
    fargs = (0.9,300.,1.)
    refs = [1210.74031058,889.446654163,278208.408750,940.175394023,
        3.43193033077e-03,-3844.20945693,187033.023884,8.07913626816e-06,
        1.09975521396e-05,8.50371537341e-04,91175.3848662,351.817577078]
    fnames = ['cp','cv','enthalpy','entropy','expansion','gibbsenergy',
        'internalenergy','kappa_s','kappa_t','lapserate','pressure',
        'soundspeed']
    argfmt = '({0:3.1f},{1:3g},{2:1g})'
    header = 'air2 functions'
    testair2_3 = Tester(funs,fargs,refs,fnames,argfmt,header=header)
    return (testair2_1, testair2_2, testair2_3)

def genair3a():
    """Generate air3a Testers.
    """
    import air3a
    funs = air3a.eq_atp
    fargs = (0.9,300.,1e5)
    refs = 1.09708772444
    fnames = 'eq_atp'
    argfmt = '({0:3.1f},{1:3g},{2:5g})'
    header = 'air3a equilibrium function'
    testair3a_1 = Tester(funs,fargs,refs,fnames,argfmt,header=header)
    
    funs = air3a.air_g
    args1 = (0.9,300.,1e5)
    fargs = [(der+args1) for der in _DERS3]
    eqfun = air3a.eq_atp
    eqkeys = ['dhum']
    eqargs = args1
    refs = [4577.93065689,-210141.953243,-911.170080461,0.911504137473,
        1415779.23410,7566.34779196,-0.483353002176,-4.15449972148,
        3.15111222847e-03,-9.14462130186e-06]
    fnames = 'air_g'
    argfmt = '({0:1d},{1:1d},{2:1d},{3:3.1f},{4:3g},{5:5g})'
    header = 'air3a air_g derivatives'
    testair3a_2 = Tester(funs,fargs,refs,fnames,argfmt,header=header,
        eqfun=eqfun,eqargs=eqargs,eqkeys=eqkeys)
    return (testair3a_1, testair3a_2)

def genair3b():
    """Generate air3b Testers.
    """
    import air3b
    funs = [air3b.contraction,air3b.compressibility_lemmon,air3b.cp,air3b.cv,
        air3b.density,air3b.enthalpy,air3b.entropy,air3b.expansion,
        air3b.gibbsenergy,air3b.internalenergy,air3b.kappa_s,air3b.kappa_t,
        air3b.lapserate,air3b.soundspeed,air3b.vappot]
    fargs = (0.9,300.,1e5)
    eqfun = air3b._eq_atp
    eqargs = fargs
    eqkeys = ['dhum']
    refs = [0.530280645260,0.997825670991,1246.34991644,920.600781012,
        1.09708772444,277928.954795,911.170080461,3.45704654420e-3,
        4577.93065689,186778.541048,7.41034505449e-6,1.00324517749e-5,
        7.58481752251e-4,350.719656182,193705.688576]
    fnames = ['contraction','compress_lemmon','cp','cv','density','enthalpy',
        'entropy','expansion','gibbsenergy','internalenergy','kappa_s',
        'kappa_t','lapserate','soundspeed','vappot']
    argfmt = '({0:3.1f},{1:3g},{2:5g})'
    header = 'air3b functions'
    testair3b = Tester(funs,fargs,refs,fnames,argfmt,header=header,
        eqfun=eqfun,eqargs=eqargs,eqkeys=eqkeys)
    return (testair3b,)

def genair3c():
    """Generate air3c Testers.
    """
    import air3c
    funs = air3c.air_h
    args1 = (0.9,900.,1e5)
    fargs = [(der+args1) for der in _DERS3]
    eqfun = air3c.eq_aep
    eqargs = args1
    eqkeys = ['temp','dhum']
    refs = [274592.611783,-229706.905765,297.403043058,0.903262695636,
        13963273.0104,1676.85098552,4.89772680822,0.223684689765,
        7.15703143992e-04,-6.78105152859e-06]
    fnames = 'air_h'
    argfmt = '({0:1d},{1:1d},{2:1d},{3:3.1f},{4:3g},{5:5g})'
    header = 'air3c air_h derivatives'
    testair3c_1 = Tester(funs,fargs,refs,fnames,argfmt,header=header,
        eqfun=eqfun,eqargs=eqargs,eqkeys=eqkeys)
    
    funs = [air3c.pottemp,air3c.potdensity,air3c.potenthalpy]
    fargs = (0.9,300.,5e4,1e5)
    eqfun = air3c.eq_pot
    eqargs = fargs
    eqkeys = ['dhum','tpot','dpot']
    refs = [363.653905688,0.903509489711,348872.568665]
    fnames = ['pottemp','potdensity','potenthalpy']
    argfmt = '({0:3.1f},{1:3g},{2:5g},{3:5g})'
    header = 'air3c adiabatic (potential) functions'
    testair3c_2 = Tester(funs,fargs,refs,fnames,argfmt,header=header,
        eqfun=eqfun,eqargs=eqargs,eqkeys=eqkeys)
    return (testair3c_1, testair3c_2)

def genair5():
    """Generate air5 Testers.
    """
    import air5
    funs = air5.lapserate_c100m
    fargs = (50.,20.,1e3)
    refs = 0.971588085046
    fnames = 'lapserate_c100m'
    argfmt = '({0:2g},{1:2g},{2:4g})'
    header = 'air5 lapse rate function'
    testair5 = Tester(funs,fargs,refs,fnames,argfmt,header=header)
    return (testair5,)


## Dictionary relating modules to functions
_GENDICT = {'convert0': gencnv0, 'air1': genair1, 'air2': genair2,
    'air3a': genair3a, 'air3b': genair3b, 'air3c': genair3c, 'air5': genair5}


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

