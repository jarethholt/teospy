"""Test accuracy of the flu modules.

This module provides tests of the fluid water thermodynamic functions.
This module can be called from the command line as::

    python testflu.py mod1 mod2 ...

which will run all tests whose modules match one of ``mod1``, ``mod2``,
etc. For example, one can use ``flu3`` to check all of the level 3
modules (``flu3a``, ``flu3b``). If no arguments are given, all available
tests are run.

The functions provided by this module generate the tests for the ``flu``
module of the same name. Each function returns a tuple of
:class:`~teospy.tests.tester.Tester` instances which include the
functions checked, values of the arguments, and tables of reference
values. Use the ``run`` method of a Tester to run the test, and
``printresults`` to print a summary.

:Available modules to test:

* :mod:`~teospy.flu1`
* :mod:`~teospy.flu2`
* :mod:`~teospy.flu3a`
* :mod:`~teospy.flu3b`
* :mod:`~teospy.liq5_f03`

"""

__all__ = ['genflu1','genflu2','genflu3a','genflu3b','genliq5']

from teospy.tests.tester import Tester, _DERS2


## Generating Tester instances
def genflu1():
    """Generate flu1 Testers.
    """
    from teospy import flu1
    funs = flu1.flu_f
    args1 = (300.,1e3)
    fargs = [(der+args1) for der in _DERS2]
    refs = [-5351.74115204,-390.904170767,7.83300135597,-13.6840204925,
        0.639359046588,2.24824656167]
    fnames = 'flu_f'
    argfmt = '({0:1d},{1:1d},{2:3g},{3:3g})'
    header = 'flu1 flu_f derivatives'
    testflu1 = Tester(funs,fargs,refs,fnames,argfmt,header=header)
    return (testflu1,)

def genflu2():
    """Generate flu2 Testers.
    """
    from teospy import flu2
    funs = [flu2.cp,flu2.cv,flu2.enthalpy,flu2.entropy,flu2.expansion,
        flu2.gibbsenergy,flu2.internalenergy,flu2.kappa_s,flu2.kappa_t,
        flu2.lapserate,flu2.pressure,flu2.soundspeed]
    fargs = (300.,1e3)
    refs = [4159.37519963,4105.20614776,119752.511434,390.904170767,
        2.82413312531e-04,2481.26020392,111919.510078,4.35960581171e-10,
        4.41713172024e-10,2.03694039833e-08,7833001.35597,1514.52479780]
    fnames = ['cp','cv','enthalpy','entropy','expansion','gibbsenergy',
        'internalenergy','kappa_s','kappa_t','lapserate','pressure',
        'soundspeed']
    argfmt = '({0:3g},{1:3g})'
    header = 'flu2 functions'
    testflu2 = Tester(funs,fargs,refs,fnames,argfmt,header=header)
    return (testflu2,)

def genflu3a():
    """Generate flu3a Testers.
    """
    from teospy import flu3a
    funs = flu3a.eq_tp_liq
    argsliq = (300.,1e5)
    fargs = argsliq
    refs = 996.556340389
    fnames = 'eq_tp_liq'
    argfmt = '({0:3g},{1:3g})'
    header = 'flu3a liquid equilibrium'
    testflu3a_1 = Tester(funs,fargs,refs,fnames,argfmt,header=header)
    
    funs = flu3a.eq_tp_vap
    argsvap = (300.,1e3)
    fargs = argsvap
    refs = 7.22603510025e-03
    fnames = 'eq_tp_vap'
    argfmt = '({0:3g},{1:3g})'
    header = 'flu3a vapour equilibrium'
    testflu3a_2 = Tester(funs,fargs,refs,fnames,argfmt,header=header)
    
    funs = flu3a.liq_g
    fargs = [(der+argsliq) for der in _DERS2]
    eqfun = flu3a.eq_tp_liq
    eqkeys = ['dliq']
    eqargs = argsliq
    refs = [-5265.05045577,-393.062433815,1.00345555938e-03,-13.9354650734,
        2.75753316815e-07,-4.52072086722e-13]
    fnames = 'liq_g'
    argfmt = '({0:1d},{1:1d},{2:3g},{3:3g})'
    header = 'flu3a liq_g derivatives'
    testflu3a_3 = Tester(funs,fargs,refs,fnames,argfmt,header=header,
        eqfun=eqfun,eqargs=eqargs,eqkeys=eqkeys)
    
    funs = flu3a.vap_g
    fargs = [(der+argsvap) for der in _DERS2]
    eqfun = flu3a.eq_tp_vap
    eqkeys = ['dvap']
    eqargs = argsvap
    refs = [-180090.341338,-9103.67940087,138.388478069,-6.24707163427,
        0.462704658818,-0.138455798864]
    fnames = 'vap_g'
    argfmt = '({0:1d},{1:1d},{2:3g},{3:3g})'
    header = 'flu3a vap_g derivatives'
    testflu3a_4 = Tester(funs,fargs,refs,fnames,argfmt,header=header,
        eqfun=eqfun,eqargs=eqargs,eqkeys=eqkeys)
    return (testflu3a_1,testflu3a_2,testflu3a_3,testflu3a_4)

def genflu3b():
    """Generate flu3b Testers.
    """
    from teospy import flu3b
    funs = [flu3b.liq_cp,flu3b.liq_cv,flu3b.liq_enthalpy,flu3b.liq_entropy,
        flu3b.liq_expansion,flu3b.liq_gibbsenergy,flu3b.liq_internalenergy,
        flu3b.liq_kappa_s,flu3b.liq_kappa_t,flu3b.liq_lapserate,
        flu3b.liq_soundspeed]
    fargs = (300.,1e5)
    eqfun = flu3b._eq_tp_liq
    eqargs = fargs
    eqkeys = ['dliq']
    refs = [4180.63952202,4130.17861503,112653.679689,393.062433815,
        2.74803716256e-4,-5265.05045577,112553.334133,4.45077521253e-10,
        4.50515304336e-10,1.97878804448e-8,1501.52041506]
    fnames = ['cp','cv','enthalpy','entropy','expansion','gibbsenergy',
        'internalenergy','kappa_s','kappa_t','lapserate','soundspeed']
    argfmt = '({0:3g},{1:3g})'
    header = 'flu3b liq functions'
    testflu3b_1 = Tester(funs,fargs,refs,fnames,argfmt,header=header,
        eqfun=eqfun,eqargs=eqargs,eqkeys=eqkeys)
    
    funs = [flu3b.vap_cp,flu3b.vap_cv,flu3b.vap_enthalpy,flu3b.vap_entropy,
        flu3b.vap_expansion,flu3b.vap_gibbsenergy,flu3b.vap_internalenergy,
        flu3b.vap_kappa_s,flu3b.vap_kappa_t,flu3b.vap_lapserate,
        flu3b.vap_soundspeed]
    fargs = (300.,1e3)
    eqfun = flu3b._eq_tp_vap
    eqargs = fargs
    eqkeys = ['dvap']
    refs = [1874.12149028,1410.22845789,2551013.47892,9103.67940087,
        3.34352010567e-3,-180090.341338,2412625.00085,7.52840457971e-4,
        1.00048646242e-3,7.40674488635e-2,428.744430495]
    fnames = ['cp','cv','enthalpy','entropy','expansion','gibbsenergy',
        'internalenergy','kappa_s','kappa_t','lapserate','soundspeed']
    argfmt = '({0:3g},{1:3g})'
    header = 'flu3b vap functions'
    testflu3b_2 = Tester(funs,fargs,refs,fnames,argfmt,header=header,
        eqfun=eqfun,eqargs=eqargs,eqkeys=eqkeys)
    return (testflu3b_1,testflu3b_2)

def genliq5():
    """Generate liq5_f03 Testers.
    """
    from teospy import liq5_f03
    funs = liq5_f03.liq_g
    args1 = (300.,1e5)
    fargs = [(der+args1) for der in _DERS2]
    refs = [-5265.05056073,-393.062597709,0.100345554745e-2,-13.9354762020,
        0.275754520492e-6,-0.452067557155e-12]
    fnames = 'liq_g'
    argfmt = '({0:3d},{1:3d})'
    header = 'Feistel 2003 g derivatives'
    testliq5_1 = Tester(funs,fargs,refs,fnames,argfmt,header=header)
    
    funs = [liq5_f03.cp,liq5_f03.density,liq5_f03.expansion,liq5_f03.kappa_t,
        liq5_f03.soundspeed]
    fargs = args1
    refs = [4180.64286060,996.556352243,0.274804919056e-3,0.450510795725e-9,
        1501.52808421]
    fnames = ['cp','density','expansion','kappa_t','soundspeed']
    argfmt = '({0:3g},{1:3g})'
    header = 'Feistel 2003 functions'
    testliq5_2 = Tester(funs,fargs,refs,fnames,argfmt,header=header)
    return (testliq5_1,testliq5_2)


## Dictionary relating modules to functions
_GENDICT = {'flu1': genflu1, 'flu2': genflu2, 'flu3a': genflu3a,
    'flu3b': genflu3b, 'liq5': genliq5}


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

