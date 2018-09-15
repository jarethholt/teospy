"""Test accuracy of the ice modules.

This module provides tests of the accuracy of the thermodynamic
functions of ice. This module can be called from the command line as

    python testice.py mod1 mod2 ...

which will run all tests on modules mod1, mod2, etc. If no arguments
are given, all available tests are run.

The functions provided by this module generate the tests for the `ice`
module of the same name. Each function returns a tuple of
:class:`~tester.Tester` instances which include the functions checked,
values of the arguments, and tables of reference values. Use the `run`
method of a Tester to run the test, and `printresults` to print a
summary.

:Available modules to test:

* :mod:`ice1`
* :mod:`ice2`

"""

__all__ = ['genice1','genice2']
from tester import Tester
_DERS2 = ((0,0),(1,0),(0,1),(2,0),(1,1),(0,2))


## Generating Tester instances
def genice1():
    """Generate ice1 Testers.
    """
    import ice1
    funs = ice1.ice_g
    args1 = (270.,1e5)
    fargs = [(der+args1) for der in _DERS2]
    refs = [-3786.74963128,1244.97335506,1.09029713624e-03,-7.67955356337,
        1.72604208697e-07,-1.27811223643e-13]
    fnames = 'ice_g'
    argfmt = '({0:1d},{1:1d},{2:3g},{3:5g})'
    header = 'ice1 ice_g derivatives'
    testice1 = Tester(funs,fargs,refs,fnames,argfmt,header=header)
    return (testice1,)

def genice2():
    """Generate ice2 Testers.
    """
    import ice2
    funs = [ice2.chempot,ice2.cp,ice2.density,ice2.enthalpy,ice2.entropy,
        ice2.expansion,ice2.helmholtzenergy,ice2.internalenergy,ice2.kappa_s,
        ice2.kappa_t,ice2.lapserate,ice2.pcoefficient,ice2.specificvolume]
    fargs = (270.,1e5)
    refs = [-3786.74963128,2073.47946211,917.181167192,-339929.555499,
        -1244.97335506,1.58309329594e-04,-3895.77934490,-340038.585212,
        1.13667916416e-10,1.17226047281e-10,2.24758128545e-08,1350462.06254,
        1.09029713624e-03]
    fnames = ['chempot','cp','density','enthalpy','entropy','expansion',
        'helmholtzenergy','internalenergy','kappa_s','kappa_t','lapserate',
        'pcoefficient','specificvolume']
    argfmt = '({0:3g},{1:5g})'
    header = 'ice2 functions'
    testice2 = Tester(funs,fargs,refs,fnames,argfmt,header=header)
    return (testice2,)


## Dictionary relating modules to functions
_GENDICT = {'ice1': genice1, 'ice2': genice2}


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

