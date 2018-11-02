"""Test accuracy of the sal modules.

This module provides tests of the accuracy of the thermodynamic
functions for salt in seawater. This module can be called from the
command line as::

    python testsal.py useext mod1 mod2 ...

which will run all tests on modules ``mod1``, ``mod2``, etc. If no
arguments are given, all available tests are run. The first argument can
also be either ``ext`` or ``noext`` to check only values with or without
the high- temperature, high-salinity extension, or ``all`` to check both
(the default behavior).

The functions provided by this module generate the tests for the ``sal``
module of the same name. Each function returns a tuple of
:class:`~tester.Tester` instances which include the functions checked,
values of the arguments, and tables of reference values. Use the ``run``
method of a Tester to run the test, and ``printresults`` to print a
summary.

:Available modules to test:

* :mod:`convert0`
* :mod:`sal1`
* :mod:`sal2`

"""

__all__ = ['gencnv0','gensal1','gensal1e','gensal2']
from tester import Tester
_DERS2 = ((0,0),(1,0),(0,1),(2,0),(1,1),(0,2))
_DERS3 = ((0,0,0),(1,0,0),(0,1,0),(0,0,1),(2,0,0),(1,1,0),(1,0,1),
    (0,2,0),(0,1,1),(0,0,2))


## Generating Tester instances
def gencnv0():
    """Generate convert0 Testers.
    """
    import convert0
    funs = convert0.sal_molality
    fargs = (0.035,)
    refs = 1.15493681893
    fnames = 'sal_molality'
    argfmt = '({0:5.3f})'
    header = 'convert0 sal functions'
    testcnv0 = Tester(funs,fargs,refs,fnames,argfmt,header=header)
    return (testcnv0,)

def gensal1():
    """Generate sal1 Testers without extension.
    """
    import sal1
    NSALTERMS = sal1.NSALTERMS
    funs = sal1.sal_g_term
    args1 = (300.,1e6)
    fargs = [((iterm,)+der+args1) for iterm in range(1,NSALTERMS+1)
        for der in _DERS2]
    refs = [79427.9694846,264.759898282,0.,0.,0.,0.,301223.934546,1558.10730393,
        -7.60101106780e-04,19.0397747694,9.35496617545e-07,1.28077803426e-12,
        -345570.781498,-1749.19911587,2.01608449881e-04,-9.32088823620,
        2.17412289837e-08,-4.78063827318e-13,1468073.64559,7741.24404962,
        -4.33952339916e-04,-6.13689642370,-3.50876195080e-06,-6.06204383305e-13,
        -3776969.31546,-15135.6522248,1.11272425476e-03,0.,0.,0.,6151235.69234,
        14157.0509933,0.,0.,0.,0.,-3734033.38866,0.,0.,0.,0.,0.]
    fnames = 'sal_g_term'
    argfmt = '({0:1d},{1:1d},{2:1d},{3:3g},{4:5g})'
    header = 'sal1 derivatives noext'
    testsal1 = Tester(funs,fargs,refs,fnames,argfmt,header=header)
    return (testsal1,)

def gensal1e():
    """Generate sal1 Testers with extension.
    """
    import sal1
    NSALTERMS = sal1.NSALTERMS
    funs = sal1.sal_g_term
    args1 = (300.,1e6)
    fargs = [((iterm,)+der+args1) for iterm in range(1,NSALTERMS+1)
        for der in _DERS2]
    refs = [79427.9694845756,264.759898281919,0.,0.,0.,0.,301219.100652623,
        1558.00684228239,-7.65480018679231e-4,19.0654571217448,
        8.2370799190862e-7,1.28077803425625e-12,-345511.73639383,
        -1748.27447798359,2.67310849508517e-4,-9.53134759747051,
        1.05063141157166e-6,-4.78063827317764e-13,1467897.23729704,
        7738.95481914288,-6.30250546547853e-4,-5.74817212818652,
        -6.05610162524257e-6,-6.06204383304926e-13,-3776969.31546023,
        -15135.6522247701,1.11272425476394e-3,0.,0.,0.,6151235.69234474,
        14157.050993291,0.,0.,0.,0.,-3734033.38866189,0.,0.,0.,0.,0.]
    fnames = 'sal_g_term'
    argfmt = '({0:1d},{1:1d},{2:1d},{3:3g},{4:5g})'
    header = 'sal1 derivatives ext'
    fkwargs = {'useext': True}
    testsal1e = Tester(funs,fargs,refs,fnames,argfmt,header=header,
        fkwargs=fkwargs)
    return (testsal1e,)

def gensal2():
    """Generate sal2 Testers.
    """
    import sal2
    funs = sal2.sal_g
    args1 = (0.035,300.,1e6)
    fargs = [(der+args1) for der in _DERS3]
    refs = [127.033640309,77949.2100395,18.6360407073,-2.55600080319e-05,
        2248200.54660,790.563810558,-7.15686520588e-04,0.597842170749,
        2.85865076269e-08,4.09543164905e-14]
    fnames = 'sal_g'
    argfmt = '({0:1d},{1:1d},{2:1d},{3:5.3f},{4:3g},{5:5g})'
    header = 'sal2 sal_g derivatives'
    testsal2_1 = Tester(funs,fargs,refs,fnames,argfmt,header=header)
    
    funs = [sal2.actcoeff,sal2.activityw,sal2.actpotential,sal2.chemcoeff,
        sal2.dilution,sal2.liqpot,sal2.osmcoeff,sal2.salpot,sal2.saltenthalpy,
        sal2.saltentropy,sal2.saltvolume]
    fargs = (0.035,300.,1e6)
    refs = [-0.527003008913,0.981388410188,-0.429940465498,2754.04566958,
        78687.0191309,-2601.18871107,0.902937456585,77949.2100395,
        -156107.959196,-532.458305922,-7.30285943768e-04]
    fnames = ['actcoeff','activityw','actpotential','chemcoeff','dilution',
        'liqpot','osmcoeff','salpot','saltenthalpy','saltentropy','saltvolume']
    argfmt = '({0:5.3f},{1:3g},{2:5g})'
    header = 'sal2 functions'
    testsal2_2 = Tester(funs,fargs,refs,fnames,argfmt,header=header)
    
    funs = [sal2.mixenthalpy,sal2.mixentropy,sal2.mixvolume]
    fargs = (0.01,0.035,0.6,300.,1e6)
    refs = [16.1539911284,0.966829422617,-5.94174956892e-08]
    fnames = ['mixenthalpy','mixentropy','mixvolume']
    argfmt = '({0:4.2f},{1:5.3f},{2:3.1f},{3:3g},{4:5g})'
    header = 'sal2 mix functions'
    testsal2_3 = Tester(funs,fargs,refs,fnames,argfmt,header=header)
    return (testsal2_1,testsal2_2,testsal2_3)


## Dictionaries relating modules to functions
_GENDICT_NOEXT = {'convert0': gencnv0, 'sal1': gensal1, 'sal2': gensal2}
_GENDICT_EXT = {'sal1': gensal1e}


## See if all values fall within the given tolerances
if __name__ == "__main__":
    # Figure out which dictionaries to include
    import sys
    testlist = list()
    if len(sys.argv) == 1:
        for (modname,genfun) in _GENDICT_NOEXT.items():
            testlist += list(genfun())
        for (modname,genfun) in _GENDICT_EXT.items():
            testlist += list(genfun())
    elif len(sys.argv) == 2 and (sys.argv[1] in ('noext','ext','all')):
        if sys.argv[1] != 'ext':
            for (modname,genfun) in _GENDICT_NOEXT.items():
                testlist += list(genfun())
        if sys.argv[1] != 'noext':
            for (modname,genfun) in _GENDICT_EXT.items():
                testlist += list(genfun())
    else:
        if sys.argv[1] in ('noext','ext','all'):
            useext = sys.argv[1]
            mods = sys.argv[2:]
        else:
            useext = 'all'
            mods = sys.argv[1:]
        for mod in mods:
            if useext != 'ext':
                for (modname,genfun) in _GENDICT_NOEXT.items():
                    testlist += list(genfun())
            if useext != 'noext':
                for (modname,genfun) in _GENDICT_EXT.items():
                    testlist += list(genfun())
    
    # Run tests
    for test in testlist:
        test.run()
        test.printresults()

