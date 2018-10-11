"""Test accuracy of the convert modules.

This module provides tests of the conversion modules. This module can be
called from the command line as

    python testcnv.py mod ...

which will run all tests on the module `mod`. There are only two
conversion modules - convert0 and convert5. If called with no arguments,
then both modules are checked.

The functions provided by this module generate the tests for the
`convert` module of the same name. Each function returns a tuple of
:class:`~tester.Tester` instances which include the functions checked,
values of the arguments, and tables of reference values. Use the `run`
method of a Tester to run the test, and `printresults` to print a
summary.

:Available modules to test:

* :mod:`convert0`
* :mod:`convert5`

"""

__all__ = ['genconvert0','genconvert5']

import constants0
from tester import Tester

_TCELS = constants0.CELSIUS_TEMPERATURE_SI
_PATM = constants0.SEALEVEL_PRESSURE_SI


## Generating Tester instances
def genconvert0():
    """Generate convert0 Testers.
    """
    import convert0
    funs = convert0.sal_asalfrompsal
    fargs = [(35.52764437773386,201.,-21.,_PATM+1023e4),
        (35.,180.,40.,_PATM+2e7), (8.,20.,57.,_PATM)]
    refs = [35.7001299401e-3,35.1890932890e-3,8.10483771429e-3]
    refs_alt = [35.7e-3,35.1888478029e-3,8.13338057143e-3,]
    fnames = 'asalfrompsal'
    argfmt = '({0: 16f},{1:3g},{2:3g},{3:8g})'
    header = 'Convert absolute to practical salinity'
    test_afp = Tester(funs,fargs,refs,fnames,argfmt,header=header,
        refs_alt=refs_alt)
    
    funs = convert0.sal_psalfromasal
    fargs = (0.0357,201.,-21.,_PATM+1023e4)
    refs = 35.5275150654
    refs_alt = 35.5276443777
    fnames = 'psalfromasal'
    argfmt = '({0:6.4f},{1:3g},{2: 3g},{3:8g})'
    header = 'Convert practical from absolute salinity'
    test_pfa = Tester(funs,fargs,refs,fnames,argfmt,header=header,
        refs_alt=refs_alt)
    return (test_afp, test_pfa)

def genconvert5():
    """Generate convert5 Testers.
    """
    import convert5
    _TSAL1 = convert5._TSAL1
    _LON0 = convert5._LON0
    _LAT0 = convert5._LAT0
    
    p_dbar = 1023.
    p_pa = _PATM + p_dbar*convert5._DBAR2PA
    funs = convert5.cnvpressure
    fargs = [
        (p_pa,'pa','dbar'), (p_pa,'pa','torr'), (p_pa,'pa','kgf'),
        (p_pa,'pa','atm'), (p_pa,'pa','lbf/in2'), (p_pa,'pa','psi'),
        (p_dbar,'dbar','pa'), (7.74913101406e4,'torr','pa'),
        (1.05350196040e5,'kgf','pa'), (1.05350196040e5,'atm','pa'),
        (1.49842272437e3,'lbf/in2','pa'), (1.49842272437e3,'psi','pa'),
        (p_dbar,'dbar','m'), (1011.94563591,'m','dbar')
    ]
    refs = [1023.,77491.3101406,105350.196040,101.962250185,1498.42272437,
        1498.42272437,10331325.,10331325.,10331325.,10674608613.8,10331325.,
        10331325.,1011.94563591,1023.]
    refs_alt = [None,None,None,105350.196040,None,None,None,None,None,
        10331325.,None,None,None,None]  # Incorrect initial values
    fnames = 'cnvpressure'
    argfmt = '({0:16f},{1:>7},{2:>7})'
    header = 'Convert pressure units'
    test_p = Tester(funs,fargs,refs,fnames,argfmt,header=header,
        refs_alt=refs_alt)
    
    t48 = 300.
    t68 = 299.991108231
    funs = convert5.cnvtemperature
    fargs = [
        (t48,'k(t48)','degf(t48)'), (t48,'k(t48)','degf(t68)'),
        (t48,'k(t48)','degf(t90)'), (t48,'k(t48)','degc(t48)'),
        (t48,'k(t48)','degc(t68)'), (t48,'k(t48)','degc(t90)'),
        (t48,'k(t48)','k(t68)')   , (t48,'k(t48)','k(t90)')   ,
        (t68,'k(t68)','degf(t68)'), (t68,'k(t68)','degf(t90)'),
        (t68,'k(t68)','degc(t68)'), (t68,'k(t68)','degc(t90)'),
        (t68,'k(t68)','k(t90)')   , (t68,'k(t68)','degf(t90)'),
        (t68,'k(t68)','degc(t90)')
    ]
    refs = [80.33,80.3139948158,80.3018730496,26.85,26.8411082310,26.8343739165,
        299.991108231,299.984373916,80.3139916338,80.3018698685,26.8411064632,
        26.8343721491,299.984372149,80.3018698685,26.8343721491]
    fnames = 'cnvtemperature'
    argfmt = '({0:13.9f},{1:s},{2:>9})'
    header = 'Convert temperature units'
    chktol = 1e-7  # Slight change for change in significant figures
    test_t = Tester(funs,fargs,refs,fnames,argfmt,header=header,chktol=chktol)
    
    s_abs = 0.0357
    s_pss = 35.52764437773386
    unitin = 'kg/kg(abs)'
    lon = 201.
    lat = -21.
    funs = convert5.cnvsalinity
    fargs = [
        (s_abs,unitin,'pss',_TSAL1,p_pa,lon,lat),
        (s_abs,unitin,'cl',_TSAL1,p_pa,lon,lat),
        (s_abs,unitin,'kn',_TSAL1,p_pa,lon,lat),
        (s_abs,unitin,'cnd',_TCELS+25.5,p_pa,lon,lat),
        (s_abs,unitin,'kg/kg(ref)',_TSAL1,p_pa,lon,lat),
        (0.0356951724471,'kg/kg(ref)',unitin,_TSAL1,p_pa,lon,lat),
        (s_pss,'pss','cl',_TSAL1,_PATM,_LON0,_LAT0),
        (s_pss,'pss','kn',_TSAL1,_PATM,_LON0,_LAT0),
        (s_pss,'pss','cnd',_TCELS+25.5,p_pa,lon,lat),
        (s_pss,'pss','kg/kg(ref)',_TSAL1,p_pa,lon,lat),
        (s_pss,'pss',unitin,_TSAL1,p_pa,lon,lat),
        (s_pss,'pss','kn',_TSAL1,p_pa,lon,lat),
        (35.5271620502,'kn','kg/kg(abs)',_TSAL1,p_pa,lon,lat),
        (1.27540928136,'cnd','pss',_TCELS+25.5,p_pa,_LON0,_LAT0),
        (1.27540928136,'cnd','pss',_TCELS+25.5,p_pa,lon,lat)
    ]
    refs = [35.5275150654,19.6659461767,35.5270328489,1.27556269128,
        35.6950425250e-3,35.7001299401e-3,19.6660177563,35.5271620502,
        1.27556680822,35.6951724471e-3,35.7001299401e-3,35.5271620502,
        35.7001299401e-3,35.5226965418,35.5226965418]
    refs_alt = [35.5276443777,19.6660177563,35.5271620502,1.27556680822,
        35.6951724471e-3,0.0357,19.6660177563,35.5271620502,None,None,0.0357,
        35.5276443777,0.0357,None,35.7000000001e-3]
    fnames = 'cnvsalinity'
    argfmt = '({0:13f},{1:>10},{2:>10})'
    header = 'Convert salinity units'
    test_s = Tester(funs,fargs,refs,fnames,argfmt,header=header,
        refs_alt=refs_alt)
    return (test_p,test_t,test_s)


## Dictionary relating modules to functions
_GENDICT = {'convert0': genconvert0, 'convert5': genconvert5}


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

