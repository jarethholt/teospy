"""Test accuracy of the GSW module.

This module provides tests of the accuracy of the Gibbs SeaWater (GSW)
toolbox functions. This module can be called from the command line as::

    python testgsw.py

which will run the tests and print a summary.

This module provides one function to generate
:class:`~teospy.tests.tester.Tester` instances. Each instance includes
the functions checked, values of the arguments, and tables of reference
values. Use the ``run`` method of a Tester to run the test, and
``printresults`` to print a summary.
"""

__all__ = ['gengsw5']

from teospy.tests.tester import Tester, _DERS3


## Generating Tester instances
def gengsw5():
    """Generate gsw5 Testers.
    """
    from teospy import gsw5
    funs = gsw5.asalfrompsal
    fargs = [(35.527515065427778,201.,-21.,1023.), (35.,180.,40.,2e3),
        (8.,20.,57.,0.)]
    refs = [35.7,35.1890932890,8.10483771429]
    refs_alt = [None,35.1888478029,8.13338057143]
    fnames = 'asalfrompsal'
    argfmt = '({0:7.4f},{1:3g},{2:3g},{3:4g})'
    header = 'GSW abs salinity'
    test_sal1 = Tester(funs,fargs,refs,fnames,argfmt,header=header,
        refs_alt=refs_alt)
    
    funs = gsw5.psalfromasal
    fargs = (35.7,201.,-21.,1023.)
    refs = 35.527515065427778
    fnames = 'psalfromasal'
    argfmt = '({0:4.1f},{1:3g},{2:3g},{3:4g})'
    header = 'GSW pract salinity'
    test_sal2 = Tester(funs,fargs,refs,fnames,argfmt,header=header)
    
    funs = gsw5.gsw_g
    args1 = (35.,26.85,0.)
    fargs = [(der+args1) for der in _DERS3]
    refs = [-5113.70064124,78.5928261339,-374.452000830,9.77858058750e-4,
        2.24755137017,0.789935187192,-7.16680931996e-7,-13.3358337534,
        3.04607508052e-7,-4.10939723950e-13]
    fnames = 'gsw_g'
    argfmt = '({0:1d},{1:1d},{2:1d},{3:2g},{4:5.2f},{5:1g})'
    header = 'GSW Gibbs function'
    test_der = Tester(funs,fargs,refs,fnames,argfmt,header=header)
    
    funs = [gsw5.alpha_t,gsw5.beta_t,gsw5.cp,gsw5.density,gsw5.enthalpy,
        gsw5.entropy,gsw5.kappa,gsw5.kappa_t,gsw5.specvol,gsw5.svel]
    fargs = (35.7,25.5,1023.)
    refs = [3.09837839319e-4,7.25729797839e-4,3974.42541260,1027.95249316,
        110776.712409,352.818797715,4.03386268546e-6,4.10403794615e-6,
        9.72807602158e-4,1552.93372863]
    fnames = ['alpha_t','beta_t','cp','density','enthalpy','entropy','kappa',
        'kappa_t','specvol','svel']
    argfmt = '({0:4.1f},{1:4.1f},{2:4g})'
    header = 'GSW in-situ functions'
    test_t = Tester(funs,fargs,refs,fnames,argfmt,header=header)
    
    funs = [gsw5.pottemp,gsw5.potdens]
    fargs = (35.7,25.5,1023.,0.)
    refs = [25.2720983155,1023.66254941]
    fnames = ['pottemp','potdens']
    argfmt = '({0:4.1f},{1:4.1f},{2:4g},{3:1g})'
    header = 'GSW potential functions 1'
    eqfun = gsw5.pottemp
    eqkeys = ['tp_cels']
    test_tp1 = Tester(funs,fargs,refs,fnames,argfmt,header=header,
        eqfun=eqfun,eqargs=fargs,eqkeys=eqkeys)
    
    funs = gsw5.pottemp
    fargs = (35.7,25.2720983155409,0.,1023.)
    refs = 25.5
    fnames = 'pottemp'
    argfmt = '({0:4.1f},{1:7.4f},{2:1g},{3:4g})'
    header = 'GSW potential functions 2'
    test_tp2 = Tester(funs,fargs,refs,fnames,argfmt,header=header)
    
    funs = gsw5.tconfromtpot
    fargs = (35.7,25.5)
    refs = 25.4805463842
    fnames = 'tconfromtpot'
    argfmt = '({0:4.1f},{1:4.1f})'
    header = 'GSW potential to conservative temp'
    test_p2c = Tester(funs,fargs,refs,fnames,argfmt,header=header)
    
    funs = gsw5.tpotfromtcon
    fargs = (35.7,25.4805463842239)
    refs = 25.5
    fnames = 'tpotfromtcon'
    argfmt = '({0:4.1f},{1:7.4f})'
    header = 'GSW conservative to potential temp'
    test_c2p = Tester(funs,fargs,refs,fnames,argfmt,header=header)
    
    funs = [gsw5.alpha_tcon,gsw5.alpha_tpot,gsw5.beta_tcon,gsw5.beta_tpot,
        gsw5.cabb_tcon_alt,gsw5.cabb_tpot_alt,gsw5.thrmb_tcon_alt,
        gsw5.thrmb_tpot_alt]
    fargs = (35.,20.,1e3)
    refs = [2.69414197935e-4,2.69750083750e-4,7.23213810967e-4,7.31582585561e-4,
        9.16324749728e-6,8.76597878464e-6,1.72728502693e-8,1.70965495606e-8]
    refs_alt = [2.69418609861e-4,2.69753733317e-4,7.23213672954e-4,
        7.31582583383e-4,8.96907383083e-6,8.75963154048e-6,1.72708365652e-8,
        1.70945045984e-8]
    fnames = ['alpha_tcon','alpha_tpot','beta_tcon','beta_tpot','cabb_tcon_alt',
        'cabb_tpot_alt','thrmb_tcon_alt','thrmb_tpot_alt']
    argfmt = '({0:2g},{1:2g},{2:4g})'
    header = 'GSW expansion/contraction coeffs'
    eqfun = gsw5.pottemp
    eqargs = fargs + (0.,)
    eqkeys = ['tp_cels']
    test_ec1 = Tester(funs,fargs,refs,fnames,argfmt,header=header,
        eqfun=eqfun,eqargs=eqargs,eqkeys=eqkeys,refs_alt=refs_alt)
    
    funs = [gsw5.cabb_tcon,gsw5.cabb_tpot,gsw5.thrmb_tcon,gsw5.thrmb_tpot]
    fargs = (35.,20.,1e3)
    refs = [9.16324749728e-6,8.76597878464e-6,1.72728502693e-8,1.70965495606e-8]
    fnames = ['cabb_tcon','cabb_tpot','thrmb_tcon','thrmb_tpot']
    argfmt = '({0:2g},{1:2g},{2:4g})'
    header = 'GSW cabb/thrmb with exact derivs'
    eqfun = gsw5.pottemp
    eqargs = fargs + (0.,)
    eqkeys = ['tp_cels']
    test_ec2 = Tester(funs,fargs,refs,fnames,argfmt,header=header,eqfun=eqfun,
        eqargs=eqargs,eqkeys=eqkeys)
    return (test_sal1, test_sal2, test_der, test_t, test_tp1, test_tp2,
        test_p2c, test_c2p, test_ec1, test_ec2)


## See if all values fall within the given tolerances
if __name__ == "__main__":
    testlist = gengsw5()
    for test in testlist:
        test.run()
        test.printresults()

