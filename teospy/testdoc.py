"""Test floating-point values in docstrings.

This module aids in using `doctest` on this library. The docstrings of
many functions contain reference values with a limited number of
significant digits, an amount that should be replicated regardless of
implementation (Python, Fortran, VB, etc.). In addition, the reference
values have a variety of formats.

This module provides a custom OutputChecker, a class used by `doctest`
to evaluate whether the expected and actual outputs match. This new
checker, :class:`fltOutputChecker`, standardizes the formatting of the
floats (including number of significant figures) before comparing the
strings.

This module can also be called as

    python testdoc.py [-v] modname

which will run `doctest` on the module `modname`. The optional flag
`-v` can be used to make the test runner verbose. If not verbose
(default) then only failed tests are printed, or a short message if all
tests were passed.
"""

__all__ = ['FltOutputChecker','testmod_flt']
import doctest


def _extractsigfmt(numstr):
    """Find a standard float formatter from a string.
    
    From a string representing a float, find the number of significant
    digits and return a string formatter to standardize its
    representation.
    """
    # Simple necessary exception: some zeroes are really zero
    if numstr == '0.0\n':
        fmt = '{:3.1f}'
        return fmt
    
    # Pull out only the digits part of an exponential-format float
    digs = numstr.lower().rstrip('\n').split('e')[0]
    # Separate decimal point, remove leading sign and zeroes
    res = digs.split('.')
    if len(res) == 1:
        l, r = res, ''
    elif len(res) > 2:
        errmsg = ('Multiple decimal points in input string {:s}').format(numstr)
        raise ValueError(errmsg)
    else:
        l, r = res
    l = l.lstrip(' +-0')
    if l == '':
        r = r.lstrip('0')
    nos = l + r
    
    # What's left should only be the significant numbers
    if not nos.isnumeric():
        errmsg = ('Splitting did not result in the expected format. '
            'Wanted only numbers, got {:s}').format(nos)
    nsig = len(nos)
    
    # Create formatter from number of significant digits
    fmt = '{:+' + '{0:d}.{1:d}'.format(nsig+6,nsig-1) + 'e}'
    return fmt

class FltOutputChecker(doctest.OutputChecker):
    """OutputChecker customized for floats.
    
    This OutputChecker overrides the standard check_output function to
    convert the expected and actual output to a common exponential
    format with the number of significant figures given by the expected
    output.
    """
    def check_output(self,want,got,optionflags):
        """Check that expected and actual output agree within sig figs.
        """
        # Exception: Sometimes empty strings are passed
        if want == '' or got == '':
            return True
        fmt = _extractsigfmt(want)
        want_std = fmt.format(float(want))
        got_std = fmt.format(float(got))
        return (want_std == got_std)

def testmod_flt(mod,verbose=False):
    """Doctest a module using float comparisons.
    """
    finder = doctest.DocTestFinder()
    dtlist = finder.find(mod)
    checker = FltOutputChecker()
    runner = doctest.DocTestRunner(checker)
    for dt in dtlist:
        runner.run(dt)
    results = runner.summarize(verbose=verbose)
    return results


if __name__ == '__main__':
    # Run a doctest on the named module
    import importlib
    import sys
    
    verbose = False
    if sys.argv[1] == '-v':
        verbose = True
    mod = importlib.import_module(sys.argv[-1])
    results = testmod_flt(mod,verbose=verbose)
    if results.failed == 0 and results.attempted > 0:
        okmsg = 'All {0:d} tests passed'.format(results.attempted)
        print(okmsg)

