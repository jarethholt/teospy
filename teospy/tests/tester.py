"""Compare results from a library to reference values.

This module provides the class :class:`~teospy.tests.tester.Tester` for
evaluating the numerical results from this library to reference values.
A :class:`~teospy.tests.tester.Tester` instance has as attributes the
functions to be tested; the values of the arguments to run the tests on;
a table of reference values; and several other quantities related to
printing out information. The :meth:`~teospy.tests.tester.Tester.run`
method evaluates the functions and calculates relative errors.
Afterwards, the :meth:`~teospy.tests.tester.Tester.printresults` method
can be used to display the results (if any) that are outside of a given
tolerance.
"""

__all__ = ['Tester']

import os
import sys
import warnings
import numpy

_CHKTOL = 1e-8
_NOPRINT = ('chkbnd','useext','chkvals')
_ZEROTOL = 1e-20
_RESFMT = '{:+9.2e}'
_NFMT = 9
_DERS2 = ((0,0),(1,0),(0,1),(2,0),(1,1),(0,2))
_DERS3 = ((0,0,0),(1,0,0),(0,1,0),(0,0,1),(2,0,0),(1,1,0),(1,0,1),
    (0,2,0),(0,1,1),(0,0,2))


## Auxiliary function
def _geterr(ref,result,zerotol=_ZEROTOL):
    """Calculate relative error when possible.
    
    Calculate the relative error (result/ref - 1) between a function
    result and a reference value. If the magnitude of the reference
    value is less than the given tolerance, the absolute error
    (result - ref) is returned instead.
    
    :arg ref: Reference values to compare to.
    :type ref: array(float)
    :arg result: Function results to be checked.
    :type result: array(float)
    :arg float zerotol: Tolerance to use when switching between absolute
        and relative tolerance (default _ZEROTOL).
    :returns: Error values.
    """
    cond = (numpy.abs(ref) >= zerotol)
    icond = numpy.invert(cond)
    err = numpy.zeros_like(ref)
    err[cond] = (result[cond]/ref[cond] - 1)
    err[icond] = result[icond] - ref[icond]
    return err


## Class
class Tester(object):
    """Class to aid automated testing.
    
    Class designed to aid automated testing of ``teospy``. Its
    attributes specify the names of the functions to test and the module
    they come from; the values of regular and keyword arguments to use;
    the reference values to compare to; and eventually the function
    results. It can also contain a reference to an equilibrium function
    to speed up testing higher-level modules.
    
    :arg funs: Functions to test.
    :type funs: function or list[function]
    :arg fargs: Standard input values for the functions.
    :type fargs: tuple(float) or list[tuple(float)]
    :arg refs: Reference values for the functions.
    :type refs: list[float] or list[list[float]]
    :arg fnames: Names of the functions to test. Used in printing.
    :type fnames: str or list[str]
    :arg str argfmt: Format string for the function arguments. Used in
        printing.
    :arg header: Header line to print at the start of the summary or
        None (default) to print nothing.
    :type header: str or None
    :arg fkwargs: Keyword arguments to pass to the functions. Use None
        (default) for no additional arguments.
    :type fkwargs: dict or None
    :arg eqfun: Function to calculate primary variable values. Used by
        modules level 3 and higher. Use None (default) for lower-level
        modules.
    :type eqfun: function or None
    :arg eqargs: Values to pass to the equilibrium function. Used by
        modules level 3 and higher. Use None (default) for lower-level
        modules.
    :type eqargs: tuple(float) or None
    :arg eqkwargs: Keyword arguments to pass to the equilibrium
        function. Used by modules level 3 and higher.
    :type eqkwargs: dict or None
    :arg eqkeys: Names of the keyword arguments for primary variables.
        Used by modules level 3 and higher. Use None (default) for
        lower-level modules.
    :type eqkeys: list[str] or None
    :arg keepkeys: Names of keyword arguments to actually pass on to
        other functions. Use None (default) to pass on all keys.
    :type keepkeys: list[str] or None
    :arg refs_alt: Alternative reference values for the functions, if
        any (default None). These are used to keep track of older
        reference values before new results are vetted.
    :type refs_alt: list[float or None] or None
    :arg chktol: A tolerance to default to for the `printresults`
        method that will be used if an explicit value is not passed to
        the method. If None, _CHKTOL is used.
    :type chktol: float or None
    :raises RuntimeWarning: If the shape of the reference values does
        not match the shape of the functions and arguments lists.
    
    :Additional attributes:
    
    * nfun (int): Number of functions to test.
    * narg (int): Number of arguments to test the functions on.
    * nvar (int): Number of variables the functions accept.
    * results (array(float)): Results from the functions of this
        library; created after using the `run` method.
    * errs (array(float)): Error between the results and reference
        values; relative error when the reference value is > _ZEROTOL.
        Created after using the `run` method.
    * nstr (int): Length of strings used to print function names.
    * nfmt (int): Length of strings used to print argument values.
    """
    
    def __init__(self,funs,fargs,refs,fnames,argfmt,header=None,
        fkwargs=None,eqfun=None,eqargs=None,eqkwargs=None,eqkeys=None,
        keepkeys=None,refs_alt=None,chktol=None):
        if isinstance(funs,list):
            self.funs = funs
            nfun = len(funs)
        else:
            self.funs = [funs]
            nfun = 1
        self.nfun = nfun
        
        if isinstance(fargs,list):
            self.fargs = fargs
            narg = len(fargs)
            nvar = len(fargs[0])
        else:
            self.fargs = [fargs]
            narg = 1
            nvar = len(fargs)
        self.narg = narg
        self.nvar = nvar
        
        refsarr = numpy.atleast_2d(numpy.array(refs))
        if (narg == 1 and nfun > 1):
            refsarr = refsarr.T
        self.refs = refsarr
        self.results = None
        self.errs = None
        
        # Include alternative reference values
        if refs_alt is not None:
            altrefsarr = numpy.atleast_2d(numpy.array(refs_alt,dtype=float))
            if (narg == 1 and nfun > 1):
                altrefsarr = altrefsarr.T
            cond = numpy.isnan(altrefsarr)
            altrefsarr[cond] = refsarr[cond]
            self.refs_alt = altrefsarr
        
        # Check that the shape of fnames matches funs
        if nfun == 1:
            assert isinstance(fnames,str)
            self.fnames = [fnames]
        else:
            self.fnames = fnames
        self.nstr = max(len(fname) for fname in self.fnames)
        self.argfmt = argfmt
        args = self.fargs[0]
        if fkwargs is not None:
            for (key,val) in fkwargs.items():
                args += (key,val)
        self.nfmt = len(argfmt.format(*args))
        self.header = header
        self.fkwargs = fkwargs
        self.chktol = chktol
        
        # Include equilibrium function and keywords
        self.eqfun = eqfun
        self.eqargs = eqargs
        self.eqkwargs = eqkwargs
        self.eqkeys = eqkeys
        self.keepkeys = keepkeys
        return None
    
    def run(self,zerotol=_ZEROTOL):
        """Run all the functions in the test.
        """
        kwargs = dict()
        if self.fkwargs is not None:
            kwargs.update(self.fkwargs)
        # Run equilibrium function if it exists
        if self.eqfun is None:
            eqdict = dict()
        else:
            if self.eqkwargs is not None:
                eqkwargs = self.eqkwargs
            else:
                eqkwargs = dict()
            
            with warnings.catch_warnings():
                warnings.filterwarnings('ignore','Step sizes are smaller',
                    RuntimeWarning)
                warnings.filterwarnings('ignore','Maximum number of iterations',
                    RuntimeWarning)
                eqvals = self.eqfun(*self.eqargs,**eqkwargs)
            
            if isinstance(eqvals,float):
                eqvals = (eqvals,)
            if self.keepkeys is None:
                keep = self.eqkeys
            else:
                keep = self.keepkeys
            eqdict = {key: val for (key,val) in zip(self.eqkeys,eqvals)
                if key in keep}
            self.eqdict = eqdict
        kwargs.update(eqdict)
        
        results = numpy.zeros((self.nfun,self.narg),dtype=float)
        with warnings.catch_warnings():
            warnings.filterwarnings('ignore','Step sizes are smaller',
                RuntimeWarning)
            warnings.filterwarnings('ignore','Maximum number of iterations',
                RuntimeWarning)
            for ifun in range(self.nfun):
                fun = self.funs[ifun]
                for iarg in range(self.narg):
                    farg = self.fargs[iarg]
                    res = fun(*farg,**kwargs)
                    results[ifun,iarg] = res
        self.results = results
        self.errs = _geterr(self.refs,results,zerotol=zerotol)
        return None
    
    def printresults(self,chktol=None):
        """Print results outside of a tolerance.
        """
        if self.results is None:
            errmsg = ('Results have not been calculated yet. Use '
                'Tester.run to get function results.')
            raise ValueError(errmsg)
        if chktol is None:
            if self.chktol is not None:
                chktol = self.chktol
            else:
                chktol = _CHKTOL
        if numpy.max(numpy.abs(self.errs)) <= chktol:
            okmsg = 'All results within tolerance {0}'.format(chktol)
            if self.header is not None:
                okmsg = self.header + ': ' + okmsg
            print(okmsg)
            return None
        
        # Print only results that are out of range
        msg = 'Results outside of tolerance {0}'.format(chktol)
        if self.header is not None:
            msg = self.header + ': ' + msg
        print(msg)
        nstr = max(self.nstr,len('Function'))
        lead = 'Function'.rjust(nstr) + '  '
        nfmt = self.nfmt
        lead += 'Args'.ljust(nfmt) + '  '
        lead += 'Ref'.ljust(_NFMT) + '  '
        lead += 'Result'.ljust(_NFMT) + '  '
        lead += 'Error'.ljust(_NFMT)
        print(lead)
        
        for ifun in range(self.nfun):
            for iarg in range(self.narg):
                if abs(self.errs[ifun,iarg]) <= chktol:
                    continue
                fname = self.fnames[ifun]
                farg = self.fargs[iarg]
                if self.fkwargs is not None:
                    for (key,val) in self.fkwargs.items():
                        farg += (key,val)
                ref = self.refs[ifun,iarg]
                res = self.results[ifun,iarg]
                err = self.errs[ifun,iarg]
                line = fname.rjust(nstr) + '  '
                line += self.argfmt.format(*farg) + '  '
                line += '  '.join(_RESFMT.format(val)
                    for val in (ref,res,err))
                print(line)
        return None

