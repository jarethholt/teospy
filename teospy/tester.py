"""Compare results from a library to reference values.

This module provides the class Tester for evaluating the numerical
results from this library to reference values. A Tester instance has as
attributes the functions to be tested; the values of the arguments to
run the tests on; a table of reference values; and several other
quantities related to printing out information. The `run` method
evaluates the functions and calculates relative errors. Afterwards, the
`printresults` method can be used to display the results (if any) that
are outside of a given tolerance.
"""

__all__ = ['Tester']

import warnings
import numpy
import constants0

_CHKTOL = constants0.CHKTOL
_NOPRINT = ('chkbnd','useext','chkvals')
_ZEROTOL = 1e-20
_RESFMT = '{:+9.2e}'
_NFMT = 9


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
    
    Class designed to aid automated testing of this library. Its
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
    :arg eqkws: Names of the keyword arguments for primary variables.
        Used by modules level 3 and higher. Use None (default) for
        lower-level modules.
    :type eqkws: list[str] or None
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
    
    :Methods:
    
    * run: Run all the functions in the Tester, record results, and
        calculate error values.
    * printresults: Print a summary with the function names, arguments,
        result, reference value, and error for all results outside of a
        given tolerance.
    
    """
    
    def __init__(self,funs,fargs,refs,fnames,argfmt,header=None,
        fkwargs=None,eqfun=None,eqargs=None,eqkws=None):
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
        
        # Check that the shape of fnames matches funs
        if nfun == 1:
            assert isinstance(fnames,str)
            self.fnames = [fnames]
        else:
            self.fnames = fnames
        self.nstr = max(len(fname) for fname in self.fnames)
        self.argfmt = argfmt
        self.nfmt = len(argfmt.format(*self.fargs[0]))
        self.header = header
        self.fkwargs = fkwargs
        
        # Include equilibrium function and keywords
        self.eqfun = eqfun
        self.eqargs = eqargs
        self.eqkws = eqkws
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
            eqvals = self.eqfun(*self.eqargs)
            if isinstance(eqvals,float):
                eqvals = (eqvals,)
            eqdict = {kw: val for (kw,val) in zip(self.eqkws,eqvals)}
            self.eqdict = eqdict
        kwargs.update(eqdict)
        
        results = numpy.zeros((self.nfun,self.narg),dtype=float)
        for ifun in range(self.nfun):
            fun = self.funs[ifun]
            for iarg in range(self.narg):
                farg = self.fargs[iarg]
                res = fun(*farg,**kwargs)
                results[ifun,iarg] = res
        self.results = results
        self.errs = _geterr(self.refs,results,zerotol=zerotol)
        return None
    
    def printresults(self,chktol=_CHKTOL):
        """Print results outside of a tolerance.
        """
        if self.results is None:
            errmsg = ('Results have not been calculated yet. Use '
                'Tester.run to get function results.')
            raise ValueError(errmsg)
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
                ref = self.refs[ifun,iarg]
                res = self.results[ifun,iarg]
                err = self.errs[ifun,iarg]
                line = fname.rjust(nstr) + '  '
                line += self.argfmt.format(*farg) + '  '
                line += '  '.join(_RESFMT.format(val)
                    for val in (ref,res,err))
                print(line)
        return None


'''
def runcheck(checkdict,printresult=True,tol=_CHKTOL,header=None):
    """Compare function results to reference values.
    
    Run a set of functions over a series of inputs and compare the results to reference values. The functions, inputs, and reference values are all part of the input dictionary checkdict. This function will add the function results and relative errors to this dictionary.
    
    :arg dict checkdict: Information needed to run the check. This dictionary will be modified by this function. It must contain the following:
        
        * :key str modname: Name of the parent module of the functions being checked.
        * :key args: 
    
    Args:
        checkdict (dict): Information needed to run the check. This dictionary
            will be modified by the function. It may contain the following:
            modname (str): Name of the parent module of the functions being
                checked.
            type (str): Either 'fun' to run a single set of inputs over several
                functions or 'der' to run a single function with several inputs.
            args (tuple): Non-keyword arguments to pass to the functions. If the
                'type' is 'der', these are the inputs *not* being varied; if all
                inputs are varied, set 'args' to 'tuple()'.
            ders (tuple): If the 'type' is 'der', these are the inputs being
                varied.
            kwargs (tuple): Keyword arguments to pass to the functions.
            funs (tuple or function): If the 'type' is 'fun', this lists the
                functions to be checked; if the 'type' is 'der', this is just
                the single function to be checked.
            names (tuple or str): Names of the functions being checked.
            geteqvals (function): If provided, this function will be used to
                calculate any auxiliary variables first, so that inversion
                routines are not re-run for every function or set of inputs.
            eqkws (tuple): If a 'geteqvals' function is given, this lists the
                names of the auxiliary variables calculated by the function.
            refs (tuple): Reference values for the functions being checked.
            tol (float): Tolerance to use for printing results. Overrides the
                'tol' argument provided to runcheck.
            header (str): Format to use when printing results for this
                dictionary. Overrides the 'header' argument provided to
                runcheck.
        printresult (boolean, optional): If True (default), any result outside
            the given tolerance is printed, along with the module, function,
            inputs, and reference values. A message is also printed if all
            results are within the given tolerance.
        tol (float, optional): Tolerance to use for printing results (default
            1e-8).
        header (str or None, optional): Format to use when printing results.
            If None (default), a header is constructed from the module and
            function names and the inputs. If a 'header' key is in checkdict,
            it will override any value given here.
    
    Returns None, but adds the following to checkdict:
        results (list): Results returned by the functions.
        stats (list): Relative error between the results and the reference
            values.
    """
    
    # Grab values from the dictionary and initialize results
    chkType = checkdict['type']
    if chkType not in ('fun','der'):
        errmsg = 'The "type" must be either "fun" or "der"'
        raise ValueError(errmsg)
    
    chkArgs = checkdict['args']
    chkKWArgs = checkdict.get('kwargs',dict())
    tol = checkdict.get('tol', tol)
    
    # Determine kwargs to print
    prntKWArgs = dict(chkKWArgs)
    for key in _NOPRINT:
        if key in prntKWArgs.keys():
            del prntKWArgs[key]
    
    # Get equilibrium values if needed
    geteqvals = checkdict.get('geteqvals',None)
    if geteqvals is not None:
        eqkws = checkdict['eqkws']
        eqvals = geteqvals(*chkArgs,**chkKWArgs)
        narg = len(chkArgs)
        chkKWArgs.update({kw: r for (kw,r) in zip(eqkws,eqvals[narg:])})
        chkKWArgs['chkvals'] = False
    
    # Run each function
    results = list()
    stats = list()
    refs = checkdict['refs']
    if chkType == 'fun':
        for (fun,ref) in zip(checkdict['funs'],refs):
            result = fun(*chkArgs,**chkKWArgs)
            stat = _getstat(ref,result)
            results.append(result)
            stats.append(stat)
    else:
        fun = checkdict['funs']
        for (der,ref) in zip(checkdict['ders'],refs):
            result = fun(*(der + chkArgs),**chkKWArgs)
            stat = _getstat(ref,result)
            results.append(result)
            stats.append(stat)
    checkdict['results'] = results
    checkdict['stats'] = stats
    
    # If not printing the results, end here
    if not printresult:
        return None
    
    # Create header if None is given
    header = checkdict.get('header',header)
    if header is None:
        if chkType == 'fun':
            h0 = checkdict['modname']
        else:
            h0 = checkdict['modname'] + '.' + checkdict['names']
        if len(chkArgs) > 0:
            h1 = '(' + ','.join(str(arg) for arg in chkArgs)
        else:
            h1 = str()
        if len(prntKWArgs) > 0:
            if len(chkArgs) > 0:
                h2 = ','
            else:
                h2 = '('
            h3 = ','.join('{0}={1}'.format(*item)
                for item in prntKWArgs.items()) + ')'
        else:
            h2 = str()
            if len(chkArgs) > 0:
                h3 = ')'
            else:
                h3 = str()
        header = h0 + h1 + h2 + h3
    
    # Are all the stats okay?
    if all(abs(stat) <= tol for stat in stats):
        msg = header + ': all functions within {0}'.format(tol)
        print(msg)
        return None
    
    # Print out all the stats outside of the tolerances
    msg = header + ': functions outside of {0}'.format(tol)
    print(msg)
    nval = len(refs)
    for ival in range(nval):
        stat = stats[ival]
        if abs(stat) > tol:
            # Print which function is wrong
            if chkType == 'fun':
                line = '\t' + checkdict['names'][ival]
            else:
                der = checkdict['ders'][ival]
                line = '\t(' + ','.join(str(d) for d in der) + ')'
            print(line)
            
            # Print the reference, result, and relative error
            ref = checkdict['refs'][ival]
            val = checkdict['results'][ival]
            print('\t\tref: {0}'.format(ref))
            print('\t\tres: {0}'.format(val))
            print('\t\terr: {0}'.format(stat))
    
    return None
'''

