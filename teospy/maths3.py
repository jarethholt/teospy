"""Root-finding methods.

This module provides a Newton solver for finding the root of a function.
A root-finding method is required for modules of level 3 and higher to
calculate the values of primary variables at thermodynamic equilibrium
from given secondary variables.

For example, fluid water properties (flu_1) are constructed from the
Helmholtz free energy, a function of temperature and density. The Gibbs
free energy (flu_3a) is a function of temperature and pressure, so the
equation

    pres = dliq^2 * d(free(temp,dliq))/d(dliq)

has to be inverted to find the liquid density dliq.

This module may be expanded in the future to provide more choices for
the inversion method (e.g. Brent, secant). However, this would require
modifying the high-level modules to allow a choice of solver, and Newton
iteration should still be the recommended and default choice. It also
currently uses numpy to solve the linearized equations at each step.
This requirement may also be worked around in the future if there is
demand for it.

Several constants are provided here as well. RYTOL is the recommended
relative tolerance to use to decide when the iteration is accurate
enough. RXTOL is the recommended relative tolerance to use to decide
when the step sizes in the iteration are too small to continue. MAXITER
is the recommended maximum number of iterations to allow before
returning and raising a warning.

:Functions:

* :func:`newton`: Equate two quantities using Newton iteration.

"""

__all__ = ['RTOL','MAXITER','newton']

import warnings
import numpy

# Default values
RYTOL = 1e-8  # Use a relative tolerance of 0.01 ppm
RXTOL = 1e-12  # Limit step sizes to >1 ppt
MAXITER = 100  # Use up to 100 iterations

# Formatting specification for printing warnings
_FLOATFMT = '{:8.4e}'
_NONELEN = 10


### Auxiliary function
def _fmttol(tols,n,stacklevel=2):
    """Format tolerances to fit a length.
    
    Format the given tolerances as a numpy array of the given length.
    The tolerances can be a single float, None, or an iterable of floats
    and Nones.
    
    :arg tols: Given tolerances to format.
    :type tols: float or None or iterable(float or None)
    :arg int n: Length of the output array needed.
    :arg int stacklevel: Controls how many levels deep to raise warnings
        from (default 2).
    :returns: Standardized tolerances of the right shape.
    :rtype: numpy.array(float)
    :raises RuntimeWarning: If tols is an iterable but does not have
        length n. The output gains the desired length either by
        truncating or by appending numpy.inf values.
    """
    # What type of input was given?
    if tols is None:
        tols_np = numpy.ones((n,), dtype=float) * numpy.inf
        return tols_np
    elif isinstance(tols,float):
        tols_np = numpy.ones((n,), dtype=float) * tols
        tols_np[tols_np < 0] = numpy.inf
        return tols_np
    
    # Assume tols is an iterable
    if len(tols) < n:
        warnmsg = ('tols has length {0} < desired length {1}. Additional '
            'components will not be checked.').format(len(tols),n)
        warnings.warn(warnmsg,RuntimeWarning,stacklevel=stacklevel)
        tols_np = numpy.ones((n,), dtype=float) * numpy.inf
        for (it,tol) in enumerate(tols):
            if tol is not None:
                tols_np[it] = tol
    elif len(tols) > n:
        warnmsg = ('tols has length {0} > desired length {1}. The array will '
            'be truncated to match the length.').format(len(tols),n)
        warnings.warn(warnmsg,RuntimeWarning,stacklevel=stacklevel)
        tols_np = numpy.array(tols[:n], dtype=float)
    else:
        tols_np = numpy.array(tols, dtype=float)
    
    # Adjust by removing negative and NaN values (NaN comes from None)
    tols_np[numpy.isnan(tols_np)] = numpy.inf
    tols_np[tols_np < 0] = numpy.inf
    return tols_np


### Main iteration routine
def newton(fun,x0,fargs=None,fkwargs=None,maxiter=MAXITER,rxtol=RXTOL,
    axtol=None,rytol=RYTOL,aytol=None,gamma=0.,stacklevel=2):
    """Equate two functions using Newton iteration.
    
    Equate two functions using Newton iteration. Functionally equivalent
    to finding the root of a function, but allows for checking relative
    as well as absolute tolerance during iteration. This routine's
    primary purpose is calculating primary variables from secondary
    variables, e.g. density from temperature and pressure.
    
    There are several ways to specify the tolerances as described below.
    For each one, different tolerances can be specified for each
    component of a multi-variable function. To bypass checking a
    tolerance, use None or numpy.inf.
    
    :arg fun: A function with the format
        `(lhs, rhs, dlhs, drhs) = fun(\\*(x+fargs),\\*\\*fkwargs)` where
        (lhs,rhs) at the two quantities to be equated; x contains the
        primary variables to be calculated; and (dlhs,drhs) are the
        Jacobians of the two quantities with respect to x.
    :type fun: function
    :arg x0: Initial guess for the values of the primary variables.
    :type x0: float or iterable(float)
    :arg fargs: Values of the secondary variables, which will be held
        constant while solving. If no secondary variables are present,
        use None (default).
    :type fargs: None or float or tuple(float)
    :arg fkwargs: Additional keyword arguments to pass to fun. Use None
        (default) for no keyword arguments.
    :type fkwargs: None or dict
    :arg int maxiter: Maximum number of iterations to allow before
        stopping (default MAXITER).
    :arg rxtol: Relative tolerance between successive values of x for
        halting the iteration (default RXTOL).
    :type rxtol: None or float or iterable(float or None)
    :arg axtol: Absolute tolerance between successive values of x for
        halting the iteration (default None). These values have implied
        units and thus should be different for each component.
    :type axtol: None or float or iterable(float or None)
    :arg rytol: Relative tolerance to use when comparing (lhs,rhs) for
        halting the iteration (default RYTOL).
    :type rytol: None or float or iterable(float or None)
    :arg aytol: Absolute tolerance to use when comparing (lhs,rhs) for
        halting the iteration (default None). These values have implied
        units and thus should be different for each component.
    :type aytol: None or float or iterable(float or None)
    :arg float gamma: Regularization parameter (default 0.). If the
        iteration oscillates without converging, try setting this to a
        positive value less than 1.
    :arg int stacklevel: Controls how many levels deep to raise warnings
        from (default 2).
    :returns: Final estimate of the primary variables x. A float is
        returned for single-variable functions and an array for multi-
        variable functions.
    :rtype: float or array(float)
    :raises RuntimeWarning: If maxiter is non-positive. In this case, no
        iteration is done and the initial value x0 is returned.
    :raises RuntimeWarning: If any of (rxtol,axtol,rytol,aytol) are
        iterables but do not have the same length as x0. The program
        will try to append or truncate the iterables to match.
    :raises RuntimeWarning: If the step sizes for x are smaller than the
        tolerances (rxtol,axtol) before the y values are within
        (rytol,aytol).
    :raises RuntimeWarning: If the maximum number of iterations is
        reached before the error is below the given tolerance.
    
    :Examples:
    
    >>> # Equate w - z^2 = z e^z for w = 0.2
    >>> import math
    >>> def fun(z,w):
    ...     lhs = w - z**2
    ...     rhs = z*math.exp(z)
    ...     dlhs = -2*z  # d(lhs)/dz
    ...     drhs = (z+1)*math.exp(z)  # d(rhs)/dz
    ...     return (lhs,rhs,dlhs,drhs)
    >>> w = 0.2
    >>> z0 = 0.0
    >>> fargs = (w,)
    >>> z1 = newton(fun,z0,fargs)
    >>> print(z1)
    0.15196574180900127
    >>> lhs, rhs, dlhs, drhs = fun(z1,w)
    >>> print(rhs/lhs-1)  # Relative error
    2.220446049250313e-16
    """
    # For no iterations, return initial guess
    if maxiter <= 0:
        warnmsg = 'Maxiter {0} is <= 0; no iteration was done'.format(maxiter)
        warnings.warn(warnmsg,RuntimeWarning,stacklevel=stacklevel)
        x1 = x0
        return x1
    
    # Reformat initial guess, fargs/fkwargs, and tolerances
    x1 = numpy.atleast_1d(x0)
    nx = x1.size
    if fargs is None:
        fargs = tuple()
    if fkwargs is None:
        fkwargs = dict()
    rxtol_np = _fmttol(rxtol,nx)
    axtol_np = _fmttol(axtol,nx)
    rytol_np = _fmttol(rytol,nx)
    aytol_np = _fmttol(aytol,nx)
    
    # Main iteration loop
    for it in range(maxiter):
        # Calculate function and its derivatives
        allargs = x1.tolist() + list(fargs)
        lhs, rhs, dlhs, drhs = fun(*allargs,**fkwargs)
        lhs = numpy.atleast_1d(lhs)
        rhs = numpy.atleast_1d(rhs)
        dlhs = numpy.atleast_2d(dlhs)
        drhs = numpy.atleast_2d(drhs)
        
        # Always do at least one iteration before checking tolerance
        if it == 0:
            dx = numpy.linalg.solve(drhs-dlhs,-(rhs-lhs))
            x1 += dx*(1-gamma)
            continue
        
        # Calculate iteration step
        dx = numpy.linalg.solve(drhs-dlhs,-(rhs-lhs))
        x1 += dx*(1-gamma)
        
        # Are the function values sufficiently close?
        ryerr = numpy.zeros_like(x1)
        cond = (numpy.abs(lhs) >= rytol)
        if numpy.any(cond):
            ryerr[cond] = numpy.abs(rhs[cond]/lhs[cond]-1)
        icond = numpy.invert(cond)
        if numpy.any(icond):
            ryerr[icond] = numpy.abs((rhs-lhs)[icond])
        ayerr = numpy.abs(rhs-lhs)
        if (numpy.all(ryerr<rytol_np) and numpy.all(ayerr<aytol_np)):
            break
        
        # Are the steps between iterations too small to continue?
        rxdiff = numpy.zeros_like(x1)
        cond = (numpy.abs(x1) >= rxtol)
        if numpy.any(cond):
            rxdiff[cond] = numpy.abs(dx[cond]/x1[cond])
        icond = numpy.invert(cond)
        if numpy.any(icond):
            rxdiff[icond] = numpy.abs(dx[icond])
        axdiff = numpy.abs(dx)
        if (numpy.all(rxdiff<rxtol_np) and numpy.all(axdiff<axtol_np)):
            warnmsg = 'Step sizes are smaller than accepted tolerance.\n'
            if numpy.all(rxdiff < rxtol_np):
                rxtol_str = [_FLOATFMT.format(tol) if tol<numpy.inf
                    else 'None'.rjust(_NONELEN) for tol in rxtol_np[:]]
                rxdif_str = [_FLOATFMT.format(dif) for dif in rxdiff[:]]
                warnmsg += ('rxtol:  ' + '  '.join(tol for tol in rxtol_str)
                    + '\n')
                warnmsg += ('rxdif:  ' + '  '.join(dif for dif in rxdif_str)
                    + '\n')
            if numpy.all(axdiff < axtol_np):
                axtol_str = [_FLOATFMT.format(tol) if tol<numpy.inf
                    else 'None'.rjust(_NONELEN) for tol in axtol_np[:]]
                axdif_str = [_FLOATFMT.format(dif) for dif in axdiff[:]]
                warnmsg += ('axtol:  ' + '  '.join(tol for tol in axtol_str)
                    + '\n')
                warnmsg += ('axdif:  ' + '  '.join(dif for dif in axdif_str)
                    + '\n')
            warnmsg = warnmsg[:-1]
            warnings.warn(warnmsg,RuntimeWarning,stacklevel=stacklevel)
            break
    else:
        # Maximum number of iterations reached while above tolerance
        warnmsg = ('Maximum number of iterations {0} reached. This might be '
            'due to oscillations; try setting gamma to a positive value '
            '(e.g. gamma=.1) to stabilize.\n').format(maxiter)
        if numpy.any(ryerr >= rytol_np):
            rytol_str = [_FLOATFMT.format(tol) if tol<numpy.inf
                else 'None'.rjust(_NONELEN) for tol in rytol_np[:]]
            ryerr_str = [_FLOATFMT.format(dif) for dif in ryerr[:]]
            warnmsg += 'rytol:  ' + '  '.join(tol for tol in rytol_str) + '\n'
            warnmsg += 'ryerr:  ' + '  '.join(tol for tol in ryerr_str) + '\n'
        if numpy.any(ayerr >= aytol_np):
            aytol_str = [_FLOATFMT.format(tol) if tol<numpy.inf
                else 'None'.rjust(_NONELEN) for tol in aytol_np[:]]
            ayerr_str = [_FLOATFMT.format(dif) for dif in ayerr[:]]
            warnmsg += 'aytol:  ' + '  '.join(tol for tol in aytol_str) + '\n'
            warnmsg += 'ayerr:  ' + '  '.join(tol for tol in ayerr_str) + '\n'
        warnmsg = warnmsg[:-1]
        warnings.warn(warnmsg,RuntimeWarning,stacklevel=stacklevel)
    
    # Reformat the result to match the input type
    if isinstance(x0,float):
        x1 = float(x1)
    elif isinstance(x0,list):
        x1 = x1.tolist()
    elif isinstance(x0,tuple):
        x1 = tuple(x1.tolist())
    return x1


