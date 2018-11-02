"""Auxiliary math module for approximation functions.

This module provides functions for approximating the Lambert `W`
function in a few different formats. The Lambert function arises when
calculating the equilibrium between two phases that have constant heat
capacities. The `W` function is the inverse of :math:`f(z) = z e^z`. The
form most often used here is aimed at the negative branch :math:`W_{-1}`
for negative values, though occassionally the principal branch
:math:`W_0` is needed. See :download:`/latex/approximations.pdf` for
more details.

:Functions:

* :func:`lamb1`
* :func:`lamb2`
* :func:`lamb0`

"""

__all__ = ['lamb1','lamb2','lamb0']
import numpy

_C_SP = 0.396166676603


def lamb1(u,alpha=.5):
    """Approximate the Lambert W function.
    
    Approximate the Lambert W function from its upper and lower bounds.
    The parameter alpha (between 0 and 1) determines how close the
    approximation is to the lower bound instead of the upper bound.
    
    :arg float u: Modified argument of the function.
    :arg float alpha: Bound parameter (default 0.5).
    :returns: (-z)-value of the Lambert function.
    :raises ValueError: If u is negative.
    :raises ValueError: If alpha is not between 0 and 1.
    """
    if u < 0:
        errmsg = 'Argument u must be positive'
        raise ValueError(errmsg)
    if alpha < 0 or alpha > 1:
        errmsg = 'Parameter alpha must be between 0 and 1'
        raise ValueError(errmsg)
    beta = (2 + alpha)/3
    negz = 1 + (2*u)**.5 + beta*u
    return negz

def lamb2(v,r):
    """Approximate the Lambert W function.
    
    Approximate the Lambert W function from its upper and lower bounds.
    The result replicates the triple point of the (idealized) system
    exactly, because lamb2(0,r) = 1.
    
    :arg float v: Modified argument of the function. Must be positive.
    :arg float r: Latent heat/heat capacity ratio. Must be positive.
    :returns: x-value in the Lambert function.
    """
    u0 = r - numpy.log(1+r)
    negz = 1 + (2*(u0-v))**.5 + (r - (2*u0)**.5)/u0*(u0-v)
    x = negz / (1+r)
    return x

def lamb0(w,gam=_C_SP):
    """Approximate Lambert function principal branch.
    
    Approximate the principal branch of the Lambert W function using a
    power law. The default exponent is the result of fitting to the
    range -1/e < w < 0.
    
    :arg float w: Argument of the Lambert function. Should be negative.
    :arg float gam: Exponent to use in the power law. Default (_C_SP)
        comes from fitting the power law to the negative range.
    :returns: (-z)-value of the Lambert function.
    """
    negz = 1 - (1 + w*numpy.exp(1))**gam
    return negz

