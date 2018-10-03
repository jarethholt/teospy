"""Auxiliary math module for approximation functions.

This module provides a single function `lamb` for calculating an
approximation of the Lambert W function in a certain form. The Lambert
function arises when calculating the equilibrium between two phases that
have constant heat capacities.
"""

__all__ = ['calalpha','lamb1','lamb2']
import numpy


def calalpha(r):
    """Calculate the alpha value for a system.
    
    Calculate the alpha value to use in the Lambert approximation that
    replicates the triple point of the system exactly. The argument here
    is the ratio of the latent heat to the heat capacity difference.
    
    :arg float r: Latent heat/heat capacity ratio.
    :returns: Value of alpha parameter.
    """
    u0 = r - numpy.log(1+r)
    alpha = 3*(r - (2*u0)**.5)/u0 - 2
    return alpha

def lamb1(u,alpha=.5):
    """Approximate the Lambert W function.
    
    Approximate the Lambert W function from its upper and lower bounds.
    The parameter alpha (between 0 and 1) determines how close the
    approximation is to the lower bound instead of the upper bound.
    
    :arg float u: Modified argument of the function.
    :arg float alpha: Bound parameter (default 0.5).
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
    
    Approximate the Lambert W function from its upper and lower bounds. The result replicates the triple point of the (idealized) system exactly. The first argument is proportional to the pressure. The second argument is the ratio of latent heat to the heat capacity difference.
    """
    u0 = r - numpy.log(1+r)
    negz = 1 + (2*(u0-v))**.5 + (r - (2*u0)**.5)/u0*(u0-v)
    x = negz / (1+r)
    return x

