.. teospy root module documentation

Root modules
============

There are several root-level modules which are not directly related to calculating the thermodynamic properties of sea/ice/air. `constants0` contains physical constants and functions for checking whether inputs are within bounds. `convert0` provides functions for converting between mole and mass fractions and from practical salinity to absolute salinity. `maths3` contains the root-finding method used by the modules of level 3 and higher. `maths4` contains a convenient function commonly used for approximating phase equilibria in the level 4 (mixture) modules.


`constants0`
------------

.. automodule:: constants0
    :members: chkdrybnds, chkflubnds, chkhumbnds, chkicebnds, chksalbnds


`convert0`
----------

.. automodule:: convert0
    :members: air_molarmass, air_molfractionvap, air_molfractiondry, air_massfractionvap, air_massfractiondry, sal_molality, gsw_safromsp, gsw_spfromsa, sal_asalfrompsal, sal_psalfromasal


`maths3`
--------

.. automodule:: maths3
    :members: newton


`maths4`
--------

The `Lambert W function <https://en.wikipedia.org/wiki/Lambert_W_function>`_ arises when calculating equilibria if the substances are approximated as having constant heat capacities. The Lambert W function is the inverse function of :math:`z \exp(z)`. In particular, we need the lower branch :math:`W_{-1}` for exclusively negative arguments.

While the Lambert W function is well-known, it cannot be expressed in terms of elementary functions. Since this function is only used for making initial guesses for the complete inversion routines, calculating the function exactly would be unnecessary. Instead, we use the bounds

.. math:: -1 - \sqrt{2u} - u < W_{-1}(-e^{-u-1}) < -1 - \sqrt{2u} - \frac{2}{3} u.

(This inequality was proved in [LWRef]_. It is available `here <https://arxiv.org/abs/1601.04895>`_; see also the `Wikipedia section <https://en.wikipedia.org/wiki/Lambert_W_function#Asymptotic_expansions>`_.) The `lamb1` function in `maths4` returns the value

.. math:: \text{lamb1}(u,\alpha=1/2) = 1 + \sqrt{2u} + \frac{2 + \alpha}{3} u.

The parameter :math:`0\leq\alpha\leq1` controls whether the result is closer to the lower bound or the upper bound. The `lamb2` function is a variant that returns

.. math:: \text{lamb2}(v,r) = \frac{1}{1+r} \left[ 1 + \sqrt{2(u_0-v)} + \beta (u_0-v) \right]

.. math:: u_0(r) = r - \ln(1+r), \qquad \beta(r) = (r - \sqrt{2u_0})/u_0.

This function is an approximation for the value of :math:`x` in the equation :math:`y = a (1-x) + b \ln(x)`, with :math:`r=a/b-1` and :math:`v=y/b`. The advantage here is that the solution correctly returns the known value :math:`x=1` when :math:`y=0`.


.. [LWRef] I. Chatzigeorgiou, "Bounds on the Lambert Function and Their Application to the Outage Analysis of User Cooperation," in IEEE Communications Letters, vol. 17, no. 8, pp. 1505-1508, August 2013. doi: 10.1109/LCOMM.2013.070113.130972.




















