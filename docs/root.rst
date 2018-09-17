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

The `Lambert W function <https://en.wikipedia.org/wiki/Lambert_W_function>`_ arises when calculating equilibria if the substances are approximated as having constant heat capacities. For a gas, the idealized Gibbs energy function is

.. math:: g_g(T,p) = c_g (T - T_0 - T \ln(T/T_0)) - L \left( \frac{T}{T_0} - 1 \right) + R_* T \ln(p/p_0)

where :math:`c_g` is the isobaric heat capacity, :math:`L` is a latent heat, :math:`R_*` is the specific gas constant, and :math:`T_0,p_0` specify the reference state. The idealized Gibbs energy Gibbs energy of the liquid state is just the first term, :math:`c_{\ell} (T - T_0 - T \ln(T/T_0))`, ignoring pressure and compressibility contributions. Equating these at a given temperature gives the well-known form for the saturation vapour pressure,

.. math:: \ln(p/p_0) = \frac{L}{R_* T_0} \left( 1 - \frac{T_0}{T} \right) + \frac{c_{\ell} - c_g}{R_*} (1 - T_0/T + \ln(T_0/T)) = (A+B) (1 - x) + B \ln(x)

where :math:`A,B` are non-dimensional parameters related to the latent heat of evaporation and the difference in heat capacities. However, trying to solve this equation for the temperature is trickier. It can be recast into the form

.. math:: -(1 + A/B) e^{-(1 + A/B)} (p/p_0)^{1/B} = -(1 + A/B) x e^{-(1 + A/B) x} = z e^z.

The Lambert W function is the inverse function of :math:`z \exp(z)`, so solving for the temperature here involves calculating this function. In particular, we need the lower branch :math:`W_{-1}` for exclusively negative arguments.

While the Lambert W function is well-known, it cannot be expressed in terms of elementary functions. Since this function is only used for making initial guesses for the complete inversion routines, calculating the function exactly would be unnecessary. Instead, we can use the bounds

.. math:: -1 - \sqrt{2u} - u < W_{-1}(-e^{-u-1}) < -1 - \sqrt{2u} - \frac{2}{3} u.

(This inequality was proved in [LWRef]_. It is available `here <https://arxiv.org/abs/1601.04895>`_; see also the `Wikipedia section <https://en.wikipedia.org/wiki/Lambert_W_function#Asymptotic_expansions>`_.) The `lamb` function in `maths4` returns the value

.. math:: lamb(u,\alpha=1/2) = 1 + \sqrt{2u} + \frac{2 + \alpha}{3} u

which, from the above equation, is an approximation for :math:`(1+A/B)T_0/T`. The parameter :math:`0\leq\alpha\leq1` controls whether the result is closer to the lower bound or the upper bound. Since many idealized phase equilibrium problems can be put in the format, this is a useful function.


.. [LWRef] I. Chatzigeorgiou, "Bounds on the Lambert Function and Their Application to the Outage Analysis of User Cooperation," in IEEE Communications Letters, vol. 17, no. 8, pp. 1505-1508, August 2013. doi: 10.1109/LCOMM.2013.070113.130972.




















