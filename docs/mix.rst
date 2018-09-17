.. teospy mixture modules documentation

Mixture modules
===============

The level 4 modules involve mixtures of multiple phases and substances, such as humid air with liquid water or seawater with ice. As such, they provide functions for calculating vapour pressure, enthalpy of evaporation/melting, etc. Because they involve equilibrium between different components, all of the mixture modules require the root-finder :func:`maths3.newton`. The modules may also rely on functions from level 3 modules to calculate quantities more easily. Thus the mixture modules are all level 4.

The basic contents of the modules are as follows:

* :mod:`liqvap4` deals with the equilibrium between the liquid and vapour phases of **pure** liquid water.


`liqvap4`
---------

.. automodule:: liqvap4
    :members: eq_tp, temperature, pressure, densityvap, densityliq, chempot, enthalpyevap, enthalpyliq, enthalpyvap, entropyevap, entropyliq, entropyvap, volumeevap, chkiapws95table8


`iceliq4`
---------

.. automodule:: iceliq4
    :members: eq_tp, temperature, pressure, densityliq, chempot, densityice, enthalpyice, enthalpyliq, enthalpymelt, entropyice, entropyliq, entropymelt, volumemelt





























