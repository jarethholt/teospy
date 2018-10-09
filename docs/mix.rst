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


`icevap4`
---------

.. automodule:: icevap4
    :members: eq_tp, temperature, pressure, densityvap, chempot, densityice, enthalpyice, enthalpyvap, entropyice, entropyvap, volumesubl, entropysubl, enthalpysubl


`iceair4a`
----------

.. automodule:: iceair4a
    :members: eq_atpe, temperature, pressure, densityair, densityvap, densityice, enthalpysubl, condensationpressure, frostpoint, massfractionair, sublimationpressure, eq_icl, icl, ict, airffromrh_wmo, rhfromairf_wmo, airffromrh_cct, rhfromairf_cct


`iceair4b`
----------

.. automodule:: iceair4b
    :members: iceair_g, cp, density, enthalpy, entropy, expansion, kappa_t, lapserate, solidfraction, vapourfraction


`iceair4c`
----------

.. automodule:: iceair4c
    :members: eq_wpte, iceair_h, cp, density, kappa_s, lapserate, temperature, eq_pot, potdensity, potenthalpy, pottemp


`liqair4a`
----------

.. automodule:: liqair4a
    :members: eq_atpe, temperature, pressure, densityair, densityvap, densityliq, enthalpyevap, entropy, condensationpressure, dewpoint, massfractionair, eq_icl, icl, ict, airffromrh_wmo, rhfromairf_wmo, airffromrh_cct, rhfromairf_cct


`liqair4b`
----------

.. automodule:: liqair4b
    :members: liqair_g, cp, density, enthalpy, entropy, expansion, kappa_t, lapserate, liquidfraction, vapourfraction


`liqair4c`
----------

.. automodule:: liqair4c
    :members: eq_wpte, liqair_h, cp, density, kappa_s, lapserate, temperature, eq_pot, potdensity, potenthalpy, pottemp


`liqiceair4`
------------

.. automodule:: liqiceair4
    :members: eq_atp, eq_wefli, eq_all, pressure, temperature, density, dryairfraction, enthalpy, entropy, liquidfraction, solidfraction, vapourfraction, iml, ifl


`sealiq4`
---------

.. automodule:: sealiq4
    :members: eq_stp, osmoticpressure


`seavap4`
---------

.. automodule:: seavap4
    :members: eq_stp, pressure, salinity, temperature, densitysea, densityvap, enthalpysea, enthalpyvap, entropysea, entropyvap, enthalpyevap, volumeevap, boilingtemperature, brinesalinity, vapourpressure, eq_seavap, seavap_g, cp, density, enthalpy, entropy, expansion, kappa_t, brinefraction


`seaice4`
---------

.. automodule:: seaice4
    :members: eq_stp, densityice, densitysea, enthalpyice, enthalpysea, entropyice, entropysea, pressure, temperature, salinity, enthalpymelt, volumemelt, brinesalinity, meltingpressure, freezingtemperature, dtfdp, dtfds, eq_seaice, seaice_g, brinefraction, cp, density, enthalpy, entropy, expansion, kappa_t



























