.. teospy mixture modules documentation

Mixture modules
===============

The level 4 modules involve mixtures of multiple phases and substances, such
as humid air with liquid water or seawater with ice. They provide functions
for calculating vapour pressure, enthalpy of evaporation/melting, etc. Because
they involve equilibrium between different components, all of the mixture
modules require the root-finder :func:`~teospy.maths3.newton`. The modules may
also rely on functions from level 3 modules, so the mixture modules are all
level 4. (There is one level 5 mixture module, :mod:`~teospy.iceflu5`, that
contains approximate fits for the melting and sublimation curves of ice.)

The basic contents of the modules are as follows:

* :mod:`~teospy.liqvap4` is for equilibrium between the liquid and vapour
  phases of **pure** liquid water.
* :mod:`~teospy.iceliq4` is for equilibrium between ice and pure liquid water.
* :mod:`~teospy.icevap4` is for equilibrium between ice and pure water vapour.
* :mod:`~teospy.iceair4a` is for equilibrium between ice and water vapour in
  humid air (ice-saturated or "icy" air). Any two of the dry air mass
  fraction, temperature, and pressure of icy air are used to calculate the
  third. It also includes functions to convert between relative humidity and
  dry fraction.
* :mod:`~teospy.iceair4b` provides the Gibbs energy and related thermodynamic
  properties for a parcel of icy air. The total water content of a parcel is
  conserved, but the water can convert from ice to vapour as temperature and
  pressure change.
* :mod:`~teospy.iceair4c` provides the enthalpy and related properties for a
  parcel of icy air. This module also provides the potential temperature of
  icy air and the moist adiabatic lapse rate.
* :mod:`~teospy.liqair4a`, :mod:`~teospy.liqair4b`, and
  :mod:`~teospy.liqair4c` are the same as the corresponding ``iceair``
  modules, but for humid air and liquid water.
* :mod:`~teospy.liqiceair4` is for equilibrium of ice, liquid water, and water
  vapour in humid air. Here, only one of the dry fraction, temperature, and
  pressure can be provided, and the other two are calculated.
* :mod:`~teospy.sealiq4` is for seawater and pure liquid water in equilibrium,
  for example across a semipermeable membrane.
* :mod:`~teospy.seavap4` is for seawater and pure water vapour in equilibrium,
  and primarily provides the vapour pressure over seawater.
* :mod:`~teospy.seaice4` is for seawater and ice in equilibrium, providing
  both the freezing temperature and melting pressure of ice over seawater.
* :mod:`~teospy.seaicevap4` is for seawater, ice, and pure water vapour in
  equilibrium. Like :mod:`~teospy.liqiceair4`, only one of the salinity,
  temperature, and pressure can be provided.
* :mod:`~teospy.seaair4` is for seawater and humid air in equilibrium, and
  describes the properties of saturated air over seawater. It also provides
  :func:`~teospy.seaair4.chempotevap` which calculates the Onsager force, a
  measure of the disequilibrium between the air and seawater.
* :mod:`~teospy.iceflu5` provides fits to the equilibrium curves between ice
  and either liquid water or water vapour, i.e. the melting pressure, melting
  temperature, and saturation vapour pressure.


``liqvap4``
-----------

.. automodule:: teospy.liqvap4
    :members: eq_tp, temperature, pressure, densityvap, densityliq, chempot,
        enthalpyevap, enthalpyliq, enthalpyvap, entropyevap, entropyliq,
        entropyvap, volumeevap, chkiapws95table8


``iceliq4``
-----------

.. automodule:: teospy.iceliq4
    :members: eq_tp, temperature, pressure, densityliq, chempot, densityice,
        enthalpyice, enthalpyliq, enthalpymelt, entropyice, entropyliq,
        entropymelt, volumemelt


``icevap4``
-----------

.. automodule:: teospy.icevap4
    :members: eq_tp, temperature, pressure, densityvap, chempot, densityice,
        enthalpyice, enthalpyvap, entropyice, entropyvap, volumesubl,
        entropysubl, enthalpysubl


``iceair4a``
------------

.. automodule:: teospy.iceair4a
    :members: eq_atpe, temperature, pressure, densityair, densityvap,
        densityice, enthalpysubl, condensationpressure, frostpoint,
        massfractionair, sublimationpressure, eq_icl, icl, ict,
        airffromrh_wmo, rhfromairf_wmo, airffromrh_cct, rhfromairf_cct


``iceair4b``
------------

.. automodule:: teospy.iceair4b
    :members: iceair_g, cp, density, enthalpy, entropy, expansion, kappa_t,
        lapserate, solidfraction, vapourfraction


``iceair4c``
------------

.. automodule:: teospy.iceair4c
    :members: eq_wpte, iceair_h, cp, density, kappa_s, lapserate, temperature,
        eq_pot, potdensity, potenthalpy, pottemp


``liqair4a``
------------

.. automodule:: teospy.liqair4a
    :members: eq_atpe, temperature, pressure, densityair, densityvap,
        densityliq, enthalpyevap, entropy, condensationpressure, dewpoint,
        massfractionair, eq_icl, icl, ict, airffromrh_wmo, rhfromairf_wmo,
        airffromrh_cct, rhfromairf_cct


``liqair4b``
------------

.. automodule:: teospy.liqair4b
    :members: liqair_g, cp, density, enthalpy, entropy, expansion, kappa_t,
        lapserate, liquidfraction, vapourfraction


``liqair4c``
------------

.. automodule:: teospy.liqair4c
    :members: eq_wpte, liqair_h, cp, density, kappa_s, lapserate, temperature,
        eq_pot, potdensity, potenthalpy, pottemp


``liqiceair4``
--------------

.. automodule:: teospy.liqiceair4
    :members: eq_atp, eq_wefli, eq_all, airfraction, pressure, temperature,
        density, dryairfraction, enthalpy, entropy, liquidfraction,
        solidfraction, vapourfraction, iml, ifl


``sealiq4``
-----------

.. automodule:: teospy.sealiq4
    :members: eq_stp, osmoticpressure


``seavap4``
-----------

.. automodule:: teospy.seavap4
    :members: eq_stp, pressure, salinity, temperature, densitysea, densityvap,
        enthalpysea, enthalpyvap, entropysea, entropyvap, enthalpyevap,
        volumeevap, boilingtemperature, brinesalinity, vapourpressure,
        eq_seavap, seavap_g, cp, density, enthalpy, entropy, expansion,
        kappa_t, brinefraction


``seaice4``
-----------

.. automodule:: teospy.seaice4
    :members: eq_stp, densityice, densitysea, enthalpyice, enthalpysea,
        entropyice, entropysea, pressure, temperature, salinity, enthalpymelt,
        volumemelt, brinesalinity, meltingpressure, freezingtemperature, dtfdp,
        dtfds, eq_seaice, seaice_g, brinefraction, cp, density, enthalpy,
        entropy, expansion, kappa_t


``seaicevap4``
--------------

.. automodule:: teospy.seaicevap4
    :members: eq_stp, densityvap, pressure, salinity, temperature


``seaair4``
-----------

.. automodule:: teospy.seaair4
    :members: eq_satp, densityair, densityvap, entropyair, enthalpyevap,
        massfractionair, vapourpressure, condensetemp, chempotevap


``iceflu5``
-----------

.. automodule:: teospy.iceflu5
    :members: liqpressure, liqtemperature, vappressure

