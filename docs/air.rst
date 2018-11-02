.. teospy air modules documentation

``air`` modules
===============

The ``air`` modules provide thermodynamic properties of both dry and humid air.
(It does not provide properties for saturated air; for those, see the level 4
modules.) For pure dry air, the primary variables are the temperature and dry
air density. For humid air, the primary variables are the dry air mass fraction,
temperature, and humid air density. For converting between mass fraction and
relative humidity, see :mod:`liqair4a` and :mod:`iceair4a`.

The basic contents of the modules are as follows:

* :mod:`air1` provides the Helmholtz free energy of pure dry air and the virial
  coefficients between dry air and water vapour.
* :mod:`air2` provides the Helmholtz free energy of humid air and thermodynamic
  properties derived therefrom.
* :mod:`air3a` and :mod:`air3b` provide properties calculated from dry air mass
  fraction, temperature, and pressure.
* :mod:`air3c` provides properties calculated from dry air mass fraction,
  entropy, and pressure. It also provides dry adiabatic ascent functions (e.g.
  potential temperature).
* :mod:`air5` provides a single function for calculating the adiabatic lapse
  rate in the more conventional units of degrees C per 100 m.


``air1``
--------

.. automodule:: air1
    :members: dry_f, air_baw, air_caaw, air_caww


``air2``
--------

.. automodule:: air2
    :members: air_f, cp, cv, enthalpy, entropy, expansion, gibbsenergy,
        internalenergy, kappa_s, kappa_t, lapserate, pressure, soundspeed,
        vappot, eq_entropy, eq_pressure, eq_vappot, chkiapws10table,
        chkiapws10table13, chkiapws10table14, chkiapws10table15


``air3a``
---------

.. automodule:: air3a
    :members: eq_atp, air_g


``air3b``
---------

.. automodule:: air3b
    :members: compressibility, compressibility_lemmon, contraction, cp, cv,
        density, enthalpy, entropy, expansion, gibbsenergy, internalenergy,
        kappa_s, kappa_t, lapserate, soundspeed, vappot, chklemmon2000


``air3c``
---------

.. automodule:: air3c
    :members: eq_aep, air_h, eq_pot, pottemp, potdensity, potenthalpy


``air5``
--------

.. automodule:: air5
    :members: lapserate_c100m

