.. teospy air modules documentation

Air modules
===========

The `air` modules provide thermodynamic properties of both dry and humid air. (It does not provide properties for saturated air; for those, see the level 4 modules.) For pure dry air, the primary variables are the temperature and dry air density. For humid air, the primary variables are the dry air mass fraction, temperature, and humid air density. For converting between mass fraction and relative humidity, see `liqair4a` and `iceair4a`.

The basic contents of the modules are as follows:

* `air1` provides the Helmholtz free energy of pure dry air and the virial coefficients between dry air and water vapour.
* `air2` provides the Helmholtz free energy of humid air and thermodynamic properties derived therefrom.
* `air3a` and `air3b` provide properties calculated from dry air mass fraction, temperature, and pressure.
* `air3c` provides properties calculated from dry air mass fraction, entropy, and pressure.
* `air5` provides a single function for calculating the moist adiabatic lapse rate in the more conventional units of degrees C per 100 m.


`air1`
------

.. automodule:: air1
    :members: dry_f, air_baw, air_caaw, air_caww


`air2`
------

.. automodule:: air2
    :members: air_f, cp, cv, enthalpy, entropy, expansion, gibbsenergy, internalenergy, kappas, kappat, lapserate, pressure, soundspeed, vappot, eq_entropy, eq_pressure, eq_vappot































