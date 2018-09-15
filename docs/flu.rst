.. teospy flu modules documentation

Flu modules
===========

The `flu` modules provide thermodynamic properties of fluid water, both liquid
water and water vapour. The primary variables are the temperature and fluid
density.

The basic contents of the modules are as follows:

* :mod:`flu1` provides the Helmholtz free energy of fluid water.
* :mod:`flu2` provides thermodynamic properties that can be derived directly
  from the Helmholtz free energy in the primary variables.
* :mod:`flu3a` provides the Gibbs free energy of liquid water and water vapour
  separately, with temperature and pressure as variables.
* :mod:`flu3b` provides other properties of liquid water and water vapour using
  temperature and pressure as variables.
* The level 5 modules provide formulations of fluid water properties based on
  Gibbs energy functions for each phase. :mod:`flu5_if97` provides Gibbs energy
  functions for both phases, whereas :mod:`liq5_f03` provides functions for the
  liquid phase only.


`flu1`
------

.. automodule:: flu1
    :members: flu_f, chkiapws95table6


`flu2`
------

.. automodule:: flu2
    :members: cp, cv, enthalpy, entropy, expansion, gibbsenergy,
        internalenergy, kappa_s, kappa_t, lapserate, pressure, soundspeed,
        eq_chempot, eq_enthalpy, eq_entropy, eq_pressure, chkiapws95table7


`flu3a`
-------

.. automodule:: flu3a
    :members: eq_tp_liq, eq_tp_vap, liq_g, vap_g


`flu3b`
-------

.. automodule:: flu3b
    :members: liq_cp, liq_cv, liq_density, liq_enthalpy, liq_entropy,
        liq_expansion, liq_gibbsenergy, liq_helmholtzenergy,
        liq_internalenergy, liq_kappa_s, liq_kappa_t, liq_lapserate,
        liq_soundspeed, vap_cp, vap_cv, vap_density, vap_enthalpy,
        vap_entropy, vap_expansion, vap_gibbsenergy, vap_helmholtzenergy,
        vap_internalenergy, vap_kappa_s, vap_kappa_t, vap_lapserate,
        vap_soundspeed


`flu5_if97`
-----------

.. automodule:: flu5_if97
    :members: liq_g, liq_cp, liq_density, liq_enthalpy, liq_entropy,
        liq_internalenergy, liq_soundspeed, liq_volume, vap_g, vap_cp,
        vap_density, vap_enthalpy, vap_entropy, vap_internalenergy,
        vap_soundspeed, vap_volume, chkiapws97table, chkiapws97table5,
        chkiapws97table15


`liq5_f03`
----------

.. automodule:: liq5_f03
    :members: liq_g, cp, density, expansion, kappa_t, soundspeed, enthalpy,
        entropy, helmholtzenergy, internalenergy, chkiapws09table6

