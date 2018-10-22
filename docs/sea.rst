.. teospy sea modules documentation

Sea modules
===========

The `sea` modules provide thermodynamic properties of seawater by combining the Gibbs free energy of salt in seawater (`sal2`) with the Gibbs free energy of liquid water (`flu3a`). The primary variables are salinity, temperature, and pressure, but even the most basic version requires solving for the liquid water density. Thus all `sea` modules are level 3 and above.

The basic contents of the modules are as follows:

* :mod:`sea3a` provides the Gibbs free energy of seawater and related properties, using salinity, temperature, and pressure as primary variables.
* :mod:`sea3b` provides the enthalpy of seawater from salinity, entropy, and pressure. It also provides adiabatic ascent functions (e.g. potential temperature) and the related expansion and contraction coefficients defined with respect to potential temperature and potential enthalpy.
* :mod:`sea3c` provides the same functions as `sea3b`, but based on the seawater entropy with salinity, enthalpy, and pressure as primary variables.
* :mod:`sea3d` provides only one function, for calculating salinity from in-situ temperature, pressure, and density.
* :mod:`sea5` provides functions for converting to and from conservative temperature, a scaled potential enthalpy that may be more closely conserved under oceanic advection than potential temperature or entropy. It also provides the cabbeling and thermobaric coefficients, quantities related to dianeutral mixing.
* :mod:`gsw5` provides many of the same functions as `sea5`, but using the polynomial version of the liquid water Gibbs energy from :mod:`liq5_f03`. This removes the need to invert liquid water density from pressure, so many quantities can be directly calculated. Since both the liquid water and salt free energies are polynomials, they support high-order derivatives, and so the cabbeling and thermobaric coefficients can be calculated directly instead of using finite differences.


`sea3a`
-------

.. automodule:: sea3a
    :members: sea_g, liqpot, salpot, contraction_t, cp, density, enthalpy, helmholtzenergy, entropy, expansion_t, gibbsenergy, internalenergy, kappa_s, kappa_t, lapserate, osmcoeff, soundspeed, temp_maxdensity, chkiapws08table8a, chkiapws08table8b, chkiapws08table8c


`sea3b`
-------

.. automodule:: sea3b
    :members: eq_sep, sea_h, temperature, contraction_t, expansion_t, eq_pot, pottemp, potdensity, potenthalpy, contraction_h, contraction_theta, expansion_h, expansion_theta



`sea3c`
-------

.. automodule:: sea3c
    :members: eq_shp, eq_pot, temperature, density, entropy, contraction_t, expansion_t, pottemp, potdensity, contraction_h, contraction_theta, expansion_h, expansion_theta


`sea3d`
-------

.. automodule:: sea3d
    :members: eq_tpd, salinity


`sea5`
------

.. automodule:: sea5
    :members: tconfromtpot, tpotfromtcon, expansion_tcon, expansion_tpot, expansion_t, contraction_tcon, contraction_tpot, contraction_t, cabb_tcon, cabb_tpot, thrmb_tcon, thrmb_tpot


`gsw5`
------

.. automodule:: gsw5
    :members: asalfrompsal, psalfromasal, gsw_g, alpha_t, beta_t, cp, density, enthalpy, entropy, kappa, kappa_t, specvol, svel, pottemp, potdens, alpha_tpot, beta_tpot, cabb_tpot, cabb_tpot_alt, thrmb_tpot, thrmb_tpot_alt, tconfromtpot, tpotfromtcon, alpha_tcon, beta_tcon, cabb_tcon, cabb_tcon_alt, thrmb_tcon, thrmb_tcon_alt



























