.. teospy sal modules documentation

Sal modules
===========

The `sal` modules provide thermodynamic properties of salt in seawater. The primary variables are salinity, temperature, and pressure. These properties are combined with liquid water to form seawater in the `sea3\*` modules. The contribution of salt to the Gibbs free energy is formulated as an expansion in terms of salinity:

    g(S,T,p) = g1(T,p) S ln(S) + sum_{i=2}^{NSALTERMS} gi(T,p) S^(i/2)

where S is the absolute salinity and NSALTERMS=7.

All functions derived from these modules include the optional boolean keyword `useext`. When set to True, the Gibbs energy function uses extended coefficients modified for high temperature and high salinity applications. The default (False) is recommended for most applications.

The basic contents of the modules are as follows:

* `sal1` provides a function that calculates each of the coefficients gi in the salinity expansion.
* `sal2` provides thermodynamic properties of seawater that depend only on the salinity, such as dilution and osmotic coefficients.


`sal1`
------

.. automodule:: sal1
    :members: sal_g_term


`sal2`
------

.. automodule:: sal2
    :members: sal_g, actcoeff, activityw, actpotential, chemcoeff, dilution,
        liqpot, osmcoeff, salpot, saltenthalpy, saltentropy, saltvolume,
        mixenthalpy, mixentropy, mixvolume, eq_liqpot, eq_enthalpy,
        eq_entropy


