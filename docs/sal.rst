.. teospy sal modules documentation

``sal`` modules
===============

The ``sal`` modules provide thermodynamic properties of salt in seawater. The
primary variables are salinity, temperature, and pressure. These properties are
combined with liquid water to form seawater in the ``sea`` modules. The
contribution of salt to the Gibbs free energy is formulated as an expansion in
terms of salinity:

.. math:: g(S,T,p) = g_1(T,p) \cdot S \ln(S) + \sum_{i=2}^N g_i(T,p) S^{i/2}

where :math:`S` is the absolute salinity and :math:`N=7` is a constant in
``sal1`` named ``NSALTERMS``.

All functions derived from these modules include the optional boolean keyword
``useext``. When set to True, the Gibbs energy function uses extended
coefficients modified for high temperature and high salinity applications. The
default (False) is recommended for most applications.

The basic contents of the modules are as follows:

* :mod:`~teospy.sal1` provides a function that calculates each of the
  coefficients :math:`g_i` in the salinity expansion.
* :mod:`~teospy.sal2` provides thermodynamic properties of seawater that depend
  only on the salinity, such as dilution and osmotic coefficients.
* :mod:`~teospy.calsalcoeffs` is a script that is used to calculate the
  coefficients of the above formula from the original coefficients, which use
  the square root of salinity. It is provided only to document how the
  coefficients were generated and is not intended to use unless the original
  coefficients are revised by IAPWS.


``sal1``
--------

.. automodule:: teospy.sal1
    :members: sal_g_term


``sal2``
--------

.. automodule:: teospy.sal2
    :members: sal_g, actcoeff, activityw, actpotential, chemcoeff, dilution,
        liqpot, osmcoeff, salpot, saltenthalpy, saltentropy, saltvolume,
        mixenthalpy, mixentropy, mixvolume, eq_liqpot, eq_enthalpy,
        eq_entropy


``calsalcoeffs``
----------------

.. automodule:: teospy.calsalcoeffs

