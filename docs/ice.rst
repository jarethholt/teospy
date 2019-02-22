.. teospy ice modules documentation

``ice`` modules
===============

The ``ice`` modules provide thermodynamic properties of ice, specifically
hexagonal ice 1 which is most relevant to geophysical applications. The primary
variables are temperature and pressure. The module :mod:`~teospy.ice1` provides
the Gibbs energy function and :mod:`~teospy.ice2` provides all other
thermodynamic properties.


``ice1``
--------

.. automodule:: teospy.ice1
    :members: ice_g


``ice2``
--------

.. automodule:: teospy.ice2
    :members: chempot, cp, density, enthalpy, entropy, expansion,
        helmholtzenergy, internalenergy, kappa_s, kappa_t, lapserate,
        pcoefficient, specificvolume, chkiapws06table6

