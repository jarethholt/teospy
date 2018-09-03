.. teospy ice modules documentation

Ice modules
===========

The `ice` modules provide thermodynamic properties of ice, specifically
hexagonal ice 1 (Ih) which is most relevant to geophysical applications. The
primary variables are temperature and pressure. The only modules dealing with
ice are `ice1` that provides the Gibbs energy function and `ice2` that provides
all other thermodynamic properties.



`ice1`
------

.. automodule:: ice1
    :members: ice_g


`ice2`
------

.. automodule:: ice2
    :members: chempot, cp, density, enthalpy, entropy, expansion,
        helmholtzenergy, internalenergy, kappas, kappat, lapserate,
        pcoefficient, specificvolume


