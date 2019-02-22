.. teospy root module documentation

Root modules
============

There are several root-level modules which are not directly related to
calculating the thermodynamic properties of sea/ice/air.
:mod:`~teospy.constants0` contains physical constants and functions for checking
whether inputs are within bounds. :mod:`~teospy.convert0` provides functions for
converting between mole and mass fractions and from practical salinity to
absolute salinity. :mod:`~teospy.maths3` contains the root-finding method used
by the modules of level 3 and higher. :mod:`~teospy.maths4` contains a few
functions that are commonly used for approximating phase equilibria in the level
4 (mixture) modules. For more details on these approximation functions, see the
notes in :download:`/latex/approximations.pdf`.



``constants0``
--------------

.. automodule:: teospy.constants0
    :members: chkdrybnds, chkflubnds, chkhumbnds, chkicebnds, chksalbnds


``convert0``
------------

.. automodule:: teospy.convert0
    :members: air_molarmass, air_molfractionvap, air_molfractiondry,
        air_massfractionvap, air_massfractiondry, sal_molality, gsw_safromsp,
        gsw_spfromsa, sal_asalfrompsal, sal_psalfromasal


``maths3``
----------

.. automodule:: teospy.maths3
    :members: newton


``maths4``
----------

.. automodule:: teospy.maths4
    :members: lamb1, lamb2, lamb0


``convert5``
------------

.. automodule:: teospy.convert5
    :members: swpres, swdepth, cnvpressure, cnvtemperature,
        sal78fromcnd, cndfromsal78, cnvsalinity

