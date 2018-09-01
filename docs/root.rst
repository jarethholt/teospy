.. teospy root module documentation

Root modules
============

There are several root-level modules which are not directly related to calculating the thermodynamic properties of sea/ice/air. `constants0` contains physical constants and functions for checking whether inputs are within bounds. `convert0` provides functions for converting between mole and mass fractions and from practical salinity to absolute salinity. `maths3` contains the root-finding method used by the modules of level 3 and higher.


`constants0`
-------------

.. automodule:: constants0
    :members: chkdrybnds, chkflubnds, chkhumbnds, chkicebnds, chksalbnds


`convert0`
-----------

.. automodule:: convert0
    :members: air_molarmass, air_molfractionvap, air_molfractiondry, air_massfractionvap, air_massfractiondry, sal_molality, gsw_safromsp, gsw_spfromsa, sal_asalfrompsal, sal_psalfromasal


`maths3`
---------

.. automodule:: maths3
    :members: newton



























