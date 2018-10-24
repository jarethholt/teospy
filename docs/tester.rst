.. teospy tester module documentation

Test modules
============

The `tester` module provides a single class :class:`Tester` designed to help
automate testing. This class is initialized with a list of functions to test,
values of the input arguments, and an array of reference values. It also takes
several arguments related to printing the results. The `run` method then
calculates results from this library and relative errors. The `printresults`
method prints the function name, arguments, reference value, result, and
relative error wherever the relative error is larger than a given tolerance.

Several tests have been compiled into modules such as :mod:`testair` to readily
separate the tests by substance. These modules correspond to modules such as
`Values_Air` in the Fortran version of TEOS. The `Tester` class is also used in
the various `chk` functions, for example :func:`~flu1.chkiapws95table6` and
:func:`~air2.chkiapws10table`.


`tester`
--------

.. automodule:: tester
    :members: Tester

.. autoclass:: Tester
    .. method:: run, printresults


`testair`
---------

.. automodule:: testair


`testcnv`
---------

.. automodule:: testcnv


`testflu`
---------

.. automodule:: testflu


`testgsw`
---------

.. automodule:: testgsw


`testice`
---------

.. automodule:: testice


`testmix`
---------

.. automodule:: testmix


`testsal`
---------

.. automodule:: testsal


`testsea`
---------

.. automodule:: testsea

