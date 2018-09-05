.. teospy tester module documentation

Tester module
=============

The `tester` module provides a single class `Tester` designed to help automate testing. The `Tester` class is initialized with a list of functions to test, values of the input arguments, and an array of reference values. It also takes several arguments related to printing the results. The `run` method then calculates results from this library and relative errors. The `printresults` method prints the function name, arguments, reference value, result, and relative error wherever the relative error is larger than a given tolerance.

The `Tester` class is used in the various `chk*` functions, such as those in `flu1` and `air2`. The class is meant for internal use, but see those examples if you want to write your own tests.


`Tester`
--------

.. automodule:: tester
    :members: Tester

