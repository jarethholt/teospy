.. teospy tester module documentation

Testing modules
===============

The :mod:`tester` module provides a single class :class:`~tester.Tester`
designed to help automate testing. This class is initialized with a list of
functions to test, values of the input arguments, and an array of reference
values. It also takes several arguments related to printing the results. The
:func:`~tester.Tester.run` method then calculates results from this library and
relative errors. The :func:`~tester.Tester.printresults` method prints the
function name, arguments, reference value, result, and relative error whenever
the relative error is larger than a given tolerance.

Several tests have been compiled into modules such as :mod:`testair` to readily
separate the tests by substance. These modules correspond to modules such as
``Values_Air`` in the Fortran version of TEOS. The :class:`~tester.Tester` class
is also used in the various ``chk`` functions, for example
:func:`flu1.chkiapws95table6` and :func:`air2.chkiapws10table`. The
corresponding modules can be called from the command line to run these checks.
The modules with such ``chk`` functions are:

.. hlist::
    :columns: 4
    
    * :mod:`air2`
    * :mod:`air3b`
    * :mod:`flu1`
    * :mod:`flu2`
    * :mod:`flu5_if97`
    * :mod:`ice2`
    * :mod:`liq5_f03`
    * :mod:`liqvap4`
    * :mod:`sea3a`

For some of the tests, there are two tables of reference values, given by the ``refs`` and ``refs_alt`` attributes of the :class:`~tester.Tester`. The ``refs`` attribute always contains the recommended reference values. When ``refs_alt`` exists, any non-``None`` values represent values contained in the original documentation or Fortran version of TEOS-10 that disagree with ``teospy``. The reasons for diverging from the previous values are summarized in the section :ref:`issues-section`.

Finally, the module :mod:`testdoc` provides a way to run a modified version of
``doctest``. The standard ``doctest`` finds and runs example code in the
docstrings of Python programs. This modification extends that functionality to
testing floating point outputs only to within the significant figures of the
reference values. This module can be called as a function to test the example
codes in any module of ``teospy``.


``tester``
----------

.. automodule:: tester
    :members: Tester


``testair``
-----------

.. automodule:: testair


``testcnv``
-----------

.. automodule:: testcnv


``testflu``
-----------

.. automodule:: testflu


``testgsw``
-----------

.. automodule:: testgsw


``testice``
-----------

.. automodule:: testice


``testmix``
-----------

.. automodule:: testmix


``testsal``
-----------

.. automodule:: testsal


``testsea``
-----------

.. automodule:: testsea


``testdoc``
-----------

.. automodule:: testdoc
    :members: FltOutputChecker, testmod_flt

