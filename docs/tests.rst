.. teospy tester module documentation

Testing modules
===============

Several modules are provided in the ``tests`` folder for testing ``teospy``.
Most of them make use of the :class:`~teospy.tests.tester.Tester` class. This
class is initialized with a list of functions to test, values of the input
arguments, and an array of reference values. It also takes several arguments
related to printing the results. Use :func:`~teospy.tests.tester.Tester.run`
to calculate results from this library and relative errors. Then
:func:`~teospy.tests.tester.Tester.printresults` prints the function names,
arguments, reference value, result, and relative errors whenever the relative
error is larger than a given tolerance.

Several tests have been compiled into modules such as
:mod:`~teospy.tests.testair` to readily separate the tests by substance. These
modules correspond to modules such as ``Values_Air`` in the Fortran version of
TEOS. The :class:`~teospy.tests.tester.Tester` class is also used in the
``chk`` functions of some modules, and those modules can be called from the
command line (i.e. as a script) to run these checks. The modules with such
``chk`` functions are:

.. hlist::
    :columns: 4
    
    * :mod:`~teospy.air2`
    * :mod:`~teospy.air3b`
    * :mod:`~teospy.flu1`
    * :mod:`~teospy.flu2`
    * :mod:`~teospy.flu5_if97`
    * :mod:`~teospy.ice2`
    * :mod:`~teospy.liq5_f03`
    * :mod:`~teospy.liqvap4`
    * :mod:`~teospy.sea3a`

For some of the tests, there are two tables of reference values, given by the
``refs`` and ``refs_alt`` attributes of the
:class:`~teospy.tests.tester.Tester`. The ``refs`` attribute always contains
the recommended reference values. When ``refs_alt`` exists, any non-``None``
values represent values contained in the original documentation or Fortran
version of TEOS-10 that disagree with ``teospy``. The reasons for diverging
from the previous values are summarized in :ref:`issues-section`.

Finally, the module :mod:`~teospy.tests.testdoc` provides docstring testing
similar to ``doctest``. The standard ``doctest`` finds and runs example code
in the docstrings of Python programs. This modification extends that
functionality to testing floating point outputs only to within the significant
figures of the reference values. This module can be called as a function to
test the example codes in any module of ``teospy``.


``tester``
----------

.. automodule:: teospy.tests.tester
    :members: Tester


``testair``
-----------

.. automodule:: teospy.tests.testair


``testcnv``
-----------

.. automodule:: teospy.tests.testcnv


``testflu``
-----------

.. automodule:: teospy.tests.testflu


``testgsw``
-----------

.. automodule:: teospy.tests.testgsw


``testice``
-----------

.. automodule:: teospy.tests.testice


``testmix``
-----------

.. automodule:: teospy.tests.testmix


``testsal``
-----------

.. automodule:: teospy.tests.testsal


``testsea``
-----------

.. automodule:: teospy.tests.testsea


``testdoc``
-----------

.. automodule:: teospy.tests.testdoc
    :members: FltOutputChecker, testmod_flt

