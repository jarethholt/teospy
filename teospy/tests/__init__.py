"""Testing modules for teospy.

These submodules provide tests of ``teospy``, typically as comparisons
to results from the original Fortran routines. These modules use the
:class:`~teospy.tests.tester.Tester` class which is just a container
holding a set of functions, inputs, and reference values to test. In
addition, some of the base modules in ``teospy`` use this class for
their ``chk`` functions -- corresponding to the Fortran ``Values_*``
programs -- which also serve as tests. Finally,
:mod:`~teospy.tests.testdoc` offers docstring testing (based on the
standard Python ``doctest``) where floating point values are only
evaluated within the given significant figures.
"""

