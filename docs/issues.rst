.. teospy documentation: known issues

.. _issues-section:

Known issues
============

This section lists known issues with ``teospy``. This includes discrepancies
with the Fortran version or stated documentation; undesirable behavior with
default settings; or missing or unintuitive implementation of functions.


Seawater standard conditions
----------------------------

``teospy`` does not replicate the values of certain seawater properties (Gibbs
energy, entropy) at standard seawater conditions: absolute salinity 0.03516504
kg/kg, temperature 273.15 K, pressure 101325 Pa. This is partially because such
values are supposed to be nearly zero by design but rely on two much larger
values cancelling. The values given in the tests in :mod:`~teospy.sea3a` are
the same as given in the Fortran version, and only have one significant figure.
The Fortran implementation recreates these values to 3 significant figures
under default settings, but not if the tolerance in the fluid equilibrium
function is reduced to 1e-12. However, these low-tolerance values are not
replicated by ``teospy`` either.

In the TEOS-10 manual, it is stated that the free parameters in the salinity
function are chosen such that the enthalpy and entropy of seawater at the
standard conditions are zero, when the Feistel (2003) formulation for the Gibbs
energy of liquid water is used. I have included a way to print salinity
coefficients compatible with this condition in the script
:mod:`~teospy.calsalcoeffs`, but doing so causes many more disagreements with
seawater properties under other conditions.

It is unclear to me how this agreement could be improved or whether it would
even be meaningful, given the nature of cancellation with floating point
arithmetic. The properties at other (salinity, temperature, pressure) values
agree with the documentation.



Comparison to Fortran implementation
------------------------------------

I am listing here discrepancies between ``teospy`` and the TEOS Fortran
implementation. The specific version I am referencing is from 31 Jan 2012; any
of the concerns raised below may be invalid for newer versions.


Enthalpy of evaporation
^^^^^^^^^^^^^^^^^^^^^^^

Several of the equilibrium modules in the original Fortran (e.g. ``Ice_Air_4a``)
use an incorrect formula for the enthalpy of evaporation/sublimation/melting
into humid air. When calculating the change in pressure with dry air mass
fraction, they have :math:`\partial p/\partial a = d_h f^h_{ad}` instead of
:math:`\partial p/\partial a = d_h^2 f^h_{ad}`. The resulting discrepancies are
typically of order 50 J/kg, or a relative error of 1e-5. The enthalpy values
used in the ``test`` modules come from a modified Fortran implementation. The
modules affected are :mod:`~teospy.iceair4a`, :mod:`~teospy.liqair4a`, and
:mod:`~teospy.seaair4`.


Gibbs energy dry fraction derivatives
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Similar to the enthalpy of evaporation, several of the mixture modules with
humid air use an incorrect formula for the latent heat term in the Gibbs energy
derivatives. Specifically, this affects the second-order derivatives with
respect to temperature and pressure. A latent term of the form :math:`-d_h a_x
a_y/a`, where `x` and `y` are either temperature or pressure, should be
:math:`-d_h a_x a_y/a^2`. The relative errors in the Gibbs energy derivatives
are typically of order 2e-4. These errors propagate to their related properties
(heat capacity, thermal expansion coefficient, compressibility, and lapse rate)
as well as to the enthalpy derivatives. The modules affected are
:mod:`~teospy.iceair4b`, :mod:`~teospy.iceair4c`, :mod:`~teospy.liqair4b`, and
:mod:`~teospy.liqair4c`.


Gibbs energy salinity derivatives
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The Gibbs energy mixed salinity-temperature and salinity-pressure derivatives
for seawater mixtures have an error. These derivatives include terms such as
:math:`g^b_{sx} + g^b_{ss} s_x`, where `x` is temperature or pressure. The
Fortran version has :math:`g^b_{sx} + g^b_s s_x` instead. This affects the
:mod:`~teospy.seavap4` and :mod:`~teospy.seaice4` modules. In particular,
notice that the thermal expansion coefficient of sea-ice parcels,
:func:`~teospy.seaice4.expansion`, should be negative with the given values.


(s,p)-derivative of ``iceair_h``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

There was a minor error in the Fortran module ``Ice_Air_4c`` enthalpy function.
In the mixed entropy-pressure derivative, the value for ``g_tp`` was not set and
hence missing from the resulting ``h_sp``. Only this value within
:func:`~teospy.iceair4c.iceair_h` is affected.


Density of liquid-ice-air
^^^^^^^^^^^^^^^^^^^^^^^^^

An incorrect formula was used for the density of liquid water-ice-humid air
mixtures in the Fortran version. The densities were added by mass fraction
rather than the specific volumes. This only affects
:func:`~teospy.liqiceair4.density` values when mass fractions or entropy and
wet fraction are given.


Direct temperature conversion
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In the Fortran module ``Convert_5``, all temperatures are first converted to
ITS-90 and then to the desired output units. However, the conversion functions
between ITS-90 and IPTS-68 are not exact inverses of each other. This introduces
errors when converting IPTS-48 or IPTS-68 temperatures to IPTS-68.
Interestingly, the values given in the documentation are correct, as if only one
conversion was done, but these values are *not* replicated by the Fortran
version.


Error in pressure conversion
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The Fortran version has an error in :mod:`~teospy.convert5.cnv_pressure` in
which converting from pressure in atmospheres uses the scaling constant for
kilogram-equivalent force.


Salinity conversion reference values
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The reference values given for the salinity conversion routines in
:mod:`~teospy.convert0` and hence :mod:`~teospy.convert5` and
:mod:`~teospy.gsw5` do not match the output from the current Fortran version.
They have been replaced by those values.


GSW toolbox reference values
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Similar to salinity, the reference values given for the :mod:`~teospy.gsw5`
functions referencing potential and conservative temperature (expansion,
contraction, cabbeling, and thermobaric coefficients) do not match the output
from the current Fortran version. They have been replaced by those values.


Increased accuracy of test values
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The original reference values for equilibrium quantities match the Fortran
implementation output with the default tolerances, which are typically 1e-7. In
practice, the resulting equilibrium state often satisfies equilibrium to an even
better degree, but not always. ``teospy`` with default settings uses a lower
tolerance (1e-8) but often reaches equilibrium agreement near 1e-12. For several
tests, I have replaced the original reference values with results from a Fortran
implementation with tolerance lowered to 1e-12. The tests affected are:

* :mod:`~teospy.iceair4a` and :mod:`~teospy.iceair4c`, temperature;
* :mod:`~teospy.seavap4`, for the enthalpy and entropy of seawater;
* :mod:`~teospy.seaicevap4`, for salinity at a given pressure;
* :mod:`~teospy.seaair4`, for condensation temperature and humid air density;
  and
* :mod:`~teospy.gsw5`, for the expansion, contraction, cabbeling, and
  thermobaric coefficients with respect to both conservative and potential
  temperature.

