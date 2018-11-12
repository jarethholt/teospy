.. teospy documentation master file, created by
   sphinx-quickstart on Tue Aug 28 09:58:03 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to teospy's documentation!
==================================

``teospy`` is a pure-Python implementation of the `Thermodynamic Equation of
Seawater (2010) <http://www.teos-10.org>`_, other implementations of which can
be found at `https://github.com/TEOS-10 <https://github.com/TEOS-10>`_. The
library provides thermodynamic functions of mixtures of dry air, water vapour,
liquid water, ice, and sea salt. Examples uses include calculating the density
of air from its temperature and pressure, or the vapour pressure above seawater
with a given salinity, temperature, and pressure.

The purpose of this pure-Python implementation is to open up the library to
interactive calculations and scripting without needing to compile other
languages. It also allows for limited integration with :mod:`numpy` for arrays
and :mod:`scipy` for root-finding and optimization.

This library is divided into modules by substance. The ``air``, ``flu``,
``ice``, and ``sea`` modules contain the functions pertaining to dry and humid
air; liquid water and water vapour (fluid water collectively); pure ice; and
seawater (liquid water and salt), respectively. Modules with more than one
substance handle equilibrium between the three phases of water, for example
``iceliq4`` for pure liquid water and ice, or ``seaair4`` for seawater and humid
air.

In addition, the modules are split into several levels as follows:

* The level 0 modules (:mod:`constants0` and :mod:`convert0`) are used by many
  different higher-level modules and across substances.
* Level 1 modules implement the primary thermodynamic quantities for each
  substance: the Helmholtz energy for dry air and fluid water and the Gibbs
  energy for ice and sea salt.
* Level 2 modules provide thermodynamic quantities that can be calculated
  directly from the primary quantities using the same primary variables.
* Level 3 modules allow properties to be calculated from secondary variables.
  They rely on the module :mod:`maths3` which provides a root-finding method.
  The liquid water and water vapour phases are differentiated in these modules,
  and liquid water and sea salt are combined to make seawater modules.
* Level 4 modules are for mixtures of substances with the phase coexistence of
  water. Here you will find functions relating to saturation vapour pressure,
  condensation pressures of ice and water, sea ice melting temperature, etc.
* Level 5 modules implement a few extra features. :mod:`convert5` provides
  conversions between non-SI units. Some modules implement alternative
  approaches to calculating properties, such as a Gibbs energy formulation
  (:mod:`liq5_f03`). Most importantly, :mod:`sea5` implements the recommended
  routines for conservative temperature and the related thermobaric, cabelling,
  thermal expansion, and haline contraction coefficients. :mod:`gsw5` provides
  the same quantities, but using the :mod:`liq5_f03` formulation and the Gibbs
  SeaWater (GSW) toolbox inputs.
* The various `test` modules provide tests of this library, comparing its
  results to the original Fortran. All of these modules make use of the
  :class:`~tester.Tester` class. This class is just a container to keep track
  of the functions to test, their inputs, reference results, and formatting.
  It also provides methods to run the functions and print summaries.


Table of Contents
=================

.. toctree::
   :maxdepth: 2
   
   root
   air
   flu
   ice
   sal
   sea
   mix
   tester
   issues


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

