.. teospy introduction documentation

Intro
=====

`teospy` is a pure-Python implementation of the `Thermodynamic Equation of Seawater (2010) <http://www.teos-10.org>`_, which can be found in other languages at `https://github.com/TEOS-10 <https://github.com/TEOS-10>`_. The library provides thermodynamic functions of mixtures of dry air, water vapour, liquid water, ice, and sea salt. Examples uses include calculating the density of air from its temperature and pressure, or the vapour pressure of water above seawater from its temperature, pressure, and salinity.

The purpose of this pure-Python implementation is to open up the library to interactive calculations and scripting without needing to compile other languages. It also allows for limited integration with `numpy` for arrays and `scipy` for root-finding and optimization.

This library is divided into several sub-modules. The `air`, `flu`, `ice`, and `sea` modules contain the functions pertaining to air (dry and humid), pure fluid water (liquid and vapour), pure ice, and seawater (liquid water and salt). The module `mix` contains functions for mixtures of the above components and focuses on calculating equilibrium between phases.

In addition, the modules are split into several levels as follows:

* The level 0 modules (`constants0` and `convert0`) are used by many different higher-level modules and across substances.
* Level 1 modules implement the primary thermodynamic quantities for each substance: the Helmholtz energy for dry air and fluid water and the Gibbs energy for ice and sea salt.
* Level 2 modules provide thermodynamic quantities that can be calculated directly from the primary quantities, i.e. using the same primary variables.
* Level 3 modules allow properties to be calculated from secondary variables. They rely on the module `maths3` which provides a root-finding method. The liquid water and water vapour phases are differentiated in these modules, and liquid water and sea salt are combined to make seawater modules.
* Level 4 modules are for mixtures of substances with the phase coexistence of water. Here you will find functions relating to saturation vapour pressure, condensation pressures of ice and water, sea ice melting temperature, etc.
* Level 5 modules implement a few extra features. Some modules provide ways to calculate quantities in non-SI units, particularly `convert5`. Some modules implement alternative approaches to water properties, such as a Gibbs energy formulation (`liq5f03`). Most importantly, `sea5` implements the recommended routines for conservative temperature and the related thermobaric, cabelling, thermal expansion, and haline contraction coefficients.
* Modules with the prefix `values` keep tables of reference values and are used to test library performance.




























