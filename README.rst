teospy
======

Python implementation of the Thermodynamic Equation of Seawater (2010) for seawater, ice, and air.

The Thermodynamic Equation of Seawater (2010), TEOS10, is a framework for calculating thermodynamic properties of seawater, pure vapour/liquid/ice water, dry air, and mixtures thereof. More information on the framework and scientific background is available at `http://www.teos-10.org/ <http://www.teos-10.org/>`_. The associated software libraries in other languages (Fortran, Visual Basic, C) are available at `https://github.com/TEOS-10 <https://github.com/TEOS-10>`_.

This library is a pure-Python implementation of the TEOS10 Seawater-Ice-Air library. Its purpose is to avoid compiling interfaces to other languages so that calculations can easily be done interactively or as part of a script. Each phase equilibrium problem is cast in the same format, so that a single root-finding framework can be used. Future developments aim at expanding support for `numpy` arrays and `scipy` optimization modules.

The documentation for this project is maintained using Sphinx. The RST files are included in the `docs` folder, and the full HTML is at `https://teospy.readthedocs.io/ <https://teospy.readthedocs.io/>`_.


Credits
-------

* `teospy` author: Jareth Holt
* David Jackett, Paul Barker, Glenn Hyland, and more for development of the GSW and SIA toolboxes
* R. Feistel, T. McDougall, J. Reissmann, and more for creating and maintaining the TEOS scientific standards
* The Joint Committee on the Properties of Seawater and its parent organizations IAPWS, SCOR, and IAPSO


License
-------

This project is licensed under the MIT license. See the LICENSE file for more details.





