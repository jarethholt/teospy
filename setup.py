"""Setup script for teospy.
"""

import setuptools

with open('README.rst','r') as readme:
    long_description = readme.read()

setuptools.setup(
    name = 'teospy',
    version = '0.0.1',
    author = 'Jareth Holt',
    author_email = 'jareth.holt@gmail.com',
    description = 'Thermodynamic library for air, ice, and seawater',
    long_description = long_description,
    url = 'https://github.com/teospy',
    packages = setuptools.find_packages(),
    data_files = [('teospy',['GSW_Data_v3_0.dat'])],
    include_package_data = True,
    classifiers = [
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 3',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Atmospheric Science',
        'Topic :: Scientific/Engineering :: Physics'
    ],
)