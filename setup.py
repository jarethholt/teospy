"""Setup script for teospy.
"""

import setuptools

with open('README.rst','r') as readme:
    long_description = readme.read()

setuptools.setup(
    name = 'teospy',
    version = '0.0.1',
    packages = setuptools.find_packages(),
    install_requires = ['numpy'],
    include_package_data = True,
    author = 'Jareth Holt',
    author_email = 'jareth.holt@gmail.com',
    description = 'Thermodynamic library for air, ice, and seawater',
    long_description = long_description,
    license = 'MIT',
    url = 'https://github.com/teospy',
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

