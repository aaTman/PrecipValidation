# -*- coding: utf-8 -*-
"""
Created on Wed Aug  1 17:43:48 2018

@author: Taylor
"""

from setuptools import setup, find_packages

setup(
    name='weightvalidate',
    version='0.1',
    packages=find_packages(exclude=['tests*']),
    license='MIT',
    description='A python package to ingest model run data (initally precipitation) and compare to verification datasets. Different parameters and weights can be placed into the package to customize the strictness of the validation approach. The main approach will utilize a precribed radius of n kilometers. Validation will weight precipitation data based on distance from the centroid, temporal disparity, and other factors which will be addressed based on forthcoming input.',
    long_description=open('README.txt').read(),
    install_requires=['numpy','pygeodesy','datetime','os','netCDF4'],
    url='https://github.com/aaTman/precipValidation',
    author='Taylor Mandelbaum',
    author_email='mandelbaum.taylor@gmail.com'
)