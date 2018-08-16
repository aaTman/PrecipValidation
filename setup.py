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
    description='An example python package',
    long_description=open('README.txt').read(),
    install_requires=['numpy'],
    url='https://github.com/BillMills/python-package-example',
    author='Bill Mills',
    author_email='myemail@example.com'
)