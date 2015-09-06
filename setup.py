#!/usr/bin/env python2

from setuptools import setup, find_packages

requirements = [
    'numpy>1.07',
    'scipy>0.10'
    ]

setup(
    name="pyrad",
    version="3.1.0a5",
    url="https://github.com/dereneaton/pyrad",

    author="Deren Eaton",
    author_email="deren.eaton@yale.edu",

    description="Assembly and analysis of RADseq data sets",
    long_description=open('README.rst').read(),

    packages=find_packages(),
    
    install_requires=[requirements],

    entry_points={
            'console_scripts': [
                'pyrad = pyrad.pyRAD:main',
            ],
    },

    license='GPL',

    classifiers=[
        'Programming Language :: Python',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
    ],
)
