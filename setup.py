#!/usr/bin/env python
from lib.__version__ import version

from setuptools import setup, find_packages
# from distutils.core import setup
from distutils.extension import Extension
#from Cython.Build import cythonize

with open('README.md') as f:
    long_description = f.read()


setup(
    name='REPcluster',
    version=version,
    description='subphaser: phase subgenomes of allopolyploidy based on repeatitive kmers',
    url='https://github.com/zhangrengang/REPcluster',
    author='Zhang, Ren-Gang and Wang, Zhao-Xuan',
    license='GPL-3.0',

    python_requires='>=3.6:',
    packages=find_packages(),
    include_package_data=True,
    scripts=[],
    entry_points={
        'console_scripts': ['REPclust = lib.cluster:main',
        ],
    },
)
