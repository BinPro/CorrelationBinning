#!/usr/bin/env python
from setuptools import setup, find_packages
import sys, os

version = '0.0'

setup(name='corrbin',
      version=version,
      description="Correlation Binning Module",
      long_description="""\
Useful functions and classes for correlation binning experiment.""",
      classifiers=[], # Get strings from http://pypi.python.org/pypi?%3Aaction=list_classifiers
      keywords='Python Scilifelab Metagenomics Binning Clustering Contig',
      author='Johannes Alneberg, Brynjar Smari Bjarnason',
      author_email='johannes.alneberg@scilifelab.se',
      url='www.github.com/binpro/corrbin',
      license='',
      packages=find_packages(exclude=['ez_setup', 'examples', 'tests']),
      scripts=[],
      include_package_data=True,
      zip_safe=False,
      install_requires=[
          # -*- Extra requirements: -*-
      ],
      entry_points="""
      # -*- Entry points: -*-
      """,
      )
