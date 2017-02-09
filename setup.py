# -*- coding: utf-8 -*-
from __future__ import print_function

try:
    from setuptools import setup, find_packages, Extension
except ImportError:
    from distutils.core import setup, Extension

    def find_packages():
        return ['sidekick']

import pkg_resources
import platform
import re
import subprocess

# Install splash
VERSION = '0.0.0.9999'

setup(name="sidekick",
      version=VERSION,
      author='Kevin Gori',
      author_email='kg8@sanger.ac.uk',
      description='Find abandoned sidekick reads in bam files',
      url='https://github.com/TransmissibleCancerGroup/sidekick.git',
      packages=find_packages(),
      
      scripts=[
          'bin/sidekick.py'
      ],
      install_requires=[
          'pysam'
      ],
      test_suite='tests',
)
