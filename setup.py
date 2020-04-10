#!/usr/bin/env python

import sys
from setuptools import setup

setup(name='pyrex',
      version='1.1.0',
      packages=['pyrex'],
      entry_points={
          'console_scripts': [
              'pyrex = pyrex.__main__:main'
          ]
      },
      )
