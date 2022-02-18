#!/usr/bin/env python3

from setuptools import find_packages, setup
import os

here = os.path.abspath(os.path.dirname(__file__))

# Get the long description from the README file
with open(os.path.join(here, 'README.md')) as f:
    long_description = f.read()

# get the dependencies and installs
with open(os.path.join(here, 'requirements.txt')) as f:
    all_reqs = f.read().split('\n')

install_requires = [x.strip() for x in all_reqs[:-1] if 'git+' not in x]

setup(name='hierarchical-chain-growth',
      version='0.1',
      description='Grow Ensembles of Disordered Biomolecules from Fragment Libraries',
      url='https://github.com/bio-phys/hierarchical-chain-growth',
      author='Lisa M. Pietrek, Lukas S. Stelzl, Gerhard Hummer',
      author_email='lisapietrek@gmail.com',
      install_requires=install_requires,
      license='GPLv3',
      packages=find_packages(exclude=['examples']),
      zip_safe=False)
