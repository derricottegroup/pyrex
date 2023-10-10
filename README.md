![pyrex](logos/pyrex_logo_2019_banner.png)

[![Travis-CI](https://travis-ci.org/WDerricotte/pyrex.svg?branch=master)](https://travis-ci.org/WDerricotte/pyrex)
[![Documentation Status](https://readthedocs.org/projects/pyrex/badge/?version=latest)](http://pyrex.readthedocs.io/en/latest/?badge=latest)
[![License](https://img.shields.io/github/license/WDerricotte/pyrex)](https://github.com/WDerricotte/pyrex/blob/master/LICENSE)
[![Commits](https://img.shields.io/github/commit-activity/m/WDerricotte/pyrex)](https://github.com/WDerricotte/pyrex/commits/master)
[![Twitter](https://img.shields.io/twitter/follow/ProfDerricotte?style=social&logo=twitter)](https://twitter.com/ProfDerricotte)

#### Code Authors: Derricotte Research Group

#### Manual (http://pyrex.readthedocs.io/)
#### Web: www.derricotteresearchgroup.com

### Overview

pyREX(Python Reaction Energy eXtension), is a free open-source implementation of reaction coordinate analysis techniques that interfaces with PSI4 ab-initio quantum chemistry software in order to streamline the process of investigating energetic/electronic properties along an intrinsic reaction coordinate. This code is currently under development within the Derricotte Research Group at Morehouse College in Atlanta, GA. Current Key features include:(1) Calculating SCF energies along the reaction coordinate, (2) Reaction Force Analysis, (3) Reaction Electronic Flux (REF) Analysis, (4) Decomposition of REF into Polarization and Transfer Components, (5) Symmetry Adapted Perturbation Theory Decomposition analysis

### Installation
#### Original Instructions
1. Obtain required software
    1. [pyREX](https://github.com/WDerricotte/pyrex) (must clone this repository; no binary install currently available)
    2. [Psi4](http://psicode.org/psi4manual/1.1/build_obtaining.html) pyREX depends intrinsically on Psi4, currently the easiest way to obtain the code is to download the available binaries available at http://vergil.chemistry.gatech.edu/psicode-download/1.1.html
    3. Export the necessary paths either through the command line or .bashrc
    ``` 
    export PATH=$PATH:/path/to/pyrex/src/directory/pyrex:
    export PYTHONPATH=$PYTHONPATH:/path/to/pyrex/src/directory/pyrex:
    ```
    4. The main "pyrex" module should now be executable from anywhere in your directory tree by simply typing "pyrex".  

#### Dom's Instructions
0. Clone this repository
1. Install Dependencies (via `conda`)
    - Create a clean environment, install required dependencies
    ```
    ~$ conda create -n rex python=3.7 jupyter seaborn
    ~$ conda activate rex
    (rex) ~$ conda install pyscf -c pyscf
    (rex) ~$ conda install psi4 -c psi4
    ```
2. Install PyREX in new environment
    - Navigate to the local clone of this repository, then execute
    ```
    (rex) /path/to/pyrex ~$ pip install -e .
    ```

