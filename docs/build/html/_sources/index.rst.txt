Overview
========

pyREX(Python Reaction Energy eXtension), is a free open-source implementation of reaction coordinate analysis techniques that interfaces with popular quantum chemistry software in order to streamline the process of investigating energetic/electronic properties along an intrinsic reaction coordinate. This code is currently under development within the Derricotte Research Group at Morehouse College in Atlanta, GA.

Capabilities
------------

pyREX implements a series of reaction energy methods based on the total reaction energy, conceptual DFT properties, derivative properties along the reaction coordinate, and energy decomposition techniques. The following methods are available in pyREX:

#. Intrinsic Reaction Coordinate (IRC) following code 
#. Activation Strain Model (ASM) Decomposition of reaction energy
#. Reaction Force 
#. Conceptual DFT Properties (chemical potential, chemical hardness, etc.)
#. Reaction Electronic Flux
#. Atomic Decomposition of Reaction Force
#. Reaction Fragility Spectrum
#. Symmetry-Adapted Perturbation Theory Decomposition of the Reaction Force

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   getting_started
   methods/activation_strain
   methods/reaction_force
   input_generator

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
