Activation Strain Model
=======================

General Theory
--------------

The activation strain model (ASM), also known as the "distortion interaction model" is a rather simple energy decomposition scheme that lays the foundation for many of the techniques implemented in pyrex. Considering a reaction coordinate :math:`\xi`, the potential energy surface :math:`\Delta E(\xi)` can be decomposed into two contributions along the coordinate: (1) the strain energy :math:`\Delta E_{\rm strain}(\xi)`, which is associated with the structural deformation of the reactants as the reaction proceeds and (2) the interaction energy :math:`\Delta E_{\rm int}(\xi)`, which is associated with the electronic interactions between the electrons of the reactants as the reaction proceeds. This yields the following decomposition::
 
    .. math:: \Delta E(\xi) = \Delta E_{\rm strain}(\xi) + \Delta E_{\rm int}(\xi),

in general the strain energy is a positive quantity, and thus destabilizing with respect to the separated reactants. Analysis of the interplay between these two energy quantities are directly related to the observed barrier height of a given reaction. This type of analysis has been used to effectively analyze barrier heights in numerous reaction mechanisms.

Example
-------

For an ASM calculation you will have to specify fragments for your system, the following example on the reaction between carbon dioxide and the hydrogen molecule specifies each reactant as a special fragment. The "do_polarization" keyword is added in order to perform the ASM calculations::

    {
      "molecule": {
        "symbols": ["C","O","O","H","H"],
        "molecular_charge": "0",
        "molecular_multiplicity": 1,
        "fragments" : [[0,1,2],[3,4]],
        "fragment_charges" : [0,0],
        "fragment_multiplicities" : [1,1]
      },
      "model": {
        "method": "mp2",
        "basis": "def2-svp"
      },
      "pyrex": {
        "nthreads": 4,
        "irc_filename": "full_irc.xyz",
        "do_energy": true,
        "do_conceptualdft": true,
        "do_polarization" : true,
        "irc_stepsize": 0.05
      }
    }

When running this calculation, the first thing pyrex will do is perform a geometry optimization on each of your specified fragments in isolation. Then is will calculate the energy of the supermolecule along the coordinate, then the energies of each isolated fragment (with ghost functions accounting for the other fragment to account for basis set superposition errors). After the calculation finishes, pyrex will produce a CSV file called "activation_strain.csv" that contains the interaction energy, strain energy, and total relative energy. We can easily plot the results of our calculation by using the following REXplot input::

    {
      "rexplot" : {
        "file" : "activation_strain.csv",
        "properties" : ["Energy", "Interaction Energy", "Strain Energy"],
        "coordinate" : "Coordinate",
        "x_label" : "Reaction Coordinate ($\\xi$)",
        "y_label" : "$\\Delta E$ (kcal/mol)",
        "scale" : 627.509,
        "fig_dims" : [9.0, 5.0],
        "plot_file" : "activation_strain.png"
      }
    }

The "scale" option is optional, and has only been added to convert the energies from Hartrees to kcal/mol. This will produce a plot that contains all of the relevant energies for an activation strain analysis. Let's investigate the plot:

.. image:: figures/activation_strain.png
   :width: 500 px
   :alt: alternate text
   :align: center

As you can see from inspecting the plot, the destabalizing interaction energy is the primary contributor to the activation energy barrier with little contribution from the strain energy. It isn't until you have significant elongation of the bond in hydrogen that you see a significant strain contribution. This simple example highlights the utility of activation strain analysis and how to run these calculations in pyrex. 
