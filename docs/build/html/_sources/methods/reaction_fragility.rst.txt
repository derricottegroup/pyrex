Reaction Fragility Spectrum
===========================


General Theory
--------------

While atomic contributions to the reaction force and reaction force constant can be beneficial in tracking the breaking/formation of chemical bonds in simple reaction systems, it is limited in its application because the changes in both are often dominated by the change in position of the atom rather than its actual participation in a bond. In other words, for more complex cases it is difficult to discern the difference between the elongation/shortening of a bond and true bond breaking and formation for different types of atoms. The reaction fragility spectrum is better able to discern this difference by making use of the coupled normal modes of the molecular Hessian to obtain a quantitative measure of the bonding status of individual atoms in a molecule. If we consider a molecular system composed of :math:`N` atoms, the Hessian matrix element for the interaction between atoms :math:`i` and :math:`j` is defined by the following mixed partial derivative:

  .. math:: k_{ij} = \left(\frac{\partial^2 E}{\partial R_i \partial R_j}\right),  

where :math:`E` is the total molecular energy while :math:`R_i` and :math:`R_j` represent the cartesian position (:math:`R_x`, :math:`R_y`, :math:`R_z`) of atoms :math:`i` and :math:`j` respectively. Given the three cartesian coordinates, this would compose a :math:`3N \ \times \ 3N` Hessian matrix (:math:`\mathbf{K}`). Every :math:`3 \ \times \ 3` square matrix block along the block diagonal of the Hessian matrix is related to the individual atoms in the molecule. The trace of the Hessian is invariant to changes in coordinate and responds to changes in bonding structure of the molecule. The trace of each individual square block is that atoms contribution to the bonding pattern, this naturally lends itself to an atomic decomposition:

  .. math:: {\rm Tr} \ \mathbf{K} = \sum^N_A {\rm Tr} \ \mathbf{k}^A,

where :math:`\mathbf{k}^A` is the square block diagonal matrix related to atom :math:`A`. Any decrease (increase) in the trace of :math:`\mathbf{k}^A` corresponds to the weakening (strengthening) of a bond. The reaction fragility spectrum (:math:`a_{\xi}`) tracks the change in this contribution as the reaction proceeds, it does this by investigating the derivative of the trace of the Hessian as the reaction progresses: 

  .. math:: a_{\xi} = \frac{\partial}{\partial \xi}({\rm Tr} \ \mathbf{K}) = \sum^N_A \frac{\partial}{\partial \xi} {\rm Tr} \ \mathbf{k}^A = \sum^N_A a^A_{\xi},

where (:math:`a^A_{\xi}`) would be the atomic fragility spectrum for atom :math:`A`. By utilizing the invariant trace of the Hessian rather than atomic contributions to the total energy, the fragility spectrum represents a tool that is less sensitive to nuclear motion and more akin to bond indices. 

Example
-------

The example below is an input file to calculate the reaction fragility spectrum for the reaction of carbon dioxide and the hydrogen molecule. The "do_fragility_spec" keyword is added to the "pyrex" block in order to calculate the fragility spectrum::

    {
      "molecule": {
        "symbols": ["C","O","O","H","H"],
        "molecular_charge": "0",
        "molecular_multiplicity": 1
      },
      "model": {
        "method": "scf",
        "basis": "sto-3g"
      },
      "pyrex": {
        "nthreads": 4,
        "irc_filename": "full_irc.xyz",
        "do_energy" : true,
        "energy_read" : "energy.csv",
        "do_fragility_spec" : true,
        "irc_stepsize": 0.2
      }
    }

the "energy_read" keyword is added in this case so that the total energies can be read in from a previously calculated file, avoiding the need to recalculate the energy. This calculation will produce a file called "frag_spec.csv", that will contain the atomic fragility spectrum for each atom in the molecule. To highlight this method, we will take a look at the fragility spectra for the two oxygens, the figure below shows the molecular structure of the reactants, transition state, and product for this reaction:

.. image:: figures/co2_h2_structs.png
   :width: 500 px
   :alt: alternate text
   :align: center

Notice that O1 forms a new bond with a hydrogen atom while O2 does not, based on our understanding of the reaction fragility spectrum we would expect a much larger change in this spectrum for O1 when compared to O2. However, we will have to investigate this to find out. Using the plotting utility, we can directly plot the two atomic fragility spectra on the same plot using the following input::

    {
      "rexplot" : {
        "file" : "frag_spec.csv",
        "properties" : ["O1", "O2"],
        "coordinate" : "Coordinate",
        "x_label" : "Reaction Coordinate ($\\xi$)",
        "y_label" : "$a_{\\xi}$",
        "fig_dims" : [9.0, 5.0],
        "plot_file" : "frag_spec.png"
      }
    }

This will produce the following plot of both atomic fragility spectra:

.. image:: figures/frag_spec.png
   :width: 500 px
   :alt: alternate text
   :align: center
