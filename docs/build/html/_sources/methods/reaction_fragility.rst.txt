Reaction Fragility Spectrum
===========================


General Theory
--------------

While atomic contributions to the reaction force and reaction force constant can be beneficial in tracking the breaking/formation of chemical bonds in simple reaction systems, it is limited in its application because the changes in both are often dominated by the change in position of the atom rather than its actual participation in a bond. In other words, for more complex cases it is difficult to discern the difference between the elongation/shortening of a bond and true bond breaking and formation for different types of atoms. The reaction fragility spectrum is better able to discern this difference by making use of the coupled normal modes of the molecular Hessian to obtain a quantitative measure of the bonding status of individual atoms in a molecule. If we consider a molecular system composed of :math:`N` atoms, the Hessian matrix element for the interaction between atoms :math:`i` and :math:`j` is defined by the following mixed partial derivative:

  .. math:: k_{ij} = \left(\frac{\partial^2 E}{\partial R_i \partial R_j}\right),  

where :math:`E` is the total molecular energy while :math:`R_i` and :math:`R_j` represent the cartesian position (:math:`R_x`, :math:`R_y`, :math:`R_z`) of atoms :math:`i` and :math:`j` respectively. Given the three cartesian coordinates, this would compose a :math:`3N \ x \ 3N` Hessian matrix (:math:`\mathbf{K}`). Every :math:`3 \ x \ 3` square matrix block along the block diagonal of the Hessian matrix is related to the individual atoms in the molecule. The trace of the Hessian is invariant to changes in coordinate and responds to changes in bonding structure of the molecule. The trace of each individual square block is that atoms contribution to the bonding pattern, this naturally lends itself to an atomic decomposition:

  .. math:: {\rm Tr}\mathbf{K} = \sum^N_A {\rm Tr} \mathbf{k}^A


