Reaction Fragility Spectrum
===========================


General Theory
--------------

While atomic contributions to the reaction force and reaction force constant can be beneficial in tracking the breaking/formation of chemical bonds in simple reaction systems, it is limited in its application because the changes in both are often dominated by the change in position of the atom rather than its actual participation in a bond. In other words, for more complex cases it is difficult to discern the difference between the elongation/shortening of a bond and true bond breaking and formation for different types of atoms. The reaction fragility spectrum is better able to discern this difference by making use of the coupled normal modes of the molecular Hessian to obtain a quantitative measure of the bonding status of individual atoms in a molecule. If we consider a molecular system composed of :math:`N` atoms, the Hessian matrix element for the interaction between atoms :math:`i` and :math:`j` is defined by the following mixed partial derivative::

    .. math:: k_{ij} = \left(\frac{\partial^2 E}{\partial R_i \partial R_j}\right),  


