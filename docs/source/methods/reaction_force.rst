Reaction Force Analysis
=======================


General Theory
--------------

Many of the methods implemented in pyrex center around the definition of a quantity known as the reaction force (:math:`F`). Analogous to the force defined in classical physics, the force for the reaction is defined as the negative gradient of the energy (:math:`E`) with respect to the reaction coordinate (:math:`\xi`)::
 
 .. math:: F(\xi) = - \frac{\partial E}{\partial \xi},

this produces a reaction force profile with a general shape defined by two critical points along the coordinate based on the force maximum and force minimum. This allows for a general partitioning of any reaction into three well-defined regions, a reactant region, transition state region, and a product region. In another analogy to classical physics one can integrate over the force in a particular region and obtain a reaction work (:math:`w`). For example, the first region of a reaction which would occur from the reactant structure (:math:`\xi_{\rm R}`) to the force minimum (:math:`\xi_{\rm min}`) would be defined as:

 .. math:: w_1 = - \int_{\xi_{\rm R}}^{\xi_{\rm min}} F(\xi) d\xi,

this is essentially the area under the curve of the force profile indicating the work done in the that region of the reaction. When interpreting the force profile, one can think of the negative forces as impeding reaction progress, while positive forces are seen as reaction driving.

Example
-------

PUT EXAMPLE HERE
 
