Getting Started
===============

Before starting a typical calculation in pyrex you will need an XYZ coordinate file that contains the geometries associated with an intrinsic reaction coordinate (IRC) or some sort of surface scan/potential energy surface. When you have that file, pyrex can compute multiple energies and properties along this coordinate. Pyrex uses the JSON file format for inputs, a standard input file will have the following format::

   {
     "molecule": {
       "symbols": ["C","O","O","H","H"],
       "molecular_charge": "0",
       "molecular_multiplicity": 1
     },
     "model": {
       "method": "b3lyp",
       "basis": "cc-pvdz"
     },
     "pyrex": {
       "nthreads": 4,
       "irc_filename": "full_irc.xyz",
       "do_energy": true,
       "irc_stepsize": 0.2
     }
   }

For this minimal case, the molecule is specified first, with the only necessary items being the symbols corresponding to each atom in your system (these should correspond with the order the atoms appear in your coordinate file). After that the molecular charge and multiplicity must be specified. The next "model" block contains the options related to the level of theory you wish to calculate the energy. The final "pyrex" block contains pyrex specific options for your calculations, "nthreads" controls the number of threads to be used, "irc_filename" lets pyrex know the name of your coordinate file, "do_energy" instructs it to calculate the energy, and "irc_stepsize" gives the stepsize of your irc. Running this output will return the reaction energy and the reaction force analysis along the coordinate. 

This is the most common use case for pyrex, however there are many other features of the code that will be covered in the other sections of the documentation.  
