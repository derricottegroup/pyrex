Input Generator
===============

In order to get started using pyREX you must have a .xyz coordinate file that contains the geometries for the reaction/process you are investigating and a .json input file that contains the relevant options for the calculation you would like to do in pyREX. There is a built-in input file generator that can help you with all the steps of getting your input file set up.

Joining IRC Files
-----------------

Often times the coordinate file you would like to run in pyREX is an intrinsic reaction coordinate (IRC) file. Most IRC codes output separate coordinates for the forward and backward trajectories for the reaction. However for calculations in pyrex you will need a continuous trajectory that begins with the last structure of your backward trajectory and ends with the last structure of your forward trajectory. You can do this very simply from the command line using the following command line utility::

    pyrex --join_irc

pyrex will then ask you the name of your forward and backward trajectory and produce a single "full_irc.xyz" file.

Generating Pyrex Input Files
----------------------------

Once you have an IRC or other coordinate file, the standard use of pyrex is to obtain energies and properties along the reaction coordinate. In order to generate a standard input for the reaction energy, you can run the following from the command line::

    pyrex --input_gen=energy

pyrex will then ask you the name of your coordinate file and generate your standard .json input. You can edit the input file from there to run your desired level of theory and other options. For example, the following input was created from an IRC file for a reaction of carbon dioxide and the hydrogen molecule::

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
        "nthreads": 1,
        "irc_filename": "full_irc.xyz",
        "do_energy": true,
        "do_conceptualdft": true,
        "irc_stepsize": 0.2
      }
    }
