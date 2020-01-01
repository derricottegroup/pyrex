Input Generator
===============

In order to get started using pyREX you must have a .xyz coordinate file that contains the geometries for the reaction/process you are investigating and a .json input file that contains the relevant options for the calculation you would like to do in pyREX. There is a built-in input file generator that can help you with all the steps of getting your input file set up.


Creating an IRC Input
---------------------

Pyrex contains an IRC algorithm that uses gradients from QM software in order to calculate a minimum energy path for the reaction. If you have an optimized transition state structure you can generate an input file to run an IRC calculation using the following command line utility::

    pyrex --input_gen=irc 

pyrex will then ask you the name of your XYZ coordinate file containing your transition state structure. You can edit the input file from there to run the IRC calculation at your preferred level of theory. For example the following input file was created from an xyz file containing the transition state for the reaction of carbon dioxide and the hydrogen molecule::

    {
    "molecule": {
      "geometry" :  [0.113595072843,-0.153880407652,0.037596022684,
                    -1.135961394234,-0.223379947195,-0.034588131338,
                    1.064242147389,0.524605985828,-0.028439251841,
                    0.404693854873,-1.556634791353,0.327550495378,
                    -0.619014433727,-1.391813384411,0.225090412236],
      "symbols": ["C","O","O","H","H"],
      "molecular_charge": "0",
      "molecular_multiplicity": 1
    },
      "model": {
      "method": "scf",
      "basis": "sto-3g"
    },
      "pyrex": {
      "nthreads": 1
    },
      "irc": {
      "direction": "forward",
      "step_size": 0.01,
      "mode": 1,
      "normal_mode_file": "normal_modes.dat"
      }
    }

where the additional options for the irc calculation are added in the IRC block. 

.. caution::
   The IRC algorithm implemented is a simple version of the Morokuma-Ishida predictor-corrector IRC. It is suggested that you use a step size of 0.01 au for all IRC calculations or else the algorithm becomes unstable. This fine of a stepsize is also recommended for computing finite-differences along the reaction coordinate. If you wish to use the IRC to simply confirm a transition state and want to converge as quickly as possible we recommend using other algorithms suited for this purpose. 

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
