Building a Coordinate
=====================

The first step to doing a calculation in PYREX is to obtain some sort of coordinate to do your calculation on. This is typically either an intrinsic reaction coordinate (IRC) or a surface scan of some kind. PYREX is general enough to take any IRC or surface scan you may have calculated in your favorite quantum chemistry software. So if you already have a coordinate and are ready to calculate properties, skip ahead to see how. However, if you need some help building a coordinate PYREX has a few tools that may be useful to get started.

Intrinsic Reaction Coordinate
-----------------------------

PYREX uses the mass-weighted cartesian coordinate definition of the IRC originally proposed by Fukui, the algorithm implemented is essentially a series of steepest-descent steps taken in the "forward" and "backward" direction downhill from the transition state. The algorithm is an implementation of the method proposed by Ishida and Morokuma which basically supplements the initial steepest-descent with a line search algortihm based on a second gradient calculation to ensure that the minimum energy path is followed. Like most steepest-descent algorithms, this one struggles with flat potential energy surfaces. However for cases where the imaginary vibrational frequency is rather large (about 1000 cm:math:`^{-1}` or higher) this algorithm should converge to the appropriate minimum energy path. A small step size of about 0.01 a.u. is recommended for smooth convergence and to obtain a sufficiently tight grid for calculations of numerical derivative properties like the reaction force. An example IRC input is shown below::

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
  
In the "molecule" block the transition state structure is provided in cartesian coordinates. The "irc" block is used to specify irc specific options. The first option is the direction of the irc, whether you will go forward (positive initial gradient) or backward (negative initial gradient). The "step_size" option specifies the step size in atomic units. The "mode" options specifies which mode in the normal mode file you want to follow (this is typically the first mode), and finally the file that contains normal modes from a vibrational frequency calculation are provided. 

Surface Scans
-------------
