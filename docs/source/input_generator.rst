Input Generator
===============

In order to get started using pyREX you must have a .xyz coordinate file that contains the geometries for the reaction/process you are investigating and a .json input file that contains the relevant options for the calculation you would like to do in pyREX. There is a built-in input file generator that can help you with all the steps of getting your input file set up.

Joining IRC Files
-----------------
Often times the coordinate file you would like to run in pyREX is an intrinsic reaction coordinate (IRC) file. Most IRC codes output separate coordinates for the forward and backward trajectories for the reaction. However for calculations in pyrex you will need a continuous trajectory that begins with the last structure of your backward trajectory and ends with the last structure of your forward trajectory. You can do this very simply from the command line using the following command line utility::
    pyrex --join_irc
pyrex will then ask you the name of your forward and backward trajectory and produce a single "full_irc.xyz" file. 
