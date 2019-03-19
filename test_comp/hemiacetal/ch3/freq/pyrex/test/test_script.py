import os
import subprocess


#subprocess.call(['java', '-jar', '/root/wderricotte/src/janpa/molden2molden.jar', '-NormalizeBF', '-cart2pure', '-i', 'r.molden', '-o', 'r_conv.molden'])
log_file = open('janpa.log', 'w+')
subprocess.call(['java', '-jar', '/root/wderricotte/src/janpa/janpa.jar', '-i', 'r_conv.molden', '-CLPO_Molden_File', 'CLPO.molden'], stdout=log_file)

