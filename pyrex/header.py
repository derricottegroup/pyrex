import os
import datetime

def header(outfile,json_data):
    output = open(outfile, 'w+')
    header ='''
-----------------------------------------------------------------------
     PYREX: Python Reaction Energy eXtension for Quantum Chemistry

                        PYREX version 1.0

                      Wallace D. Derricotte
                    Derricotte Research Group
                        Morehouse College
                         Atlanta,Georgia
-----------------------------------------------------------------------
'''
    pid = os.getpid()
    output.write(header)
    datetime_now = datetime.datetime.now().strftime('pyREX Run on %A %B %d, %Y at %H:%M%p')
    output.write("\n")
    output.write(datetime_now)
    output.write("\n")
    output.write("Process ID: %d" %pid)

    # Write input
    output.write("\n\n======================> JSON Input File <======================\n\n")
    output.write(json_data)
    output.write("\n\n===============================================================\n\n")
    output.close()
