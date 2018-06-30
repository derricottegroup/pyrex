import os
import datetime

def header(outfile):
    output = open(outfile, 'w+')
    header ='''
-----------------------------------------------------------------------
     pyREX: Python Reaction Energy eXtension for Quantum Chemistry

                         pyREX beta 0.1

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
    output.close()
