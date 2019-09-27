import os
import time

_log_output = False
def log(msg):
    global _log_output
    if _log_output:
        print(msg)

orca_path = "orca"
default_output_root_dir = os.path.dirname(os.path.realpath(__file__))+os.sep+"output"+os.sep
def run_orca(xyzfile=None, xyzstring=None, jobname=None, orca_input=None,
    overwrite=False, output_root_dir=None):
    global orca_path
    global default_output_root_dir
    if output_root_dir is None:
        output_root_dir = default_output_root_dir
    if jobname is None:
        jobname_suggestion = "job_"+str(time.strftime("%d_%m_%Y"))
        jobname = jobname_suggestion
        jobname_counter = 2
        while os.path.exists(os.path.join(output_root_dir, jobname) + os.sep):
            jobname = jobname_suggestion + "_j" + str(jobname_counter)
            jobname_counter += 1
    assert(orca_input is not None)
    explicit_xyz = (xyzfile is not None or xyzstring is not None)

    output_dir = os.path.join(output_root_dir, jobname) + os.sep
    if os.path.isdir(output_dir) and (not overwrite):
        log("output directory {} already exists. Doing nothing.".format(output_dir)+
            " Use overwrite=True to overwrite")
        return None
    log("orca executable path is {}".format(orca_path))
    if orca_path == "orca":
        log("assuming orca is on your path")
    log("output directory is {}".format(output_dir))
    if not os.path.exists(output_dir):
        log("creating output directory")
        os.makedirs(output_dir)
    log("cleaning output directory from previous job files")
    previous_wd = os.getcwd()
    os.chdir(output_dir)
    #os.system("rm job*")
    os.chdir(previous_wd)
    if xyzfile:
        with open(xyzfile, 'r') as fin:
            xyzstring = fin.read()
    if explicit_xyz:
        xyzstring += '\n\n'
        xyzinputfile = output_dir + "guess.xyz"
        log("creating xyz geometry input file")
        with open(xyzinputfile, 'w') as fout:
            fout.write(xyzstring)
    else:
        log("neither xyz string nor xyz file specified!")
        log("assuming xyz coordinates are explicitly specified")
        log("in the orca input file ...")
    orca_input += '\n\n'
    jobinputfile = output_dir + "{}.in".format(jobname)
    log("creating job file")
    with open(jobinputfile, 'w') as fout:
        fout.write(orca_input)
    joboutputfile = output_dir + "{}.out".format(jobname)
    log("running job (this could take a while)")
    previous_wd = os.getcwd()
    os.chdir(output_dir)
    #print(joboutputfile)
    os.system("{} {} > {}".format(orca_path, jobinputfile, joboutputfile))
    log("running of job finished")
    log("changing back current working directory")
    os.chdir(previous_wd)
    return ORCAReporter(joboutputfile)


def reporter_by_name(jobname, output_root_dir=None):
    global default_output_root_dir
    if output_root_dir is None:
        output_root_dir = default_output_root_dir
    flepath = os.path.join(os.path.join(output_root_dir, jobname), jobname+".out")
    return ORCAReporter(flepath)


class ORCAReporter(object):
    def __init__(self, joboutputfile):
        self.joboutputfile = joboutputfile
        self._output = None
        self._trajoutput = None

    def output_lines(self):
        lst = []
        for lne in self.output().split('\n'):
            lst.append(lne)
        return lst

    def output(self):
        if self._output is not None:
            return self._output
        else:
            self._output = None
            with open(self.joboutputfile, 'r') as fin:
                self._output = fin.read()
        assert self._output is not None, 'error reading in file {}'.format(self.joboutputfile)
        return self._output























#
