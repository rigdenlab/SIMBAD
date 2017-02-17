'''
Various miscellaneous functions
Largely stolen from ample
'''

import logging
import os
import shlex
import subprocess
import sys
import tempfile
import warnings

CCP4_VERSION=None
EXE_EXT = '.exe' if sys.platform.startswith('win') else ''
SCRIPT_HEADER = '' if sys.platform.startswith('win') else '#!/bin/bash'


_logger = logging.getLogger(__name__)

def ccp4_version():
    """Return the CCP4 version as a tuple"""
    global CCP4_VERSION
    if CCP4_VERSION is None:
        # Currently there seems no sensible way of doing this other then running a program and grepping the output
        pdbcur = 'pdbcur' + EXE_EXT
        log_fname = tmp_file_name(delete=False)
        run_command([pdbcur], stdin="", logfile=log_fname)
        tversion = None

        with open(log_fname, 'r') as logfh:
            for i, line in enumerate(logfh.readlines()):
                if i > 20:
                    break
                if line.startswith(' ### CCP4'):
                    tversion = line.split()[2].rstrip(':')
                    break

        if not tversion:
            raise RuntimeError("Cannot determine CCP4 version")
        vsplit = tversion.split('.')
        if len(vsplit) == 2:
            major = int(vsplit[0])
            minor = int(vsplit[1])
            rev = '-1'
        elif len(vsplit) == 3:
            major = int(vsplit[0])
            minor = int(vsplit[1])
            rev = int(vsplit[2])
        else:
            raise RuntimeError("Cannot split CCP4 version: {0}".format(tversion))
    os.unlink(log_fname)
    return (major, minor, rev)

def filename_append(filename=None, astr=None,directory=None, separator="_"):
    """Append astr to filename, before the suffix, and return the new filename."""
    dirname, fname = os.path.split( filename )
    name, suffix = os.path.splitext( fname )
    name  =  name + separator + astr + suffix
    if directory is None:
        directory = dirname
    return os.path.join( directory, name )

def is_exe(fpath):
    return fpath and os.path.exists(fpath) and os.access(fpath, os.X_OK)

def run_command(cmd, logfile=None, directory=None, dolog=True, stdin=None, check=False, env=None):
    """Execute a command and return the exit code.

    We take care of outputting stuff to the logs and opening/closing logfiles

    Args:
    cmd - command to run as a list
    stdin - a string to use as stdin for the command
    logfile (optional) - the path to the logfile
    directory (optional) - the directory to run the job in (cwd assumed)
    dolog: bool - whether to output info to the system log
    """
    assert type(cmd) is list, "run_command needs a list!"
    if check:
        if not is_exe(cmd[0]): raise RuntimeError,"run_command cannot find executable: {0}".format(cmd[0])

    if not directory:  directory = os.getcwd()
    if dolog: _logger.debug("In directory {0}\nRunning command: {1}".format(directory, " ".join(cmd)))
    file_handle=False
    if logfile:
        if type(logfile)==file:
            file_handle=True
            logf=logfile
            logfile=os.path.abspath(logf.name)
        else:
            logfile = os.path.abspath(logfile)
            logf = open(logfile, "w")
        if dolog: _logger.debug("Logfile is: {0}".format(logfile))
    else:
        logf = tempfile.TemporaryFile()
        
    if stdin != None:
        stdinstr = stdin
        stdin = subprocess.PIPE

    # Windows needs some special treatment
    kwargs = {}
    if os.name == "nt":
        kwargs = { 'bufsize': 0, 'shell' : "False" }
    p = subprocess.Popen(cmd, stdin=stdin, stdout=logf, stderr=subprocess.STDOUT, cwd=directory, env=env, **kwargs)

    if stdin != None:
        p.stdin.write( stdinstr )
        p.stdin.close()
        if dolog: _logger.debug("stdin for cmd was: {0}".format( stdinstr ) )

    p.wait()
    if not file_handle: logf.close()
    return p.returncode
def tmpFileName():
    """Return a filename for a temporary file

    See Also
    --------
    tmp_file_name

    Warnings
    --------
    This function was deprecated and will be removed in future releases. Please use ``tmp_file_name()`` instead.

    """
    msg = "This function was deprecated and will be removed in future release"
    warnings.warn(msg, DeprecationWarning, stacklevel=2)
    return tmp_file_name()
def run_job(command_line, logfile, key=""):
    """ Generic job runner """

    if os.name == "nt":
        process_args = shlex.split(command_line, posix=False)
        p = subprocess.Popen(process_args, bufsize=0, shell="False", stdin = subprocess.PIPE,
                             stdout = subprocess.PIPE, stderr = subprocess.STDOUT)
    else:
        process_args = shlex.split(command_line)
        p = subprocess.Popen(process_args, stdin = subprocess.PIPE,
                             stdout = subprocess.PIPE, stderr = subprocess.STDOUT)


    (child_stdin, child_stdout_and_stderr) = (p.stdin, p.stdout)

    # Write the keyword input
    child_stdin.write(key)
    child_stdin.close()

    # Watch the output for successful termination
    out=child_stdout_and_stderr.readline()

    log=open(logfile, "w")

    while out:
        out=child_stdout_and_stderr.readline()
        log.write(out)

    child_stdout_and_stderr.close()
    log.close()

    return

def tmpFileName():
    """Return a filename for a temporary file

    See Also
    --------
    tmp_file_name

    Warnings
    --------
    This function was deprecated and will be removed in future releases. Please use ``tmp_file_name()`` instead.

    """
    msg = "This function was deprecated and will be removed in future release"
    warnings.warn(msg, DeprecationWarning, stacklevel=2)
    return tmp_file_name()

def tmp_file_name(delete=True, directory=None, suffix=""):
    """Return a filename for a temporary file

    Parameters
    ----------
    delete : bool, optional
       Flag whether the temporary file should be deleted [default: True]
    directory : str, optional
       Path to a directory to write the files to.
    suffix : str, optional
       A suffix to the temporary filename

    """
    directory = os.getcwd() if not directory else directory
    t = tempfile.NamedTemporaryFile(dir=directory, delete=delete, suffix=suffix)
    tmp1 = t.name
    t.close()
    return tmp1

def generate_mr_input_file(optd, model, stage):
    '''Create an input file for MR'''

    # create input file path
    if stage == 'lattice':
        input_file = os.path.join(optd.d['work_dir'], 'MR_LATTICE', model, 'input.txt')
    elif stage == 'contaminant':
        input_file = os.path.join(optd.d['work_dir'], 'MR_CONTAMINANT', model, 'input.txt')
    elif stage == 'full':
        input_file = os.path.join(optd.d['work_dir'], 'MR_FULL', model, 'input.txt')

    # Assign variables
    DIRE = optd.d['work_dir']
    SGIN = optd.d['space_group']
    HKL1 = optd.d['mtz']
    HKLR = '{0}_refinement_input.mtz'.format(model)
    HKLO = '{0}_phaser_output.mtz'.format(model)
    PDBO = '{0}_mr_output.pdb'.format(model)
    MRLO = '{0}_mr.log'.format(model)
    REFH = '{0}_refinement_output.mtz'.format(model)
    REFP = '{0}_refinement_output.pdb'.format(model)
    ENAN = optd.d['enan']
    FPIN = optd.d['F']
    SIGF = optd.d['SIGF']
    FREE = optd.d['FREE']
    SOLV = optd.d['solvent']
    RESO = optd.d['resolution']
    PDBI = os.path.join(optd.d['work_dir'], 'lattice_input_models', model)


    # Write input file
    with open(input_file, 'w') as f:
        f.write("DIRE {0}\n".format(DIRE))
        f.write("SGIN {0}\n".format(SGIN))
        f.write("HKL1 {0}\n".format(HKL1))
        f.write("HKLR {0}\n".format(HKLR))
        f.write("PDBO {0}\n".format(PDBO))
        f.write("MRLO {0}\n".format(MRLO))
        f.write("REFH {0}\n".format(REFH))
        f.write("REFP {0}\n".format(REFP))
        f.write("ENAN {0}\n".format(ENAN))
        f.write("FPIN {0}\n".format(FPIN))
        f.write("SIGF {0}\n".format(SIGF))
        f.write("FREE {0}\n".format(FREE))
        f.write("SOLV {0}\n".format(SOLV))
        f.write("RESO {0}\n".format(RESO))
        f.write("PDBI {0}\n".format(PDBI))
        if optd.d['MR_program'] == "phaser":
            f.write("HKLO {0}\n".format(HKLO))

    return


# ======================================================================
# Some default string messages that we need during the program to inform
# the user of certain information
# ======================================================================

header = """
#####################################################################################################
#####################################################################################################
#####################################################################################################
# CCP4: SIMBAD - Sequence Independent Molecular replacement Based on Available Database             #
#####################################################################################################

"""


