"""Various miscellaneous functions"""

__author__ = "Adam Simpkin, Felix Simkovic & Jens Thomas"
__date__ = "05 May 2017"
__version__ = "1.0"

import logging
import os
import subprocess
import sys
import tempfile

SCRIPT_EXT = '.bat' if sys.platform.startswith('win') else '.sh'
EXE_EXT = '.exe' if sys.platform.startswith('win') else ''
SCRIPT_HEADER = '' if sys.platform.startswith('win') else '#!/bin/bash'

logger = logging.getLogger(__name__)


def filename_append(filename=None, astr=None, directory=None, separator="_"):
    """Append astr to filename, before the suffix, and return the new filename."""
    dirname, fname = os.path.split(filename)
    name, suffix = os.path.splitext(fname)
    name = name + separator + astr + suffix
    if directory is None:
        directory = dirname
    return os.path.join(directory, name)


def run_job(cmd, logfile=None, directory=None, stdin=None):
    """Execute a command and return the exit code.

    Parameters
    ----------
    cmd : list
       Command to run as a list
    stdin : str, optional
       Stdin for the command
    logfile : str, optional
       The path to the logfile
    directory : str, optional
       The directory to run the job in (cwd assumed)
    stdin : str
       Additional keys to go into STDIN 

    Returns
    -------
    returncode : int
       Subprocess exit code

    Notes
    -----
    We take care of outputting stuff to the logs and opening/closing logfiles

    """
    if not directory:
        directory = os.getcwd()
    
    logger.debug("Running job in %s\n%s\n%s", directory, " ".join(cmd), stdin)

    kwargs = {"bufsize":0, "shell":"False"} if os.name == "nt" else {}
    p = subprocess.Popen(
        cmd, stdin=subprocess.PIPE, cwd=directory, stdout=subprocess.PIPE, 
        stderr=subprocess.STDOUT, **kwargs
    )

    # Write the keyword input
    if stdin is not None:
        p.stdin.write(stdin)
        p.stdin.close()

    # Watch the output for successful termination
    if logfile is None:
        logfile = os.devnull

    with open(logfile, "w") as f_out:
        for line in p.stdout:
            f_out.write(line)

    p.wait()
    p.stdout.close()

    return p.returncode


def molecular_weight(model):
    """Function to run ``rwcontents`` to get the molecular weight of a model

    Parameters
    ----------
    model : str
       Path to input model

    Returns
    -------
    float
       Molecular weight of input model

    """
    cmd = ['rwcontents', 'xyzin', model]
    logfile = 'rwcontents_{0}.log'.format(os.path.basename(model).rsplit('.', 1)[0])
    run_job(cmd, logfile=logfile, stdin="")

    # Exctract molecular weight from log file
    molecular_weight = None
    with open(logfile, 'r') as f:
        for line in f:
            if line.startswith(" Molecular Weight of protein"):
                molecular_weight = float(line.split()[-1])
    if molecular_weight is None:
        msg = "Cannot find Molecular weight in logfile {0}".format(logfile)
        logger.debug(msg)
        raise RuntimeError(msg)

    os.remove(logfile)
    return molecular_weight


def tmp_file_name(delete=True, directory=None, prefix="", suffix=""):
    """Return a filename for a temporary file

    Parameters
    ----------
    delete : bool, optional
       Flag whether the temporary file should be deleted [default: True]
    directory : str, optional
       Path to a directory to write the files to.
    prefix : str, optional
       A prefix to the temporary filename
    suffix : str, optional
       A suffix to the temporary filename

    """
    directory = os.getcwd() if not directory else directory
    return tempfile.NamedTemporaryFile(dir=directory, delete=delete, prefix=prefix, suffix=suffix).name

