"""Various miscellaneous functions"""

__author__ = "Adam Simpkin, Felix Simkovic & Jens Thomas"
__date__ = "05 May 2017"
__version__ = "1.0"

from distutils.version import StrictVersion

import logging
import os
import shlex
import subprocess
import sys
import tempfile
import warnings

CCP4_VERSION = None
EXE_EXT = '.exe' if sys.platform.startswith('win') else ''
SCRIPT_HEADER = '' if sys.platform.startswith('win') else '#!/bin/bash'

logger = logging.getLogger(__name__)


# Migh consider moving this to simbad/command_line/__init__.py
def ccp4_version():
    """Identify the CCP4 version

    Returns
    -------
    obj
       A :obj:`StrictVersion` object containing the CCP4 version
    
    Raises
    ------
    RuntimeError
       Cannot determine CCP4 version
    Run
    """
    if CCP4_VERSION is None:
        # Currently there seems no sensible way of doing this other then running a program and grepping the output
        pdbcur = 'pdbcur' + EXE_EXT
        log_fname = tmp_file_name(delete=False)
        run_job([pdbcur], logfile=log_fname)
        tversion = None
        for line in open(log_fname, 'r'):
            if line.startswith(' ### CCP4'):
                tversion = line.split()[2].rstrip(':')
        if not tversion:
            raise RuntimeError("Cannot determine CCP4 version")
    os.unlink(log_fname)
    # Create the version as StrictVersion to make sure it's valid and allow for easier comparisons
    return StrictVersion(tversion)


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

    p.stdout.close()
    return p.returncode


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

