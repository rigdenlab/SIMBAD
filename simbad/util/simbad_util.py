"""Various miscellaneous functions"""

__author__ = "Adam Simpkin, Felix Simkovic & Jens Thomas"
__date__ = "05 May 2017"
__version__ = "1.0"

import logging
import os

import mbkit.dispatch.cexectools

logger = logging.getLogger(__name__)


def filename_append(filename=None, astr=None, directory=None, separator="_"):
    """Append astr to filename, before the suffix, and return the new filename."""
    dirname, fname = os.path.split(filename)
    name, suffix = os.path.splitext(fname)
    name = name + separator + astr + suffix
    if directory is None:
        directory = dirname
    return os.path.join(directory, name)


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
    stdout = mbkit.dispatch.cexectools.cexec(cmd)
    # Exctract molecular weight from log file
    molecular_weight = None
    for line in stdout.split(os.linesep):
        if line.startswith(" Molecular Weight of protein"):
            molecular_weight = float(line.split()[-1])
    if molecular_weight is None:
        msg = "Cannot find Molecular weight in logfile {0}".format(logfile)
        logger.debug(msg)
        raise RuntimeError(msg)
    return molecular_weight

