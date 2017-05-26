"""Various miscellaneous functions"""

__author__ = "Adam Simpkin, Felix Simkovic & Jens Thomas"
__date__ = "05 May 2017"
__version__ = "1.0"

import logging
import os
import string
import subprocess
import sys

import iotbx.pdb
import iotbx.pdb.amino_acid_codes

import mbkit.chemistry
import mbkit.dispatch.cexectools

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


def molecular_weight(pdbin):
    """Function to calculate the molecular weight of a model

    Parameters
    ----------
    pdbin : str
       Path to input pdb

    Returns
    -------
    float
       Molecular weight of input model
    """

    pdb_input = iotbx.pdb.pdb_input(file_name=pdbin)
    hierarchy = pdb_input.construct_hierarchy()

    # Define storage variables
    mw = 0
    hydrogen_atoms = 0

    # Collect all the data using the hierarchy
    for m in hierarchy.models():
        for c in m.chains():
            for rg in c.residue_groups():
                resseq = None
                for ag in rg.atom_groups():
                    if ag.resname in iotbx.pdb.amino_acid_codes.one_letter_given_three_letter and resseq != rg.resseq:
                        resseq = rg.resseq
                        try:
                            hydrogen_atoms += mbkit.chemistry.atomic_composition[ag.resname].H
                        except AttributeError:
                            logger.debug('Ignoring non-standard amino acid: %s', ag.resname)
                    for atom in ag.atoms():
                        if ag.resname.strip() == 'HOH' or ag.resname.strip() == 'WAT':
                            # Ignore water atoms
                            pass
                        else:
                            # Be careful, models might not have the last element column
                            if atom.element.strip():
                                aname = atom.element.strip()
                            else:
                                aname = atom.name.strip()
                                aname = aname.translate(None, string.digits)[0]
                            try:
                                mw += mbkit.chemistry.periodic_table[aname].atomic_mass * atom.occ
                            except AttributeError:
                                try:
                                    aname = ''.join([i for i in aname if not i.isdigit()])
                                    mw += mbkit.chemistry.periodic_table[aname].atomic_mass * atom.occ
                                except AttributeError:
                                    logger.debug('Ignoring non-standard atom type: %s', aname)

    mw += hydrogen_atoms * mbkit.chemistry.periodic_table['H'].atomic_mass

    return mw

