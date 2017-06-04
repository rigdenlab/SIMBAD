#!/usr/bin/env ccp4-python
"""Module to run phaser on a model"""

__author__ = "Adam Simpkin"
__date__ = "02 May 2017"
__version__ = "1.0"

import os
import shutil

from pyjob.dispatch import cexec


class Phaser(object):
    """Class to run PHASER

    Attributes
    ----------
    enant : bool
        Specify whether to try enantimorphic space groups
    f : str
        The column label for F
    hklin : str
        Path to the input hkl file
    hklout : str
        Path to the output hkl file
    logfile : str
        Path to the output log file
    pdbin : str
        Path to the input pdb file
    pdbout : str
        Path to the output pdb file
    sigf : str
        The column label for SIGF
    solvent : int float
        The estimated solvent content of the crystal
    work_dir : str
        Path to the working directory were you want PHASER to run

    Examples
    --------
    >>> from simbad.mr.phaser_mr import Phaser
    >>> phaser = Phaser('<enant>', '<f>', '<hklin>', '<hklout>', '<logfile>', '<pdbin>', '<pdbout>', '<sigf>',
    >>>                 '<solvent>', '<workdir>')
    >>> phaser.run()

    Files relating to the PHASER run will be contained within the work_dir however the location of the output hkl, pdb
    and logfile can be specified.
    """

    def __init__(self, enant, f, hklin, hklout, logfile, pdbin, pdbout, sigf, solvent, work_dir):
        self._enant = None
        self._f = None
        self._hklin = None
        self._hklout = None
        self._logfile = None
        self._pdbout = None
        self._pdbout = None
        self._sigf = None
        self._solvent = None
        self._work_dir = None

        self.enant = enant
        self.f = f
        self.hklin = hklin
        self.hklout = hklout
        self.logfile = logfile
        self.pdbin = pdbin
        self.pdbout = pdbout
        self.sigf = sigf
        self.solvent = solvent
        self.work_dir = work_dir

    @property
    def enant(self):
        """Whether to check for enantimophic space groups"""
        return self._enant

    @enant.setter
    def enant(self, enant):
        """Define whether to check for enantiomorphic space groups"""
        self._enant = enant

    @property
    def f(self):
        """The F label from the input hkl"""
        return self._f

    @f.setter
    def f(self, f):
        """Define the F label from the input hkl"""
        self._f = f

    @property
    def hklin(self):
        """The input hkl file"""
        return self._hklin

    @hklin.setter
    def hklin(self, hklin):
        """Define the input hkl file"""
        self._hklin = hklin

    @property
    def hklout(self):
        """The output hkl file"""
        return self._hklout

    @hklout.setter
    def hklout(self, hklout):
        """Define the output hkl file"""
        self._hklout = hklout

    @property
    def logfile(self):
        """The logfile output"""
        return self._logfile

    @logfile.setter
    def logfile(self, logfile):
        """Define the output logfile"""
        self._logfile = logfile

    @property
    def pdbin(self):
        """The input pdb file"""
        return self._pdbin

    @pdbin.setter
    def pdbin(self, pdbin):
        """Define the input pdb file"""
        self._pdbin = pdbin

    @property
    def pdbout(self):
        """The output pdb file"""
        return self._pdbout

    @pdbout.setter
    def pdbout(self, pdbout):
        """Define the output pdb file"""
        self._pdbout = pdbout

    @property
    def sigf(self):
        """The SIGF label from the input hkl"""
        return self._sigf

    @sigf.setter
    def sigf(self, sigf):
        """Define the SIGF label from the input hkl"""
        self._sigf = sigf

    @property
    def solvent(self):
        """The estimated solvent content of the crystal"""
        return self._solvent

    @solvent.setter
    def solvent(self, solvent):
        """Define the estimated solvent content of the crystal"""
        self._solvent = solvent

    def run(self):
        """Function to run molecular replacement using PHASER

        Returns
        -------
        file
            Output hkl file
        file
            Output pdb file
        file
            Output log file
        """
        
        # Make a note of the current working directory
        current_work_dir = os.getcwd()

        # Change to the PHASER working directory
        if os.path.exists(self.work_dir):
            os.chdir(self.work_dir)
        else:
            os.makedirs(self.work_dir)
            os.chdir(self.work_dir)

        # Copy hklin and pdbin to working dire for efficient running of PHASER
        hklin = os.path.join(self.work_dir, os.path.basename(self.hklin))
        shutil.copyfile(self.hklin, hklin)
        pdbin = os.path.join(self.work_dir, os.path.basename(self.pdbin))
        shutil.copyfile(self.pdbin, pdbin)

        if self.enant:
            sgalternative = "ALL"
        else:
            sgalternative = "HAND"

        key = """#---PHASER COMMAND SCRIPT GENERATED BY SIMBAD---
        MODE MR_AUTO
        ROOT "phaser_mr_output"'
        #---DEFINE DATA---
        HKLIN {0}
        LABIN F={1} SIGF={2}
        SGALTERNATIVE SELECT {3}
        #---DEFINE ENSEMBLES---
        ENSEMBLE ensemble1 &
            PDB "{4}" RMS 0.6
        #---DEFINE COMPOSITION---
        COMPOSITION BY SOLVENT
        COMPOSITION PERCENTAGE {5}
        #---SEARCH PARAMETERS---
        SEARCH ENSEMBLE ensemble1 NUMBER 1"""
        key = key.format(hklin, self.f, self.sigf, sgalternative, pdbin, self.solvent)

        Phaser.phaser(self.logfile, key)

        # Move the output hkl file to specified filename
        if os.path.isfile(os.path.join(self.work_dir, 'phaser_mr_output.1.mtz')):
            shutil.move(os.path.join(self.work_dir, 'phaser_mr_output.1.mtz'), self.hklout)
        # Move output pdb file to specified filename
        if os.path.isfile(os.path.join(self.work_dir, 'phaser_mr_output.1.pdb')):
            shutil.move(os.path.join(self.work_dir, 'phaser_mr_output.1.pdb'), self.pdbout)

        # Return to original working directory
        os.chdir(current_work_dir)

        # Delete any files copied across
        if os.path.isfile(os.path.join(self.work_dir, os.path.basename(self.hklin))):
            os.remove(os.path.join(self.work_dir, os.path.basename(self.hklin)))
        if os.path.isfile(os.path.join(self.work_dir, os.path.basename(self.pdbin))):
            os.remove(os.path.join(self.work_dir, os.path.basename(self.pdbin)))

    @staticmethod
    def phaser(logfile, key):
        """Function to run molecular replacement using PHASER

        Parameters
        ----------
        logfile : str
            Path to the output log file
        key : str
            PHASER keywords

        Returns
        -------
        file
            Output hkl file
        file
            Output pdb file
        file
            Output log file
        """
        cmd = ["phaser"]
        stdout = cexec(cmd, stdin=key)
        with open(logfile, 'w') as f_out:
            f_out.write(stdout)
       

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description='Runs MR using PHASER', prefix_chars="-")

    group = parser.add_argument_group()
    group.add_argument('-enant', type=bool,
                       help="Try enantimorph space groups <True/False>")
    group.add_argument('-f', type=str,
                       help="The column label for F")
    group.add_argument('-hklin', type=str,
                       help="Path the input hkl file")
    group.add_argument('-hklout', type=str,
                       help="Path the output hkl file")
    group.add_argument('-logfile', type=str,
                       help="Path to the ouput log file")
    group.add_argument('-pdbin', type=str,
                       help="Path to the input pdb file")
    group.add_argument('-pdbout', type=str,
                       help="Path to the output pdb file")
    group.add_argument('-sigf', type=str,
                       help="The column label for SIGF")
    group.add_argument('-solvent',
                       help="The estimated solvent content of the crystal")
    group.add_argument('-work_dir', type=str,
                       help="Path to the working directory")
    args = parser.parse_args()

    phaser = Phaser(args.enant, args.f, args.hklin, args.hklout, args.logfile, 
                    args.pdbin, args.pdbout, args.sigf, args.solvent, args.work_dir)
    phaser.run()

