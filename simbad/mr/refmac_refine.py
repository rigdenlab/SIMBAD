#!/usr/bin/env ccp4-python
"""Module to run REFMAC on a model"""

__author__ = "Adam Simpkin"
__date__ = "02 May 2017"
__version__ = "1.0"

import os

from pyjob.dispatch import cexec


class Refmac(object):
    """Class to run refmac

    Attributes
    ----------
    hklin : str
        Path to the input hkl file
    hklout : str
        Path to the output hkl file
    pdbin : str
        Path to the input pdb file
    pdbout : str
        Path to the output pdb file
    logfile : str
        Path to the output logfile
    key : str
        REFMAC key words
    work_dir : str
        Path to the working directory were you want REFMAC to run
    ncyc : int float
        The number of cycles of refinement to perform [default : 30]

    Examples
    --------
    >>> from simbad.mr.refmac_refine import Refmac
    >>> refmac = Refmac('<hklin>', '<hklout>', '<logfile>', '<pdbin>', '<pdbout>', '<work_dir>')
    >>> refmac.run('<ncyc>')

    Files relating to the REFMAC run will be contained within the work_dir however the location of the output hkl, pdb
    and logfile can be specified.
    """

    def __init__(self, hklin, hklout, logfile, pdbin, pdbout, work_dir):
        self._hklin = None
        self._hklout = None
        self._logfile = None
        self._pdbout = None
        self._pdbout = None
        self._work_dir = None

        self.hklin = hklin
        self.hklout = hklout
        self.logfile = logfile
        self.pdbin = pdbin
        self.pdbout = pdbout
        self.work_dir = work_dir

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
    def work_dir(self):
        """The path to the working directory"""
        return self._work_dir

    @work_dir.setter
    def work_dir(self, work_dir):
        """Define the working directory"""
        self._work_dir = work_dir

    def run(self, ncyc=30):
        """Function to run refinement using REFMAC

        Parameters
        ----------
        ncyc : int float
            The number of cycles of refinement to perform [default : 30]

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
        
        # Change to the REFMAC working directory
        if os.path.exists(self.work_dir):
            os.chdir(self.work_dir)
        else:
            os.makedirs(self.work_dir)
            os.chdir(self.work_dir)

        key = "ncyc {0}".format(ncyc)
        Refmac.refmac(self.hklin, self.hklout, self.pdbin, self.pdbout, self.logfile, key)
        
        # Return to original working directory
        os.chdir(current_work_dir)

    @staticmethod
    def refmac(hklin, hklout, pdbin, pdbout, logfile, key):
        """Function to run refinement using REFMAC

        Parameters
        ----------
        hklin : str
            Path to the input hkl file
        hklout : str
            Path to the output hkl file
        pdbin : str
            Path to the input pdb file
        pdbout : str
            Path to the output pdb file
        logfile : str
            Path to the output logfile
        key : str
            REFMAC key words

        Returns
        -------
        file
            Output hkl file
        file
            Output pdb file
        file
            Output log file
        """
        cmd = ['refmac5', 'hklin', hklin, 'hklout', hklout,      
               'xyzin', pdbin, 'xyzout', pdbout]
        stdout = cexec(cmd, stdin=key)
        with open(logfile, 'w') as f_out:
            f_out.write(stdout)


if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description='Runs refinement using REFMAC', prefix_chars="-")

    group = parser.add_argument_group()
    group.add_argument('-hklin', type=str,
                       help="Path the input hkl file")
    group.add_argument('-hklout', type=str,
                       help="Path the output hkl file")
    group.add_argument('-logfile', type=str,
                       help="Path to the ouput log file")
    group.add_argument('-ncyc', type=int,
                       help="Number of cycles of refinement to run")
    group.add_argument('-pdbin', type=str,
                       help="Path to the input pdb file")
    group.add_argument('-pdbout', type=str,
                       help="Path to the output pdb file")
    group.add_argument('-work_dir', type=str,
                       help="Path to the working directory")
    args = parser.parse_args()
    
    refmac = Refmac(args.hklin, args.hklout, args.logfile, args.pdbin, args.pdbout, args.work_dir)
    refmac.run(args.ncyc)
