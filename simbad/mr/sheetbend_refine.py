#!/usr/bin/env ccp4-python
"""Module to run sheetbend on a model"""

__author__ = "Adam Simpkin"
__date__ = "05 Aug 2018"
__version__ = "1.0"

import os

from simbad.util import mtz_util
from simbad.mr.refmac_refine import Refmac
from pyjob import cexec


class SheetBend(object):
    """Class to run sheetbend"""

    def __init__(self, hklin, hklout, logfile, pdbin, pdbout, work_dir):
        self._hklin = None
        self._hklout = None
        self._logfile = None
        self._pdbout = None
        self._pdbout = None
        self._work_dir = None

        # Temporary path for testing
        self.exe = "/data1/opt/devtoolsTrunk/install/bin/csheetbend"
        self.hklin = hklin
        self.hklout = hklout
        self.logfile = logfile
        self.pdbin = pdbin
        self.pdbout = pdbout
        self.work_dir = work_dir

        self.check_sheetbend_exe()

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

    def check_sheetbend_exe(self):
        if not os.path.isfile(self.exe):
            msg = "Sheetbend executable {0} not found".format(self.exe)
            raise RuntimeError(msg)

    def run(self, ncyc=100):

        # Make a note of the current working directory
        current_work_dir = os.getcwd()

        # Change to the sheetbend working directory
        if os.path.exists(self.work_dir):
            os.chdir(self.work_dir)
        else:
            os.makedirs(self.work_dir)
            os.chdir(self.work_dir)

        tmp_pdb = os.path.join(self.work_dir, 'sheetbend.pdb')
        SheetBend.sheetbend(self.exe, self.hklin, self.pdbin, tmp_pdb, ncyc, self.logfile)

        # Perform a cycle of Refmac to get output hkl
        key = "ncyc 0"
        Refmac.refmac(self.hklin, self.hklout, tmp_pdb, self.pdbout, self.logfile, key)

        # Return to original working directory
        os.chdir(current_work_dir)

    @staticmethod
    def sheetbend(exe, hklin, pdbin, pdbout, ncyc, logfile):
        """Function to run refinement using sheetbend

        Parameters
        ----------
        hklin : str
            Path to the input hkl file
        pdbin : str
            Path to the input pdb
        pdbout : str
            Path to the output pdb
        ncyc : int
            Number of cycles to run
        logfile : str
            Path to the output log

        Returns
        -------
        file
            Output pdb file
        file
            Output log file
        """

        mtz_labels = mtz_util.GetLabels(hklin)
        colin = "{0},{1}".format(mtz_labels.f, mtz_labels.sigf)

        cmd = [exe, '--pdbin', pdbin, '--mtzin', hklin, '--pdbout',  pdbout, '--colin-fo', colin,
               '-cycles', str(ncyc), '-resolution-by-cycle', '6,6,3']
        stdout = cexec(cmd)
        with open(logfile, 'w') as f_out:
            f_out.write(stdout)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Runs refinement using sheetbend', prefix_chars="-")

    group = parser.add_argument_group()
    group.add_argument('-hklin', type=str,
                       help="Path the input hkl file")
    group.add_argument('-hklout', type=str,
                       help="Path the output hkl file")
    group.add_argument('-logfile', type=str,
                       help="Path to the ouput log file")
    group.add_argument('-ncyc', type=int, default=100,
                       help="Number of cycles of refinement to run")
    group.add_argument('-pdbin', type=str,
                       help="Path to the input pdb file")
    group.add_argument('-pdbout', type=str,
                       help="Path to the output pdb file")
    group.add_argument('-work_dir', type=str,
                       help="Path to the working directory")
    args = parser.parse_args()

    sheetbend = SheetBend(args.hklin, args.hklout, args.logfile, args.pdbin, args.pdbout, args.work_dir)
    sheetbend.run(args.ncyc)