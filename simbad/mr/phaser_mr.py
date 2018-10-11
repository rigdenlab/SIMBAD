#!/usr/bin/env ccp4-python
"""Module to run phaser on a model"""

__author__ = "Adam Simpkin"
__date__ = "24 March 2018"
__version__ = "1.0"

import os
import shutil

from simbad.mr.options import SGAlternatives
from simbad.util import mtz_util

from phaser import InputMR_DAT, runMR_DAT, InputMR_AUTO, runMR_AUTO


class Phaser(object):
    """Class to run PHASER

    Attributes
    ----------
    hklin : str
        Path to the input hkl file
    f : str
        The column label for F
    i : str
        The column label for I
    logfile : str
        Path to the output log file
    pdbin : str
        Path to the input pdb file
    pdbout : str
        Path to the output pdb file
    sgalternative : str
        Specify whether to try alternative space groups (all | enant)
    sigf : str
        The column label for SIGF
    sigi : str
        The column label for SIGI
    solvent : int float
        The estimated solvent content of the crystal
    work_dir : str
        Path to the working directory were you want PHASER to run
    autohigh : str
        The high resolution limit in Angstroms for final high resolution refinement in MR_AUTO mode
    hires : str
        The high resolution limit of data used to find/refine this solution


    Examples
    --------
    >>> from simbad.mr.phaser_mr import Phaser
    >>> phaser = Phaser('<hklin>', '<hklout>', '<f>', '<i>', '<logfile>', '<nmol>', '<pdbin>', '<pdbout>',
    >>>                 '<sgalternative>', '<sigf>', '<sigi>', '<solvent>', '<timeout>', '<workdir>', '<autohigh>',
    >>>                 '<hires>')
    >>> phaser.run()

    Files relating to the PHASER run will be contained within the work_dir however the location of the output hkl, pdb
    and logfile can be specified.
    """

    def __init__(self, hklin, hklout, f, i, logfile, nmol, pdbin, pdbout, sgalternative, sigf, sigi, solvent, timeout,
                 work_dir, hires, autohigh):
        self._f = None
        self._i = None
        self._autohigh = None
        self._hires = None
        self._hklin = None
        self._logfile = None
        self._nmol = None
        self._pdbin = None
        self._pdbout = None
        self._sgalternative = None
        self._sigf = None
        self._sigi = None
        self._solvent = None
        self._timeout = None
        self._work_dir = None

        self.sgalternative = sgalternative
        self.f = f
        self.i = i
        self.autohigh = autohigh
        self.hires = hires
        self.hklin = hklin
        self.hklout = hklout
        self.logfile = logfile
        self.nmol = nmol
        self.pdbin = pdbin
        self.pdbout = pdbout
        self.sigf = sigf
        self.sigi = sigi
        self.solvent = solvent
        self.timeout = timeout
        self.work_dir = work_dir

    @property
    def f(self):
        """The F label from the input hkl"""
        return self._f

    @f.setter
    def f(self, f):
        """Define the F label from the input hkl"""
        self._f = f

    @property
    def i(self):
        """The I label from the input hkl"""
        return self._i

    @i.setter
    def i(self, i):
        """Define the I label from the input hkl"""
        self._i = i

    @property
    def autohigh(self):
        """The high resolution limit in Angstroms for final high resolution refinement in MR_AUTO mode"""
        return self._autohigh

    @autohigh.setter
    def autohigh(self, autohigh):
        """Define the high resolution limit in Angstroms for final high resolution refinement in MR_AUTO mode"""
        self._autohigh = autohigh

    @property
    def hires(self):
        """The high resolution limit of data used to find/refine this solution"""
        return self._hires

    @hires.setter
    def hires(self, hires):
        """Define the high resolution limit of data used to find/refine this solution"""
        self._hires = hires

    @property
    def hklin(self):
        """The input hkl file"""
        return self._hklin

    @hklin.setter
    def hklin(self, hklin):
        """Define the input hkl file"""
        self._hklin = hklin

    @property
    def logfile(self):
        """The logfile output"""
        return self._logfile

    @logfile.setter
    def logfile(self, logfile):
        """Define the output logfile"""
        self._logfile = logfile

    @property
    def nmol(self):
        """The number of molecules to look for"""
        return self._nmol

    @nmol.setter
    def nmol(self, nmol):
        """Define the number of molecules to look for"""
        self._nmol = nmol

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
    def sgalternative(self):
        """Whether to check for alternative space groups"""
        return self._sgalternative

    @sgalternative.setter
    def sgalternative(self, sgalternative):
        """Define whether to check for alternative space groups"""
        self._sgalternative = sgalternative.lower()

    @property
    def sigf(self):
        """The SIGF label from the input hkl"""
        return self._sigf

    @sigf.setter
    def sigf(self, sigf):
        """Define the SIGF label from the input hkl"""
        self._sigf = sigf

    @property
    def sigi(self):
        """The SIGI label from the input hkl"""
        return self._sigi

    @sigi.setter
    def sigi(self, sigi):
        """Define the SIGI label from the input hkl"""
        self._sigi = sigi

    @property
    def solvent(self):
        """The estimated solvent content of the crystal"""
        return self._solvent

    @solvent.setter
    def solvent(self, solvent):
        """Define the estimated solvent content of the crystal"""
        self._solvent = solvent

    @property
    def timeout(self):
        """The time in minutes before phaser is killed"""
        return self._timeout

    @timeout.setter
    def timeout(self, timeout):
        """Define the time in minutes before phaser should be killed"""
        self._timeout = timeout

    def run(self):
        """Function to run molecular replacement using PHASER

        Returns
        -------
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

        i = InputMR_DAT()
        i.setHKLI(hklin)

        if self.hires:
            i.setHIRES(self.hires)
        if self.autohigh:
            i.setRESO_AUTO_HIGH(self.autohigh)
        if self.i != "None" and self.sigi != "None":
            i.setLABI_I_SIGI(self.i, self.sigi)
        elif self.f != "None" and self.sigf != "None":
            i.setLABI_F_SIGF(self.f, self.sigf)
        else:
            msg = "No flags for intensities or amplitudes have been provided"
            raise RuntimeError(msg)
        i.setSGAL_SELE(SGAlternatives[self.sgalternative].value)
        i.setMUTE(True)
        r = runMR_DAT(i)

        if r.Success():
            i = InputMR_AUTO()
            i.setJOBS(1)
            i.setREFL_DATA(r.getREFL_DATA())
            i.setROOT("phaser_mr_output")
            i.addENSE_PDB_ID("PDB", pdbin, 0.7)
            i.setCOMP_BY("SOLVENT")
            i.setCOMP_PERC(self.solvent)
            i.addSEAR_ENSE_NUM('PDB', self.nmol)
            i.setSGAL_SELE(SGAlternatives[self.sgalternative].value)
            if self.timeout != 0:
                i.setKILL_TIME(self.timeout)
            i.setMUTE(True)
            del(r)
            r = runMR_AUTO(i)

            with open(self.logfile, 'w') as f:
                f.write(r.summary())

            shutil.move(r.getTopPdbFile(), self.pdbout)

            # Output original mtz with a change of basis if needed
            space_group, _, _ = mtz_util.crystal_data(r.getTopMtzFile())
            mtz_util.reindex(self.hklin, self.hklout, space_group)

        # Return to original working directory
        os.chdir(current_work_dir)

        # Delete any files copied across
        if os.path.isfile(os.path.join(self.work_dir, os.path.basename(self.hklin))):
            os.remove(os.path.join(self.work_dir, os.path.basename(self.hklin)))
        if os.path.isfile(os.path.join(self.work_dir, os.path.basename(self.pdbin))):
            os.remove(os.path.join(self.work_dir, os.path.basename(self.pdbin)))

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description='Runs MR using PHASER', prefix_chars="-")

    group = parser.add_argument_group()
    group.add_argument('-autohigh', type=float, default=None,
                       help="The high resolution limit in Angstroms for final high resolution refinement in MR_AUTO "
                            "mode")
    group.add_argument('-hires', type=float, default=None,
                       help="The high resolution limit of data used to find/refine this solution")
    group.add_argument('-hklin', type=str,
                       help="Path the input hkl file")
    group.add_argument('-hklout', type=str,
                       help="Path the output hkl file")
    group.add_argument('-f', type=str,
                       help="The column label for F")
    group.add_argument('-i', type=str,
                       help="The column label for I")
    group.add_argument('-logfile', type=str,
                       help="Path to the ouput log file")
    group.add_argument('-nmol', type=int,
                       help="The predicted number of molecules to build")
    group.add_argument('-pdbin', type=str,
                       help="Path to the input pdb file")
    group.add_argument('-pdbout', type=str,
                       help="Path to the output pdb file")
    group.add_argument('-sgalternative', choices=SGAlternatives.__members__.keys(),
                       help="Try alternative space groups")
    group.add_argument('-sigf', type=str,
                       help="The column label for SIGF")
    group.add_argument('-sigi', type=str,
                       help="The column label for SIGI")
    group.add_argument('-solvent', type=float,
                       help="The estimated solvent content of the crystal")
    group.add_argument('-timeout', type=int,
                       help="The time in mins before phaser will kill a job")
    group.add_argument('-work_dir', type=str,
                       help="Path to the working directory")
    args = parser.parse_args()

    phaser = Phaser(args.hklin, args.hklout, args.f, args.i, args.logfile, args.nmol, args.pdbin, args.pdbout,
                    args.sgalternative, args.sigf, args.sigi, args.solvent, args.timeout, args.work_dir,
                    args.hires, args.autohigh)
    phaser.run()

