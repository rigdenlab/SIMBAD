#!/usr/bin/env ccp4-python
"""Module to run phaser rotation search on a model"""

__author__ = "Adam Simpkin"
__date__ = "12 April 2018"
__version__ = "1.0"

import os

from phaser import InputMR_DAT, runMR_DAT, InputMR_FRF, runMR_FRF


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
    pdbin : str
        Path to the input pdb file
    pdbout : str
        Path to the output pdb file
    sigf : str
        The column label for SIGF
    sigi : str
        The column label for SIGI
    solvent : int float
        The estimated solvent content of the crystal
    work_dir : str
        Path to the working directory were you want PHASER to run
    hires : str
        The high resolution limit of data used to find/refine this solution


    Examples
    --------
    >>> from simbad.rotsearch.phaser_rotation_search import Phaser
    >>> phaser = Phaser('<hklin>', '<f>', '<i>', '<logfile>', '<nmol>', '<pdbin>', '<pdbout>', '<sgalternative>',
    >>>                 '<sigf>', '<sigi>', '<solvent>', '<timeout>', '<workdir>', '<autohigh>', '<hires>')
    >>> phaser.run()

    Files relating to the PHASER run will be contained within the work_dir however the location of the output hkl, pdb
    and logfile can be specified.
    """

    def __init__(self, hklin, f, i, logfile, nmol, pdbin, sigf, sigi, solvent, timeout, work_dir, hires):
        self._f = None
        self._hires = None
        self._hklin = None
        self._i = None
        self._logfile = None
        self._nmol = None
        self._pdbin = None
        self._sigf = None
        self._sigi = None
        self._solvent = None
        self._timeout = None
        self._work_dir = None

        self.f = f
        self.hires = hires
        self.hklin = hklin
        self.i = i
        self.logfile = logfile
        self.nmol = nmol
        self.pdbin = pdbin
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
    def i(self):
        """The I label from the input hkl"""
        return self._i

    @i.setter
    def i(self, i):
        """Define the I label from the input hkl"""
        self._i = i

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
        """Function to run rotation search using PHASER"""

        current_work_dir = os.getcwd()
        if os.path.exists(self.work_dir):
            os.chdir(self.work_dir)
        else:
            os.makedirs(self.work_dir)
            os.chdir(self.work_dir)

        i = InputMR_DAT()
        i.setHKLI(self.hklin)

        if self.hires:
            i.setHIRES(self.hires)
        if self.i != "None" and self.sigi != "None":
            i.setLABI_I_SIGI(self.i, self.sigi)
        elif self.f != "None" and self.sigf != "None":
            i.setLABI_F_SIGF(self.f, self.sigf)
        else:
            msg = "No flags for intensities or amplitudes have been provided"
            raise RuntimeError(msg)
        i.setMUTE(True)
        run_mr_data = runMR_DAT(i)

        if run_mr_data.Success():
            i = InputMR_FRF()
            i.setJOBS(1)
            i.setREFL_DATA(run_mr_data.getREFL_DATA())
            i.setSPAC_HALL(run_mr_data.getSpaceGroupHall())
            i.setCELL6(run_mr_data.getUnitCell())
            i.setROOT("phaser_mr_output")
            i.addENSE_PDB_RMS("PDB", self.pdbin, 0.6)
            i.setCOMP_BY("SOLVENT")
            i.setCOMP_PERC(self.solvent)
            i.addSEAR_ENSE_NUM('PDB', self.nmol)
            if self.timeout != 0:
                i.setKILL_TIME(self.timeout)
            run_mr_rot = runMR_FRF(i)

            with open(self.logfile, 'w') as f:
                f.write(run_mr_rot.summary())

        os.chdir(current_work_dir)

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Runs rotation search using PHASER', prefix_chars="-")

    group = parser.add_argument_group()
    group.add_argument('-f', type=str,
                       help="The column label for F")
    group.add_argument('-hires', type=float, default=None,
                       help="The high resolution limit of data used to find/refine this solution")
    group.add_argument('-hklin', type=str,
                       help="Path the input hkl file")
    group.add_argument('-i', type=str,
                       help="The column label for I")
    group.add_argument('-logfile', type=str,
                       help="Path to the ouput log file")
    group.add_argument('-nmol', type=int,
                       help="The predicted number of molecules to build")
    group.add_argument('-pdbin', type=str,
                       help="Path to the input pdb file")
    group.add_argument('-sigf', type=str,
                       help="The column label for SIGF")
    group.add_argument('-sigi', type=str,
                       help="The column label for SIGI")
    group.add_argument('-solvent', type=float,
                       help="The estimated solvent content of the crystal")
    group.add_argument('-timeout', type=int, default=0,
                       help="The time in mins before phaser will kill a job")
    group.add_argument('-work_dir', type=str,
                       help="Path to the working directory")
    args = parser.parse_args()

    phaser = Phaser(args.hklin, args.f, args.i, args.logfile, args.nmol, args.pdbin, args.sigf, args.sigi, args.solvent,
                    args.timeout, args.work_dir, args.hires)
    phaser.run()