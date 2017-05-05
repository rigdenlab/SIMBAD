"""Module to run refmac on a model"""

__author__ = "Adam Simpkin"
__date__ = "02 May 2017"
__version__ = "0.1"

import os

import simbad.util.simbad_util


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

        key = "ncyc {0}".format(ncyc)

        self.refmac(self.hklin, self.hklout, self.pdbin, self.pdbout, self.logfile, key)
        return

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

        cmd = ['refmac5',
               'hklin', hklin,
               'hklout', hklout,
               'xyzin', pdbin,
               'xyzout', pdbout]
        simbad.util.simbad_util.run_job(cmd, logfile=logfile, stdin=key)

