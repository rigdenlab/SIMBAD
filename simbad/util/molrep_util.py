"""Module to run molrep on a model"""

import os
import simbad_util
import shutil

__author__ = "Adam Simpkin"
__date__ = "02 May 2017"
__version__ = "0.1"


class Molrep(object):
    """Class to run Molrep

    Attributes
    ----------
    enant : bool
        Specify whether to try enantimorphic space groups
    hklin : str
        Path to input hkl file
    logfile : str
        Path to output log file
    pdbin : str
        Path to the input pdb file
    pdbout : str
        Path to the output pdb file
    space_group : str
        The space group of the input hkl file
    work_dir : str
        Path to the working directory were you want MOLREP to run

    Example
    -------
    >>> from simbad.util.molrep_util import Molrep
    >>> molrep = Molrep('<enant>', '<hklin>', '<logfile>', '<pdbin>', '<pdbout>', '<space_group>', '<work_dir>')
    >>> molrep.run()

    Files relating to the MOLREP run will be contained within the work_dir however the location of the output pdb and
    logfile can be specified.
    """

    def __init__(self, enant, hklin, logfile, pdbin, pdbout, space_group, work_dir):

        self._enant = None
        self._hklin = None
        self._logfile = None
        self._pdbin = None
        self._pdbout = None
        self._work_dir = None
        self._space_group = None

        self.enant = enant
        self.hklin = hklin
        self.logfile = logfile
        self.pdbin = pdbin
        self.pdbout = pdbout
        self.space_group = space_group
        self.work_dir = work_dir

        self.sg_codes = {"P31": "144",
                         "P32": "145",
                         "P3112": "151",
                         "P3212": "153",
                         "P3121": "152",
                         "P3221": "154",
                         "P41": "76",
                         "P43": "78",
                         "P4122": "91",
                         "P4322": "95",
                         "P41212": "92",
                         "P43212": "96",
                         "P61": "169",
                         "P65": "170",
                         "P62": "171",
                         "P64": "172",
                         "P6122": "178",
                         "P6522": "179",
                         "P6222": "180",
                         "P6422": "181",
                         "P4332": "212",
                         "P4132": "213"}

        self.enant_sg = {"144": "145",
                         "145": "144",
                         "151": "153",
                         "153": "151",
                         "152": "154",
                         "154": "152",
                         "76": "78",
                         "78": "76",
                         "91": "95",
                         "95": "91",
                         "92": "96",
                         "96": "92",
                         "169": "170",
                         "170": "169",
                         "171": "172",
                         "172": "171",
                         "178": "179",
                         "179": "178",
                         "180": "181",
                         "181": "180",
                         "212": "213",
                         "213": "212"}

    @property
    def enant(self):
        """Whether to check for enantimophic space groups"""
        return self._enant

    @enant.setter
    def enant(self, enant):
        """Define whether to check for enantiomorphic space groups"""
        self._enant = enant

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

    def check_contrast(self, logfile):
        """Check the logfile of the job for the contrast value

        Parameters
        ----------
        logfile : str
            Path to the logfile

        Returns
        -------
        float
            Contrast score in log file
        """

        with open(logfile, 'r') as f:
            for line in f:
                if "Contrast =" in line:
                    return float(line.split()[-1])
                else:
                    return 0.0

    def run(self):
        """Function to run molecular replacement using MOLREP

        Returns
        -------
        file
            The output pdb from MOLREP
        """

        # Make a note of the current working directory
        current_work_dir = os.getcwd()

        # Change to the MOLREP working directory
        os.chdir(self.work_dir)

        # Copy hklin and pdbin to working dire for efficient running of MOLREP
        hklin = os.path.join(self.work_dir, os.path.basename(self.hklin))
        shutil.copyfile(self.hklin, hklin)
        pdbin = os.path.join(self.work_dir, os.path.basename(self.pdbin))
        shutil.copyfile(self.pdbin, pdbin)
        logfile = os.path.join(self.work_dir, 'molrep_out_{0}.log'.format(self.space_group))
        key = ''
        self.molrep(hklin, pdbin, key, logfile)

        # Move output pdb to specified name
        if os.path.isfile(os.path.join(self.work_dir, "molrep.pdb")):
            shutil.move(os.path.join(self.work_dir, "molrep.pdb"),
                        os.path.join(self.work_dir, 'molrep_out_{0}.pdb'.format(self.space_group)))

        if self.enant:
            hklin_sg_code = self.sg_codes[self.space_group]
            enant_sg_code = self.enant_sg[hklin_sg_code]
            key = 'NOSG {0}'.format(enant_sg_code)
            logfile = os.path.join(self.work_dir, 'molrep_out_{0}.log'.format(enant_sg_code))

            self.molrep(hklin, pdbin, key, logfile)

            # Move output pdb to specified name
            if os.path.isfile(os.path.join(self.work_dir, "molrep.pdb")):
                shutil.move(os.path.join(self.work_dir, "molrep.pdb"),
                            os.path.join(self.work_dir, 'molrep_out_{0}.pdb'.format(enant_sg_code)))

            contrast_1 = self.check_contrast(os.path.join(self.work_dir, 'molrep_out_{0}.log'.format(self.space_group)))
            contrast_2 = self.check_contrast(logfile)

            if contrast_1 > contrast_2:
                # Move output pdb to specified name
                if os.path.isfile(os.path.join(self.work_dir, 'molrep_out_{0}.pdb'.format(self.space_group))):
                    shutil.move(os.path.join(self.work_dir, 'molrep_out_{0}.pdb'.format(self.space_group)), self.pdbout)
                # Move log file to specified name
                if os.path.isfile(os.path.join(self.work_dir, 'molrep_out_{0}.log'.format(self.space_group))):
                    shutil.move(os.path.join(self.work_dir, 'molrep_out_{0}.log'.format(self.space_group)),
                                self.logfile)
            elif contrast_2 > contrast_1:
                # Move output pdb to specified name
                if os.path.isfile(os.path.join(self.work_dir, 'molrep_out_{0}.pdb'.format(enant_sg_code))):
                    shutil.move(os.path.join(self.work_dir, 'molrep_out_{0}.pdb'.format(enant_sg_code)), self.pdbout)
                # Move log file to specified name
                if os.path.isfile(os.path.join(self.work_dir, 'molrep_out_{0}.log'.format(enant_sg_code))):
                    shutil.move(os.path.join(self.work_dir, 'molrep_out_{0}.log'.format(enant_sg_code)), self.logfile)

        else:
            # Move output pdb to specified name
            if os.path.isfile(os.path.join(self.work_dir, 'molrep_out_{0}.pdb'.format(self.space_group))):
                shutil.move(os.path.join(self.work_dir, 'molrep_out_{0}.pdb'.format(self.space_group)), self.pdbout)
            # Move log file to specified name
            if os.path.isfile(os.path.join(self.work_dir, 'molrep_out_{0}.log'.format(self.space_group))):
                shutil.move(os.path.join(self.work_dir, 'molrep_out_{0}.log'.format(self.space_group)), self.logfile)

        # Return to original working directory
        os.chdir(current_work_dir)

        # Delete any files copied across
        if os.path.isfile(os.path.join(self.work_dir, os.path.basename(self.hklin))):
            os.remove(os.path.join(self.work_dir, os.path.basename(self.hklin)))
        if os.path.isfile(os.path.join(self.work_dir, os.path.basename(self.pdbin))):
            os.remove(os.path.join(self.work_dir, os.path.basename(self.pdbin)))
        return

    @staticmethod
    def molrep(hklin, pdbin, key, logfile):
        """Function to run molecular replacement using MOLREP

        Parameters
        ----------
        hklin : str
            Path to input hkl file
        pdbin : str
            Path to input pdb file
        key : str
            MOLREP keywords
        logfile :
            Path to output log

        Returns
        -------
        file
            The output pdb from MOLREP
        file
            The output logfile
        """

        cmd = ["molrep",
               "-f", hklin,
               "-m", pdbin]
        simbad_util.run_job(cmd, logfile=logfile, stdin=key)

