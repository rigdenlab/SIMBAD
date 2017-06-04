#!/usr/bin/env ccp4-python
"""Module to run molrep on a model"""

__author__ = "Adam Simpkin"
__date__ = "02 May 2017"
__version__ = "1.0"

import os
import shutil

from pyjob.dispatch import cexec


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
    >>> from simbad.mr.molrep_mr import Molrep
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
    def space_group(self):
        """The input space group"""
        return self._space_group
    
    @space_group.setter
    def space_group(self, space_group):
        """Define the input space group"""
        self._space_group = space_group

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
        if os.path.exists(self.work_dir):
            os.chdir(self.work_dir)
        else:
            os.makedirs(self.work_dir)
            os.chdir(self.work_dir)

        # Copy hklin and pdbin to working dire for efficient running of MOLREP
        hklin = os.path.join(self.work_dir, os.path.basename(self.hklin))
        shutil.copyfile(self.hklin, hklin)
        pdbin = os.path.join(self.work_dir, os.path.basename(self.pdbin))
        shutil.copyfile(self.pdbin, pdbin)
        logfile = os.path.join(self.work_dir, 'molrep_out_{0}.log'.format(self.space_group))
        key = ''
        Molrep.molrep(hklin, pdbin, key, logfile)

        # Move output pdb to specified name
        if os.path.isfile(os.path.join(self.work_dir, "molrep.pdb")):
            shutil.move(os.path.join(self.work_dir, "molrep.pdb"),
                        os.path.join(self.work_dir, 'molrep_out_{0}.pdb'.format(self.space_group)))

        if self.enant and self.space_group in self.sg_codes:
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
    
    @ staticmethod
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
            Path to output log file

        Returns
        -------
        file
            The output pdb from MOLREP
        file
            The output log file
        """
        cmd = ["molrep", "-f", hklin, "-m", pdbin]
        stdout = cexec(cmd, stdin=key)
        with open(logfile, 'w') as f_out:
            f_out.write(stdout)


if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser(description='Runs MR using MOLREP', prefix_chars="-")

    group = parser.add_argument_group()
    group.add_argument('-enant',
                       help="Try enantimorph space groups <True|False>")
    group.add_argument('-hklin', type=str,
                       help="Path the input hkl file")
    group.add_argument('-logfile', type=str,
                       help="Path to the ouput log file")
    group.add_argument('-pdbin', type=str,
                       help="Path to the input pdb file")
    group.add_argument('-pdbout', type=str,
                       help="Path to the output pdb file")
    group.add_argument('-space_group', type=str,
                       help="The input space group")
    group.add_argument('-work_dir', type=str,
                       help="Path to the working directory")
    args = parser.parse_args()
    
    if args.enant.lower() == 'true':
        enant = True
    elif args.enant.lower() == 'false':
        enant = False
    else:
        raise RuntimeError("Incorrect input for '-enant', use 'True' or 'False'")
    
    molrep = Molrep(enant, args.hklin, args.logfile, args.pdbin, args.pdbout, args.space_group, args.work_dir)
    molrep.run()
