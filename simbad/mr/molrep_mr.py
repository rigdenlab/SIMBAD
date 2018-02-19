#!/usr/bin/env ccp4-python
"""Module to run molrep on a model"""

__author__ = "Adam Simpkin"
__date__ = "02 May 2017"
__version__ = "1.0"

import os
import operator
import shutil

from pyjob import cexec


def check_contrast(logfile):
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
            if "Contrast" in line:
                fields = line.split()
                if len(fields) > 3:
                    return 0.0
                else:
                    return float(fields[-1])


class Molrep(object):
    """Class to run Molrep

    Attributes
    ----------
    hklin : str
        Path to input hkl file
    logfile : str
        Path to output log file
    pdbin : str
        Path to the input pdb file
    pdbout : str
        Path to the output pdb file
    sgalternative : str
        Specify whether to try alternative space groups (all | enant)
    space_group : str
        The space group of the input hkl file
    work_dir : str
        Path to the working directory were you want MOLREP to run

    Example
    -------
    >>> from simbad.mr.molrep_mr import Molrep
    >>> molrep = Molrep('<hklin>', '<logfile>', '<nmol>', '<pdbin>', '<pdbout>', '<sgalternative>', '<space_group>',
    >>>                 '<work_dir>')
    >>> molrep.run()

    Files relating to the MOLREP run will be contained within the work_dir however the location of the output pdb and
    logfile can be specified.
    """

    def __init__(self, hklin, logfile, nmol, pdbin, pdbout, sgalternative, space_group, work_dir):

        self._hklin = None
        self._logfile = None
        self._nmol = None
        self._pdbin = None
        self._pdbout = None
        self._work_dir = None
        self._sgalternative = None
        self._space_group = None

        self.hklin = hklin
        self.logfile = logfile
        self.nmol = nmol
        self.pdbin = pdbin
        self.pdbout = pdbout
        self.sgalternative = sgalternative
        self.space_group = space_group
        self.work_dir = work_dir

        self.all_sg_codes = {"P2": "3",
                             "P21": "4",
                             "C2": "5",
                             "P222": "16",
                             "P2221": "17",
                             "P21212": "18",
                             "P212121": "19",
                             "C2221": "20",
                             "C222": "21",
                             "F222": "22",
                             "I222": "23",
                             "I212121": "24",
                             "P4": "75",
                             "P41": "76",
                             "P42": "77",
                             "P43": "78",
                             "I4": "79",
                             "I41": "80",
                             "P422": "89",
                             "P4212": "90",
                             "P4122": "91",
                             "P41212": "92",
                             "P4222": "93",
                             "P42212": "94",
                             "P4322": "95",
                             "P43212": "96",
                             "I422": "97",
                             "I4122": "98",
                             "P3": "143",
                             "P31": "144",
                             "P32": "145",
                             "R3": "146",
                             "P312": "149",
                             "P321": "150",
                             "P3112": "151",
                             "P3121": "152",
                             "P3212": "153",
                             "P3221": "154",
                             "R32": "155",
                             "P6": "168",
                             "P61": "169",
                             "P65": "170",
                             "P62": "171",
                             "P64": "172",
                             "P63": "173",
                             "P622": "177",
                             "P6122": "178",
                             "P6522": "179",
                             "P6222": "180",
                             "P6422": "181",
                             "P6322": "182",
                             "P23": "195",
                             "F23": "196",
                             "I23": "197",
                             "P213": "198",
                             "I213": "199",
                             "P432": "207",
                             "P4232": "208",
                             "F432": "209",
                             "F4132": "210",
                             "I432": "211",
                             "P4332": "212",
                             "P4132": "213",
                             "I4132": "214"}

        self.enant_sg_codes = {"P31": "144",
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

        self.all_alt_sg = {"3": ["4", "5"],
                           "4": ["3", "5"],
                           "5": ["3", "4"],
                           "16": ["17", "18", "19", "20", "21", "22", "23", "24"],
                           "17": ["16", "18", "19", "20", "21", "22", "23", "24"],
                           "18": ["16", "17", "19", "20", "21", "22", "23", "24"],
                           "19": ["16", "17", "18", "20", "21", "22", "23", "24"],
                           "20": ["16", "17", "18", "19", "21", "22", "23", "24"],
                           "21": ["16", "17", "18", "19", "20", "22", "23", "24"],
                           "22": ["16", "17", "18", "19", "20", "21", "23", "24"],
                           "23": ["16", "17", "18", "19", "20", "21", "22", "24"],
                           "24": ["16", "17", "18", "19", "20", "21", "22", "23"],
                           "75": ["76", "77", "78", "79", "80"],
                           "76": ["75", "77", "78", "79", "80"],
                           "77": ["75", "76", "78", "79", "80"],
                           "78": ["75", "76", "77", "79", "80"],
                           "79": ["75", "76", "77", "78", "80"],
                           "80": ["75", "76", "77", "78", "79"],
                           "89": ["90", "91", "92", "93", "94", "95", "96", "97", "98"],
                           "90": ["89", "91", "92", "93", "94", "95", "96", "97", "98"],
                           "91": ["89", "90", "92", "93", "94", "95", "96", "97", "98"],
                           "92": ["89", "90", "91", "93", "94", "95", "96", "97", "98"],
                           "93": ["89", "90", "91", "92", "94", "95", "96", "97", "98"],
                           "94": ["89", "90", "91", "92", "93", "95", "96", "97", "98"],
                           "95": ["89", "90", "91", "92", "93", "94", "96", "97", "98"],
                           "96": ["89", "90", "91", "92", "93", "94", "95", "97", "98"],
                           "97": ["89", "90", "91", "92", "93", "94", "95", "96", "98"],
                           "98": ["89", "90", "91", "92", "93", "94", "95", "96", "97"],
                           "143": ["144", "145", "146"],
                           "144": ["143", "145", "146"],
                           "145": ["143", "144", "146"],
                           "146": ["143", "144", "145"],
                           "149": ["150", "151", "152", "153", "154", "155"],
                           "150": ["149", "151", "152", "153", "154", "155"],
                           "151": ["149", "150", "152", "153", "154", "155"],
                           "152": ["149", "150", "151", "153", "154", "155"],
                           "153": ["149", "150", "151", "152", "154", "155"],
                           "154": ["149", "150", "151", "152", "153", "155"],
                           "155": ["149", "150", "151", "152", "153", "154"],
                           "168": ["169", "170", "171", "172", "173"],
                           "169": ["168", "170", "171", "172", "173"],
                           "170": ["168", "169", "171", "172", "173"],
                           "171": ["168", "169", "170", "172", "173"],
                           "172": ["168", "169", "170", "171", "173"],
                           "173": ["168", "169", "170", "171", "172"],
                           "177": ["178", "179", "180", "181", "182"],
                           "178": ["177", "179", "180", "181", "182"],
                           "179": ["177", "178", "180", "181", "182"],
                           "180": ["177", "178", "179", "181", "182"],
                           "181": ["177", "178", "179", "180", "182"],
                           "182": ["177", "178", "179", "180", "181"],
                           "195": ["196", "197", "198", "199"],
                           "196": ["195", "197", "198", "199"],
                           "197": ["195", "196", "198", "199"],
                           "198": ["195", "196", "197", "199"],
                           "199": ["195", "196", "197", "198"],
                           "207": ["208", "209", "210", "211", "212", "213", "214"],
                           "208": ["207", "209", "210", "211", "212", "213", "214"],
                           "209": ["207", "208", "210", "211", "212", "213", "214"],
                           "210": ["207", "208", "209", "211", "212", "213", "214"],
                           "211": ["207", "208", "209", "210", "212", "213", "214"],
                           "212": ["207", "208", "209", "210", "211", "213", "214"],
                           "213": ["207", "208", "209", "210", "211", "212", "214"],
                           "214": ["207", "208", "209", "210", "211", "212", "213"]}

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
        if sgalternative:
            self._sgalternative = sgalternative.lower()
        else:
            self._sgalternative = sgalternative
        
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
        template_key = """
        FILE_F {0}
        FILE_M {1}
        NMON {2}
        {3}
        END"""
        key = template_key.format(os.path.relpath(hklin), os.path.relpath(pdbin), self.nmol, "")
        self.molrep(key, logfile)

        # Move output pdb to specified name
        if os.path.isfile(os.path.join(self.work_dir, "molrep.pdb")):
            shutil.move(os.path.join(self.work_dir, "molrep.pdb"),
                        os.path.join(self.work_dir, 'molrep_out_{0}.pdb'.format(self.space_group)))

        if self.sgalternative == "enant" and self.space_group in self.enant_sg_codes:
            hklin_sg_code = self.enant_sg_codes[self.space_group]
            enant_sg_code = self.enant_sg[hklin_sg_code]
            contrast = check_contrast(os.path.join(self.work_dir, 'molrep_out_{0}.log'.format(self.space_group)))
            contrasts = {self.space_group: contrast}
            key = template_key.format(os.path.relpath(hklin), os.path.relpath(pdbin), self.nmol,
                                      "NOSG {0}".format(enant_sg_code))
            logfile = os.path.join(self.work_dir, 'molrep_out_{0}.log'.format(enant_sg_code))

            self.molrep(key, logfile)
            contrasts[enant_sg_code] = (check_contrast(logfile))

            if os.path.isfile(os.path.join(self.work_dir, "molrep.pdb")):
                shutil.move(os.path.join(self.work_dir, "molrep.pdb"),
                            os.path.join(self.work_dir, 'molrep_out_{0}.pdb'.format(enant_sg_code)))
            self.evaluate_results(contrasts)

        elif self.sgalternative == "all" and self.space_group in self.all_sg_codes:
            hklin_sg_code = self.all_sg_codes[self.space_group]
            all_alt_sg_codes = self.all_alt_sg[hklin_sg_code]
            contrast = check_contrast(os.path.join(self.work_dir, 'molrep_out_{0}.log'.format(self.space_group)))
            contrasts = {self.space_group: contrast}
            for alt_sg_code in all_alt_sg_codes:
                key = template_key.format(os.path.relpath(hklin), os.path.relpath(pdbin), self.nmol,
                                          "NOSG {0}".format(alt_sg_code))
                logfile = os.path.join(self.work_dir, 'molrep_out_{0}.log'.format(alt_sg_code))
                self.molrep(key, logfile)
                contrasts[alt_sg_code] = (check_contrast(logfile))
                if os.path.isfile(os.path.join(self.work_dir, "molrep.pdb")):
                    shutil.move(os.path.join(self.work_dir, "molrep.pdb"),
                                os.path.join(self.work_dir, 'molrep_out_{0}.pdb'.format(alt_sg_code)))
            self.evaluate_results(contrasts)

            
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

    def evaluate_results(self, results):
        """Function to evaluate molrep results and move the result with the best contrast score to the output pdb

        Parameters
        ----------
        results : dict
            Dictionary containing space group code with the corresponding contrast score

        Returns
        -------
        file
            The output pdb for the best result
        file
            The output log for the best result
        """

        top_sg_code = max(results.iteritems(), key=operator.itemgetter(1))[0]
        if os.path.isfile(os.path.join(self.work_dir, 'molrep_out_{0}.pdb'.format(top_sg_code))):
            shutil.move(os.path.join(self.work_dir, 'molrep_out_{0}.pdb'.format(top_sg_code)), self.pdbout)
        if os.path.isfile(os.path.join(self.work_dir, 'molrep_out_{0}.log'.format(top_sg_code))):
            shutil.move(os.path.join(self.work_dir, 'molrep_out_{0}.log'.format(top_sg_code)), self.logfile)
    
    @staticmethod
    def molrep(key, logfile):
        """Function to run molecular replacement using MOLREP

        Parameters
        ----------
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

        cmd = ["molrep"]
        stdout = cexec(cmd, stdin=key)
        with open(logfile, 'w') as f_out:
            f_out.write(stdout)


if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser(description='Runs MR using MOLREP', prefix_chars="-")

    group = parser.add_argument_group()
    group.add_argument('-hklin', type=str,
                       help="Path the input hkl file")
    group.add_argument('-logfile', type=str,
                       help="Path to the ouput log file")
    group.add_argument('-nmol', type=int,
                       help="The predicted number of molecules to build")
    group.add_argument('-pdbin', type=str,
                       help="Path to the input pdb file")
    group.add_argument('-pdbout', type=str,
                       help="Path to the output pdb file")
    group.add_argument('-sgalternative', default=None,
                       help="Try alternative space groups <all/enant>")
    group.add_argument('-space_group', type=str,
                       help="The input space group")
    group.add_argument('-work_dir', type=str,
                       help="Path to the working directory")
    args = parser.parse_args()
    
    molrep = Molrep(args.hklin, args.logfile, args.nmol, args.pdbin, args.pdbout, args.sgalternative, args.space_group,
                    args.work_dir)
    molrep.run()
