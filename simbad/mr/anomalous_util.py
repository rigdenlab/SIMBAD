#!/usr/bin/env ccp4-python

"""Class to run an anomalous phased fourier on MR results"""

__author__ = "Adam Simpkin"
__date__ = "16 April 2018"
__version__ = "0.2"

import logging
import os
import shutil

from pyjob.script import EXE_EXT
from pyjob import cexec

from simbad.core.anode_score import AnomScore
from simbad.parsers import anode_parser, mtz_parser

logger = logging.getLogger(__name__)


class AnodeSearch(object):
    """An anomalous phased fourier running class

    Attributes
    ----------
    mtz : str
        The path to the input MTZ
    work_dir : str
        The path to the work directory

    Example
    -------
    >>> from simbad.mr.anomalous_util import AnodeSearch
    >>> anomalous_search = AnodeSearch("<mtz>", "<work_dir>")
    >>> anomalous_search.run("<model>")
    """

    def __init__(self, mtz, work_dir):
        self._mtz = None
        self._mtz_obj = None
        self._output_dir = None

        self.name = None
        self.mtz = mtz
        self.work_dir = work_dir

    @property
    def mtz(self):
        """The input MTZ file"""
        return self._mtz

    @mtz.setter
    def mtz(self, mtz):
        """Define the input MTZ file"""
        self._mtz = os.path.abspath(mtz)
        self._mtz_obj = mtz_parser.MtzParser(mtz)
        self._mtz_obj.parse()

    @property
    def mtz_obj(self):
        """Column object containing info on input mtz"""
        return self._mtz_obj

    @property
    def work_dir(self):
        """The path to the working directory"""
        return self._work_dir

    @work_dir.setter
    def work_dir(self, work_dir):
        """Define the working directory"""
        self._work_dir = work_dir

    def run(self, input_model, cleanup=True):
        """Function to run ANODE to create phased anomalous fourier map"""
        if not os.path.isdir(self.work_dir):
            os.mkdir(self.work_dir)

        self.name = os.path.basename(input_model).split(".")[0]
        cwd = os.getcwd()
        os.chdir(self.work_dir)
        self.mtz2sca()
        self.shelxc()
        self.anode(input_model)
        if cleanup:
            self.cleanup()
        os.chdir(cwd)

    def search_results(self):
        """Function to extract search results"""
        lsa_file = os.path.join(self.work_dir, "{0}.lsa".format(self.name))
        ap = anode_parser.AnodeParser(lsa_file)

        score = AnomScore(dano_peak_height=ap.peak_height, nearest_atom=ap.nearest_atom)
        return score

    def mtz2sca(self):
        sca_out = os.path.join(self.work_dir, "{0}.sca".format(self.name))
        cmd = ["mtz2sca" + EXE_EXT, self.mtz, sca_out]

        if self.mtz_obj.i_plus:
            cmd += ["-p", self.mtz_obj.i_plus, "-P", self.mtz_obj.sigi_plus,
                    "-m", self.mtz_obj.i_minus, "-M", self.mtz_obj.sigi_minus]
        elif self.mtz_obj.f_plus:
            cmd += ["-p", self.mtz_obj.f_plus, "-P", self.mtz_obj.sigf_plus,
                    "-m", self.mtz_obj.f_minus, "-M", self.mtz_obj.sigf_minus]
        cexec(cmd)

    def shelxc(self):
        cmd = ["shelxc" + EXE_EXT, self.name]
        stdin = """CELL {0} {1} {2} {3} {4} {5}
SPAG {6}
SAD {7}"""
        sca_out = os.path.join(self.work_dir, "{0}.sca".format(self.name))

        stdin = stdin.format(self.mtz_obj.cell.a,
                             self.mtz_obj.cell.b,
                             self.mtz_obj.cell.c,
                             self.mtz_obj.cell.alpha,
                             self.mtz_obj.cell.beta,
                             self.mtz_obj.cell.gamma,
                             "".join(self.mtz_obj.spacegroup_symbol.split()),
                             os.path.relpath(sca_out))
        cexec(cmd, stdin=stdin)

    def anode(self, input_model):
        shutil.copyfile(input_model, os.path.join(self.work_dir, self.name + ".pdb"))
        cmd = ["anode" + EXE_EXT, self.name]
        cexec(cmd)

    def cleanup(self):
        for i in ["{0}_fa.hkl", "{0}_fa.ins", "{0}_fa.res", "{0}.hkl", "{0}.sca"]:
            f = os.path.join(self.work_dir, i.format(self.name))
            if os.path.isfile(f):
                os.remove(f)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Wrapper for running ANODE", prefix_chars="-")
    group = parser.add_argument_group()

    group.add_argument("-xyzin", type=str, help="Path to the input xyz file")
    group.add_argument("-hklin", type=str, help="Path to the input hkl file")
    group.add_argument("-work_dir", type=str, help="Path to the working directory", default=os.getcwd())
    group.add_argument("--cleanup", default=False, action="store_true", help="Delete all none essential files")

    args = parser.parse_args()

    anomalous_search = AnodeSearch(os.path.abspath(args.hklin), os.path.abspath(args.work_dir))
    anomalous_search.run(os.path.abspath(args.xyzin), cleanup=args.cleanup)
