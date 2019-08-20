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

import simbad.util.mtz_util
import simbad.parsers.anode_parser
from simbad.core.anode_score import AnomScore

logger = logging.getLogger(__name__)


class AnodeSearch(object):
    """An anomalous phased fourier running class

    Attributes
    ----------
    mtz : str
        The path to the input MTZ
    output_dir : str
        The path to the output directory

    Example
    -------
    >>> from simbad.mr.anomalous_util import AnodeSearch
    >>> anomalous_search = AnodeSearch("<mtz>", "<output_dir>")
    >>> anomalous_search.run("<model>")
    """

    def __init__(self, mtz, work_dir):
        self._mtz = None
        self._space_group = None
        self._resolution = None
        self._cell_parameters = None
        self._output_dir = None

        self.name = None
        self.mtz_labels = None
        self.mtz = mtz
        self.work_dir = work_dir

    @property
    def mtz(self):
        """The input MTZ file"""
        return self._mtz

    @mtz.setter
    def mtz(self, mtz):
        """Define the input MTZ file"""
        self._mtz = mtz

    @property
    def work_dir(self):
        """The path to the working directory"""
        return self._work_dir

    @work_dir.setter
    def work_dir(self, work_dir):
        """Define the working directory"""
        self._work_dir = work_dir

    def run(self, input_model, cleanup=True):
        """Function to run SFALL/CAD/FFT to create phased anomalous fourier map"""
        if not os.path.isdir(self.work_dir):
            os.mkdir(self.work_dir)

        self._space_group, self._resolution, cell_parameters = simbad.util.mtz_util.crystal_data(self.mtz)
        self._cell_parameters = " ".join(map(str, cell_parameters))
        self.mtz_labels = simbad.util.mtz_util.GetLabels(self.mtz)

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
        anode_parser = simbad.parsers.anode_parser.AnodeParser(lsa_file)

        score = AnomScore(dano_peak_height=anode_parser.peak_height, nearest_atom=anode_parser.nearest_atom)
        return score

    def mtz2sca(self):
        sca_out = os.path.join(self.work_dir, "{0}.sca".format(self.name))
        cmd = ["mtz2sca" + EXE_EXT, self.mtz, sca_out]

        if self.mtz_labels.iplus:
            cmd += ["-p", self.mtz_labels.iplus, "-P", self.mtz_labels.sigiplus, "-m", self.mtz_labels.iminus, "-M", self.mtz_labels.sigiminus]
        elif self.mtz_labels.fplus:
            cmd += ["-p", self.mtz_labels.fplus, "-P", self.mtz_labels.sigfplus, "-m", self.mtz_labels.fminus, "-M", self.mtz_labels.sigfminus]
        cexec(cmd)

    def shelxc(self):
        cmd = ["shelxc" + EXE_EXT, self.name]
        stdin = """CELL {0}
SPAG {1}
SAD {2}"""
        sca_out = os.path.join(self.work_dir, "{0}.sca".format(self.name))
        stdin = stdin.format(self._cell_parameters, self._space_group, os.path.relpath(sca_out))
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
