"""Class to run an anomalous phased fourier on MR results"""

__author__ = "Adam Simpkin"
__date__ = "16 April 2018"
__version__ = "0.2"

import logging
import os
import shutil

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

    def __init__(self, mtz, output_dir, mr_program):
        self._mtz = None
        self._mr_program = None
        self._space_group = None
        self._resolution = None
        self._cell_parameters = None
        self._output_dir = None

        self.name = None
        self.mtz_labels = None
        self.mr_program = mr_program
        self.mtz = mtz
        self.output_dir = output_dir
        self.work_dir = None

    @property
    def mr_program(self):
        """The molecular replacement program to use"""
        return self._mr_program

    @mr_program.setter
    def mr_program(self, mr_program):
        """Define the molecular replacement program to use"""
        self._mr_program = mr_program

    @property
    def mtz(self):
        """The input MTZ file"""
        return self._mtz

    @mtz.setter
    def mtz(self, mtz):
        """Define the input MTZ file"""
        self._mtz = mtz

    @property
    def output_dir(self):
        """The path to the working directory"""
        return self._output_dir

    @output_dir.setter
    def output_dir(self, output_dir):
        """Define the working directory"""
        self._output_dir = output_dir

    def run(self, model):
        """Function to run SFALL/CAD/FFT to create phased anomalous fourier map"""
        self.work_dir = os.path.join(self.output_dir, model.pdb_code, "anomalous")
        os.mkdir(self.work_dir)

        self._space_group, self._resolution, cell_parameters = simbad.util.mtz_util.crystal_data(self.mtz)
        self._cell_parameters = " ".join(map(str, cell_parameters))
        self.mtz_labels = simbad.util.mtz_util.GetLabels(self.mtz)

        input_model = os.path.join(self.output_dir, model.pdb_code, "mr",
                                   self.mr_program, "{0}_mr_output.pdb".format(model.pdb_code))
        self.name = model.pdb_code

        cwd = os.getcwd()
        os.chdir(self.work_dir)
        self.mtz2sca()
        self.shelxc()
        self.anode(input_model)
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
        cmd = ["mtz2sca", self.mtz, sca_out]

        if self.mtz_labels.iplus:
            cmd += ['-p', self.mtz_labels.iplus,
                    '-P', self.mtz_labels.sigiplus,
                    '-m', self.mtz_labels.iminus,
                    '-M', self.mtz_labels.sigiminus]
        elif self.mtz_labels.fplus:
            cmd += ['-p', self.mtz_labels.fplus,
                    '-P', self.mtz_labels.sigfplus,
                    '-m', self.mtz_labels.fminus,
                    '-M', self.mtz_labels.sigfminus]
        cexec(cmd)

    def shelxc(self):
        cmd = ["shelxc", self.name]
        stdin = """CELL {0}
SPAG {1}
SAD {2}"""
        sca_out = os.path.join(self.work_dir, "{0}.sca".format(self.name))
        stdin = stdin.format(self._cell_parameters, self._space_group, os.path.relpath(sca_out))
        cexec(cmd, stdin=stdin)

    def anode(self, input_model):
        shutil.copyfile(input_model, os.path.join(self.work_dir, self.name + ".pdb"))
        cmd = ["anode", self.name]
        cexec(cmd)

    def cleanup(self):
        for i in ["{0}_fa.hkl", "{0}_fa.ins", "{0}_fa.res", "{0}.hkl", "{0}.pha", "{0}.sca"]:
            os.remove(os.path.join(self.work_dir, i.format(self.name)))

