"""Class to run an anomalous phased fourier on MR results"""

__author__ = "Adam Simpkin"
__date__ = "17 Mar 2017"
__version__ = "0.1"

import logging
import numpy
import os

from pyjob import cexec
from scipy import stats

import simbad.util.mtz_util

logger = logging.getLogger(__name__)


class _AnomScore(object):
    """An anomalous phased fourier scoring class"""

    __slots__ = ("dano_peak_height", "dano_z_score")

    def __init__(self, dano_peak_height, dano_z_score):
        self.dano_peak_height = dano_peak_height
        self.dano_z_score = dano_z_score

    def __repr__(self):
        return "{0}(dano_peak_height={1} dano_z_score={2} ".format(self.__class__.__name__,
                                                                   self.dano_peak_height,
                                                                   self.dano_z_score)
    def _as_dict(self):
        """Convert the :obj:`_MrScore <simbad.mr._MrScore>`
        object to a dictionary"""
        return {k: getattr(self, k) for k in self.__slots__}


class AnomSearch(object):
    """An anomalous phased fourier running class

    Attributes
    ----------
    mtz : str
        The path to the input MTZ
    output_dir : str
        The path to the output directory
    model : class
        Class object containing the PDB code for the input model

    Example
    -------
    >>> from simbad.mr.anomalous_util import AnomSearch
    >>> anomalous_search = AnomSearch("<mtz>", "<output_dir>")
    >>> anomalous_search.run("<model>")
    """

    def __init__(self, mtz, output_dir, mr_program):
        self._mtz = None
        self._mr_program = None
        self._f = None
        self._sigf = None
        self._dano = None
        self._sigdano = None
        self._free = None
        self._space_group = None
        self._resolution = None
        self._cell_parameters = None
        self._output_dir = None

        self.name = None
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
        # Make output directory
        self.work_dir = os.path.join(self.output_dir, model.pdb_code, "anomalous")
        os.mkdir(self.work_dir)

        self._f, self._sigf, _, _, self._dano, self._sigdano, self._free = simbad.util.mtz_util.get_labels(self.mtz)
        self._space_group, self._resolution, cell_parameters = simbad.util.mtz_util.crystal_data(self.mtz)
        self._cell_parameters = " ".join(map(str, cell_parameters))

        # Create path to the placed mr solution
        input_model = os.path.join(self.output_dir, model.pdb_code, "mr",
                                   self.mr_program, "{0}_mr_output.pdb".format(model.pdb_code))
        self.name = model.pdb_code

        # Run programs
        self.sfall(input_model)
        self.cad()
        self.fft()
        self.peakmax()

    def search_results(self):
        """Function to extract search results"""
        heavy_atom_file = os.path.join(self.work_dir, "peakmax_{0}.ha".format(self.name))
        all_peaks = []

        with open(heavy_atom_file, 'r') as f:
            for line in f:
                if line.startswith('ATOM'):
                    peak = line.split()[5]
                    all_peaks.append(peak)

        z_score = stats.zscore(numpy.array(all_peaks, dtype=numpy.float)).max()

        score = _AnomScore(dano_peak_height=all_peaks[0],
                           dano_z_score=z_score)
        return score

    def sfall(self, model):
        """Function to run SFALL to calculated structure factors for the placed MR model

        Parameters
        ----------
        model : str
            path to placed model from MR
        self.name : str
            unique identifier for the input model set by :obj:`AnomSearch.run`
        self.mtz : str
            mtz file input to :obj:`AnomSearch`
        self.work_dir : str
            working directory set by :obj:`AnomSearch.run`
        self._f : str
            f column label set by :obj: `AnomSearch`
        self._sigf : str
            sigf column label set by :obj: `AnomSearch`
        self._free : str
            free column label set by :obj: `AnomSearch`

        Returns
        -------
        file
            mtz file containing FCalc and PHICalc columns

        """
        cmd = ["sfall", "HKLOUT", os.path.join(self.work_dir, "sfall_{0}.mtz".format(self.name)),
               "XYZIN", model, "HKLIN", self.mtz]
        stdin = os.linesep.join([
            "LABIN  FP={0} SIGFP={1} FREE={2}",
            "labout -",
            "FC=FCalc PHIC=PHICalc",
            "MODE SFCALC -",
            "   XYZIN -",
            "   HKLIN",
            "symmetry '{3}'",
            "badd 0.0",
            "vdwr 2.5",
            "end",
        ])
        stdin = stdin.format(self._f, self._sigf, self._free, self._space_group)
        cexec(cmd, stdin=stdin)

    def cad(self):
        """Function to run CAD to combine the calculated structure factors and the anomalous signal

        Parameters
        ----------
        self.name : str
            unique identifier for the input model set by :obj:`AnomSearch.run`
        self.mtz : str
            mtz file input to :obj: `AnomSearch`
        self.work_dir : str
            working directory set by :obj:`AnomSearch.run`
        self._f : str
            f column label set by :obj: `AnomSearch`
        self._sigf : str
            sigf column label set by :obj: `AnomSearch`
        self._free : str
            free column label set by :obj: `AnomSearch`
        self._dano : str
            dano column label set by :obj: `AnomSearch`
        self._sigdano : str
            sigdano column label set by :obj: `AnomSearch`
        self._resolution : float
            mtz resolution set by :obj: `AnomSearch`

        Returns
        -------
        file
            mtz file containing FCalc, PHICalc, DANO and SIGDANO columns

        """
        cmd = ["cad", "HKLIN1", self.mtz,
               "HKLIN2", os.path.join(self.work_dir, "sfall_{0}.mtz".format(self.name)),
               "HKLOUT", os.path.join(self.work_dir, "cad_{0}.mtz".format(self.name))]
        stdin = os.linesep.join([
            "monitor BRIEF",
            "labin file 1 -",
            "    E1 = {0} -",
            "    E2 = {1} -",
            "    E3 = {2} -",
            "    E4 = {3} -",
            "    E5 = {4}",
            "labout file 1 -",
            "    E1 = {0} -",
            "    E2 = {1} -",
            "    E3 = {2} -",
            "    E4 = {3} -",
            "    E5 = {4}",
            "ctypin file 1 -",
            "    E1 = F -",
            "    E2 = Q -",
            "    E3 = I -",
            "    E4 = D -",
            "    E5 = Q",
            "resolution file 1 50 {5}",
            "labin file 2 -",
            "    E1 = FCalc -",
            "    E2 = PHICalc",
            "labout file 2 -",
            "    E1 = FCalc -",
            "    E2 = PHICalc",
            "ctypin file 2 -",
            "    E1 = F -",
            "    E2 = P",
        ])
        stdin = stdin.format(self._f, self._sigf, self._free, self._dano, self._sigdano, self._resolution)
        cexec(cmd, stdin=stdin)

    def fft(self):
        """Function to run FFT to create phased anomalous fourier map

        Parameters
        ----------
        self.name : str
            unique identifier for the input model set by :obj:`AnomSearch.run`
        self.work_dir : str
            working directory set by :obj:`AnomSearch.run`

        Returns
        -------
        file
            anomalous phased fourier map file
        file
            log file containing the results from the anomalous phased fourier

        """
        cmd = ["fft", "HKLIN", os.path.join(self.work_dir, "cad_{0}.mtz".format(self.name)),
               "MAPOUT", os.path.join(self.work_dir, "fft_{0}.map".format(self.name))]
        stdin = os.linesep.join(["xyzlim asu", "scale F1 1.0", "labin -",
                                 "    DANO={0} SIG1={1} PHI=PHICalc", "end"])
        stdin = stdin.format(self._dano, self._sigdano)
        cexec(cmd, stdin=stdin)

    def peakmax(self):
        """Function to run peakmax to return the peaks from FFT

        Parameters
        ----------
        self.name : str
            unique identifier for the input model set by :obj:`AnomSearch.run`
        self.work_dir : str
            working directory set by :obj:`AnomSearch.run`

        Returns
        -------
        file
            PDB file containing peaks
        file
            HA file containing peaks
        file
            log file containing the peaks identified by the anomalous phased fourier

        """
        cmd = ["peakmax", "MAPIN", os.path.join(self.work_dir, "fft_{0}.map".format(self.name)),
               "XYZOUT", os.path.join(self.work_dir, "peakmax_{0}.pdb".format(self.name)),
               "XYZFRC", os.path.join(self.work_dir, "peakmax_{0}.ha".format(self.name))]
        stdin = os.linesep.join(["numpeaks 8000", "THRESHOLD RMS 0", "output brookhaven frac"])
        cexec(cmd, stdin=stdin)

