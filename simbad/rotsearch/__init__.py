"""Module to run the rotation search"""

__author__ = "Adam Simpkin & Felix Simkovic"
__date__ = "03 June 2020"
__version__ = "0.4"

import abc
import simbad.parsers.mtz_parser

ABC = abc.ABCMeta('ABC', (object,), {})


class _RotationSearch(ABC):
    def __init__(self, mtz, mr_program, tmp_dir, work_dir, max_to_keep=20, skip_mr=False, process_all=False, **kwargs):
        self.max_to_keep = max_to_keep
        self.mr_program = mr_program
        self.mtz = mtz
        self.mtz_obj = simbad.parsers.mtz_parser.MtzParser(mtz)
        self.mtz_obj.parse()
        self.skip_mr = skip_mr
        self.process_all = process_all
        self.tmp_dir = tmp_dir
        self.work_dir = work_dir

        self.columns = []
        self.score_column = None
        self.simbad_dat_files = None
        self.solution = False
        self.submit_qtype = None
        self.submit_queue = None
        self._search_results = None
        self.tested = []

    # ------------------ Abstract methods and properties ------------------

    @abc.abstractmethod
    def run(self, models_dir, **kwargs):
        """ Abstract method to run rotation function"""
        pass

    @abc.abstractmethod
    def generate_script(self, dat_model):
        """ Abstract method to generate run script"""
        pass

    @abc.abstractmethod
    def rot_succeeded_log(self, log):
        """ Abstract method to check if rotation function has succeeded from the log"""
        pass

    # ------------------ Some general methods ------------------

    @property
    def search_results(self):
        return sorted(self._search_results,
                      key=lambda x: float(x[x._fields.index(self.score_column)]), reverse=True)[: self.max_to_keep]

    def summarize(self, csv_file):
        """Summarize the search results
        Parameters
        ----------
        csv_file : str
            The path for a backup CSV file
        columns : list
            List containing column labels
        Raises
        ------
            No results found
        """
        from simbad.util import summarize_result

        summarize_result(self.search_results, csv_file=csv_file, columns=self.columns)


def rotation_search_factory(method):
    if method == "amore":
        from simbad.rotsearch.amore_search import AmoreRotationSearch

        return AmoreRotationSearch
    elif method == "phaser":
        from simbad.rotsearch.phaser_search import PhaserRotationSearch

        return PhaserRotationSearch
    else:
        raise ValueError("Unrecognised program entered to perform the rotation search: %s", method)


def get_chunk_size(total, size):
    return total if size == 0 else size


def get_total_chunk_cycles(total, step):
    total_chunk_cycles, remainder = divmod(total, step)
    if remainder > 0:
        return total_chunk_cycles + 1
    else:
        return total_chunk_cycles
