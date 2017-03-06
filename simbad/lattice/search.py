"""Module to skim the PDB for similar unit cells"""

from __future__ import division

__author__ = "Adam Simpkin & Felix Simkovic"
__date__ = "05 Mar 2017"
__version__ = "0.1"

import Bio.PDB
import cctbx.crystal
import cctbx.uctbx
import cPickle
import gzip
import logging
import numpy
import os

logger = logging.getLogger(__name__)


class _LatticeParameterScore(object):
    """A basic lattice parameter scoring class"""

    __slots__ = ('pdb_code', 'unit_cell', 'penalty_score', 'length_penalty', 'angle_penalty')

    def __init__(self, pdb_code, unit_cell, penalty_score, length_penalty, angle_penalty):
        self.pdb_code = pdb_code
        self.unit_cell = unit_cell
        self.penalty_score = penalty_score
        self.length_penalty = length_penalty
        self.angle_penalty = angle_penalty

    def __repr__(self):
        return "{0}(pdb_code={1} unit_cell={2} penalty_score={3} length_penalty={4} angle_penalty={5}".format(
            self.__class__.__name__, self.pdb_code, self.unit_cell,
            self.penalty_score, self.length_penalty, self.angle_penalty
        )

    def _for_csv(self):
        a, b, c, alpha, beta, gamma = self.unit_cell
        return ','.join(
            map(str,
                [self.pdb_code, a, b, c, alpha, beta, gamma,
                 self.penalty_score, self.length_penalty, self.angle_penalty]
            )
        )


class LatticeSearch(object):
    """A class to do a search for PDB entries with similar unit cell dimensions

    """
    def __init__(self, unit_cell, space_group, lattice_db_fname, max_to_keep=20):
        """

        Parameters
        ----------
        unit_cell : str, list, tuple
           The unit cell parameters
        space_group : str
           The space group
        lattice_db_fname : str
           The path to the lattice database [pickle format]
        max_to_keep : int, optional
           The maximum number of results to keep [default: 20]

        """

        self._unit_cell = None
        self._space_group = None

        self._niggli_cell = None
        self._lattice_db = None
        self._lattice_db_fname = None
        self._search_results = None
        self._max_to_keep = 0

        self.unit_cell = unit_cell
        self.space_group = space_group
        self.lattice_db_fname = lattice_db_fname
        self.max_to_keep = max_to_keep

    @property
    def lattice_db_fname(self):
        """The path to the lattice database"""
        return self._lattice_db_fname

    @lattice_db_fname.setter
    def lattice_db_fname(self, lattice_db_fname):
        """Define the lattice database filename"""
        self._lattice_db_fname = lattice_db_fname
        self._lattice_db = cPickle.load(open(lattice_db_fname))
        logger.info('Lattice database successfully cached')

    @property
    def max_to_keep(self):
        """The maximum number of results to keep"""
        return self._max_to_keep

    @max_to_keep.setter
    def max_to_keep(self, max_to_keep):
        """Define the maximum number of results to keep"""
        self._max_to_keep = int(max_to_keep)

    @property
    def niggli_cell(self):
        """The parameters of the Niggli cell"""
        return LatticeSearch.calculate_niggli_cell(self.unit_cell, self.space_group)

    @property
    def search_results(self):
        """The results from the lattice search"""
        return sorted(self._search_results, key=lambda x: float(x.penalty_score), reverse=False)[:self._max_to_keep]

    @property
    def space_group(self):
        """The space group"""
        return self._space_group

    @space_group.setter
    def space_group(self, space_group):
        """Define the space group"""
        # Check for known anomalies, may need to add to this
        if space_group == "B2":
            space_group = "B112"
        elif space_group == "C1211":
            space_group = "C2"
        elif space_group == "P21212A":
            space_group = "P212121"
        elif space_group == "R3":
            space_group = "R3:R"
        elif space_group == "C4212":
            space_group = "P422"
        self._space_group = space_group
        self._search_results = None

    @property
    def total_db_files(self):
        """The number of files in the lattice database"""
        return len(self._lattice_db[0])

    @property
    def unit_cell(self):
        """The unit cell parameters"""
        return self._unit_cell

    @unit_cell.setter
    def unit_cell(self, unit_cell):
        """Define the unit cell parameters"""
        if isinstance(unit_cell, str):
            unit_cell = unit_cell.split(' ')
        elif isinstance(unit_cell, list) or isinstance(unit_cell, tuple):
            unit_cell = unit_cell
        else:
            msg = "Unit cell parameters need to be of type str, list or tuple"
            raise TypeError(msg)
        self._unit_cell = unit_cell
        self._search_results = None

    def copy_results(self, pdb_db, dir=os.getcwd()):
        """Copy the results from a local copy of the PDB

        Parameters
        ----------
        pdb_db : str
           The path to the local copy of the PDB
        dir : str
           The path to save results to [default: .]

        Raises
        ------
        ValueError
           No search results found/available
        ValueError
           Output directory exists

        """
        if not self.search_results:
            msg = "No search results found/available"
            raise ValueError(msg)

        out_dir = os.path.join(dir, 'lattice_input_models')
        if not os.path.isdir(out_dir):
            os.mkdir(out_dir)
        else:
            msg = "Output directory exists: {0}".format(out_dir)
            raise ValueError(msg)

        for count, result in enumerate(self.search_results):
            if count <= self.max_to_keep:
                f_name = os.path.join(pdb_db, '{0}', 'pdb{1}.ent.gz'.format(result.pdb_code[1:3], result.pdb_code))
                with gzip.open(f_name, 'rb') as f_in:
                    f_name_out = os.path.join(out_dir, '{0}.pdb'.format(result.pdb_code))
                    with open(f_name_out, 'w') as f_out:
                        f_out.write(f_in.read())

    def download_results(self, dir=os.getcwd()):
        """Download the results directly from the PDB

        Parameters
        ----------
        dir : str
           The path to save results to [default: .]

        Raises
        ------
        ValueError
           No search results found/available
        ValueError
           Output directory exists

        """
        if not self.search_results:
            msg = "No search results found/available"
            raise ValueError(msg)

        out_dir = os.path.join(dir, 'lattice_input_models')
        if not os.path.isdir(out_dir):
            os.mkdir(out_dir)
        else:
            msg = "Output directory exists: {0}".format(out_dir)
            raise ValueError(msg)

        pdbl = Bio.PDB.PDBList()
        for count, result in enumerate(self.search_results):
            if count < self.max_to_keep:
                # Download PDB file
                pdbl.retrieve_pdb_file(result.pdb_code, pdir=out_dir)
                # Rename the PDB file as appropriate
                os.rename(
                    os.path.join(out_dir, 'pdb{0}.ent'.format(result.pdb_code.lower())),
                    os.path.join(out_dir, '{0}.pdb'.format(result.pdb_code))
                )

    def search(self):
        """Search for similar Niggli cells"""
        niggli_cell = numpy.asarray(self.niggli_cell)
        tol_niggli_cell = niggli_cell / 100 * 5

        results = []
        for i, db_cell in enumerate(self._lattice_db[1]):
            pdb_code = self._lattice_db[0][i]

            if LatticeSearch.cell_within_tolerance(niggli_cell, db_cell, tol_niggli_cell):
                total_pen, length_pen, angle_pen = LatticeSearch.calculate_penalty(niggli_cell, db_cell)
                score = _LatticeParameterScore(pdb_code, db_cell, total_pen, length_pen, angle_pen)
                results.append(score)

        results = sorted(results, key=lambda x: float(x.penalty_score), reverse=False)

        self._search_results = results

    def summarize(self):
        """Summarize the search results

        Raises
        ------
        RuntimeError : No results found - search was unsuccessful

        """
        search_results = self._search_results
        if not search_results:
            msg = "No results found - lattice search was unsuccessful"
            raise RuntimeError(msg)

        header=['PDB_CODE', 'A', 'B', 'C', 'ALPHA', 'BETA', 'GAMMA', 'LENGTH_PENALTY', 'ANGLE_PENALTY', 'TOTAL_PENALTY']
        matrix = []
        pdb_codes = []
        for result in search_results:
            # Create a CSV for reading later
            with open('lattice.csv', 'a') as f:
                f.write(result._for_csv() + os.linesep)

            list = []
            pdb_codes.append(result.pdb_code)
            list.append(result.unit_cell[0])
            list.append(result.unit_cell[1])
            list.append(result.unit_cell[2])
            list.append(result.unit_cell[3])
            list.append(result.unit_cell[4])
            list.append(result.unit_cell[5])
            list.append(result.length_penalty)
            list.append(result.angle_penalty)
            list.append(result.penalty_score)
            matrix.append(list)

        summary_table = self.format_matrix(header, pdb_codes, matrix, '{:^{}}', '{:<{}}', '{:>{}.3f}', '\n', ' | ')

        logger.info("The lattice parameter search found the following structures:")

        logger.info(summary_table)

    @staticmethod
    def format_matrix(header, pdb_codes, matrix, top_format, left_format, cell_format, row_delim, col_delim):
        """Code to format output of search"""
        table = [header] + [[name] + row for name, row in zip(pdb_codes, matrix)]
        table_format = [['{:^{}}'] + len(header) * [top_format]] \
                       + len(matrix) * [[left_format] + len(header) * [cell_format]]

        col_widths = [max(
            len(format.format(cell, 0))
            for format, cell in zip(col_format, col))
            for col_format, col in zip(zip(*table_format), zip(*table))]
        return row_delim.join(
            col_delim.join(
                format.format(cell, width)
                for format, cell, width in zip(row_format, row, col_widths))
            for row_format, row in zip(table_format, table))

    @staticmethod
    def calculate_niggli_cell(unit_cell, space_group):
        """Calculate the parameters of the Niggli cell

        Parameters
        ----------
        unit_cell : list, tuple
           The parameters of the unit cell
        space_group : str
           The space group

        Returns
        -------
        list
           The Niggli cell parameters

        """
        unit_cell = cctbx.uctbx.unit_cell(' '.join(unit_cell))

        xs = cctbx.crystal.symmetry(
            unit_cell=unit_cell,
            space_group=space_group
        )
        niggli_cell = xs.change_basis(xs.change_of_basis_op_to_niggli_cell()).unit_cell()
        niggli_cell = numpy.array(eval(str(niggli_cell)))
        logger.info("Niggli cell calculated as: {0}".format(niggli_cell))
        return niggli_cell.tolist()

    @staticmethod
    def calculate_penalty(query, reference):
        """Calculate the linear cell variation between unit cells

        Parameters
        ----------
        query : list, tuple
           The query cell parameters
        reference : list, tuple
           The reference cell parameters

        Returns
        -------

        """
        def penalty(q, r):
            delta = abs(numpy.asarray(q) - numpy.asarray(r))
            return delta[:3].sum().item(), delta[3:].sum().astype(int).item()

        length_penalty, angle_penalty = penalty(query, reference)
        total_penalty = length_penalty + angle_penalty
        return total_penalty, length_penalty, angle_penalty

    @staticmethod
    def cell_within_tolerance(query, reference, tolerance):
        """Compare two cells and determine if ``query`` is within ``reference`` cell parameter tolerance

        Parameters
        ----------
        query : list, tuple
           The query cell parameters
        reference : list, tuple
           The reference cell parameters
        tolerance : list, tuple
           The tolerance cell parameter values

        Returns
        -------
        bool

        """
        for q, r, t in zip(query, reference, tolerance):
            if numpy.allclose(q, r, atol=t):
                continue
            else:
                return False
        return True
