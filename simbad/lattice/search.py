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
import pandas

logger = logging.getLogger(__name__)


class _LatticeParameterScore(object):
    """A basic lattice parameter scoring class"""

    __slots__ = ('pdb_code', 'unit_cell', 'total_penalty', 'length_penalty', 'angle_penalty')

    def __init__(self, pdb_code, unit_cell, total_penalty, length_penalty, angle_penalty):
        self.pdb_code = pdb_code
        self.unit_cell = unit_cell
        self.total_penalty = total_penalty
        self.length_penalty = length_penalty
        self.angle_penalty = angle_penalty

    def __repr__(self):
        return "{0}(pdb_code={1} unit_cell={2} total_penalty={3} length_penalty={4} angle_penalty={5}".format(
            self.__class__.__name__, self.pdb_code, self.unit_cell,
            self.total_penalty, self.length_penalty, self.angle_penalty
        )

    def _as_dict(self):
        """Convert the :obj:`_LatticeParameterScore <simbad.lattice.search._LatticeParameterScore>`
        object to a dictionary"""
        dict = {}
        for k in self.__slots__:
            if k == 'unit_cell':
                dict['a'], dict['b'], dict['c'], dict['alpha'], dict['beta'], dict['gamma'] = self.unit_cell
            else:
                dict[k] = getattr(self, k)
        return dict


class LatticeSearch(object):
    """A class to do a search for PDB entries with similar unit cell dimensions

    Attributes
    ----------
    unit_cell : str, list, tuple
       The unit cell parameters
    space_group : str
       The space group
    lattice_db_fname : str
       The path to the lattice database [pickle format]
    max_to_keep : int, optional
       The maximum number of results to keep [default: 20]
    ...

    Examples
    --------
    >>> from simbad.lattice.search import LatticeSearch
    >>> lattice_search = LatticeSearch('< a b c alpha beta gamma>', '<space group>', '<path to database>')
    >>> lattice_search.search()
    >>> lattice_search.summarize()

    If any results are found, the :func:`summarize` function will print a summary table to the screen
    and write the same information to a comma-separated file called ``lattice.csv``.

    """
    def __init__(self, unit_cell, space_group, lattice_db_fname, max_to_keep=20):
        """Initialize a new Lattice Search class

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
        if not self._search_results:
            return None
        return sorted(self._search_results, key=lambda x: float(x.total_penalty), reverse=False)[:self._max_to_keep]

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
        search_results = self.search_results
        if not search_results:
            msg = "No search results found/available"
            raise ValueError(msg)

        out_dir = os.path.join(dir, 'lattice_input_models')
        if not os.path.isdir(out_dir):
            os.mkdir(out_dir)
        else:
            msg = "Output directory exists: {0}".format(out_dir)
            raise ValueError(msg)

        for count, result in enumerate(search_results):
            if count <= self.max_to_keep:
                f_name = os.path.join(pdb_db, '{0}', 'pdb{1}.ent.gz').format(result.pdb_code[1:3], result.pdb_code)
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
        search_results = self.search_results
        if not search_results:
            msg = "No search results found/available"
            raise ValueError(msg)

        out_dir = os.path.join(dir, 'lattice_input_models')
        if not os.path.isdir(out_dir):
            os.mkdir(out_dir)
        else:
            msg = "Output directory exists: {0}".format(out_dir)
            raise ValueError(msg)

        pdbl = Bio.PDB.PDBList()
        for count, result in enumerate(search_results):
            if count < self.max_to_keep:
                # Download PDB file
                pdbl.retrieve_pdb_file(result.pdb_code, pdir=out_dir)
                # Rename the PDB file as appropriate
                os.rename(
                    os.path.join(out_dir, 'pdb{0}.ent'.format(result.pdb_code.lower())),
                    os.path.join(out_dir, '{0}.pdb'.format(result.pdb_code))
                )

    def search(self, tolerance=0.05):
        """Search for similar Niggli cells

        Parameters
        ----------
        tolerance : int, float, optional
           The tolerance applied for Niggli cell comparison [default: 0.05]

        """
        niggli_cell = numpy.asarray(self.niggli_cell)
        tol_niggli_cell = niggli_cell * tolerance

        results = []
        for i, db_cell in enumerate(self._lattice_db[1]):
            pdb_code = self._lattice_db[0][i]

            if LatticeSearch.cell_within_tolerance(niggli_cell, db_cell, tol_niggli_cell):
                total_pen, length_pen, angle_pen = LatticeSearch.calculate_penalty(niggli_cell, db_cell)
                score = _LatticeParameterScore(pdb_code, db_cell, total_pen, length_pen, angle_pen)
                results.append(score)

        self._search_results = results

    def summarize(self, csv_file='lattice.csv'):
        """Summarize the search results

        Parameters
        ----------
        csv_file : str
           The path for a backup CSV file

        Raises
        ------
        RuntimeError : No results found - search was unsuccessful

        """
        search_results = self.search_results
        if not search_results:
            msg = "No results found - search was unsuccessful"
            raise RuntimeError(msg)

        df = pandas.DataFrame(
            [r._as_dict() for r in search_results],
            index=[r.pdb_code for r in search_results],
            columns=['a', 'b', 'c', 'alpha', 'beta', 'gamma', 'length_penalty', 'angle_penalty', 'total_penalty'],
        )
        # Create a CSV for reading later
        df.to_csv(csv_file)
        # Display table in stdout
        summary_table = """
The lattice parameter search found the following structures:

%s
"""
        logger.info(summary_table % df.to_string())

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
        niggli_cell = numpy.array(eval(str(niggli_cell))).tolist()
        logger.info("Niggli cell calculated as: {0}".format(niggli_cell))
        return niggli_cell

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
            delta = abs(numpy.asarray(q, dtype=numpy.float64) - numpy.asarray(r, dtype=numpy.float64))
            return delta[:3].sum().item(), delta[3:].sum().item()

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
