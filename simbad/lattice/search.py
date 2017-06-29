"""Module to skim the PDB for similar unit cells"""

from __future__ import division

__author__ = "Adam Simpkin & Felix Simkovic"
__date__ = "05 Mar 2017"
__version__ = "0.1"

import argparse
import ast
import cctbx.crystal
import cctbx.uctbx
import datetime
import iotbx.pdb.fetch
import gzip
import logging
import numpy
import os
import pandas

from pyjob.misc import StopWatch

import simbad
import simbad.command_line
import simbad.exit
import simbad.util.mtz_util

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
        dictionary = {}
        for k in self.__slots__:
            if k == 'unit_cell':
                for k, v in zip(['a', 'b', 'c', 'alpha', 'beta', 'gamma'], self.unit_cell):
                    dictionary[k] = v
            else:
                dictionary[k] = getattr(self, k)
        return dictionary


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
        # Check how old the database is, if older than 90 days print
        # message to suggest updating
        timestamp = os.path.getmtime(lattice_db_fname)
        if (datetime.date.today() - datetime.date.fromtimestamp(timestamp)).days > 90:
            logger.info('Lattice database is older than 90 days, consider updating!\n'
                        'Use the "simbad-create-lattice-db" script in your Terminal')
        self._lattice_db_fname = lattice_db_fname

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
        if self._search_results is None:
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
        sg_conversion = {
            'A1': 'P1', 'B2': 'B112', 'C1211': 'C2', 'F422': 'I422',
            'I21': 'I2', 'I1211': 'I2', 'P21212A': 'P212121',
            'R3': 'R3:R', 'C4212': 'P422',
        }

        space_group = sg_conversion.get(space_group, space_group)
        self._space_group = space_group
        self._search_results = None

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

    def copy_results(self, pdb_db, directory=os.getcwd()):
        """Copy the results from a local copy of the PDB

        Parameters
        ----------
        pdb_db : str
           The path to the local copy of the PDB
        directory : str
           The path to save results to [default: .]

        Raises
        ------
        ValueError
           No search results found/available
        ValueError
           Output directory does not exist
        IOError
            Search result not found in installed PDB
        """
        search_results = self.search_results
        if search_results is None:
            msg = "No search results found/available"
            raise ValueError(msg)

        if not os.path.isdir(directory):
            msg = "Output directory does not exist: {0}".format(directory)
            raise ValueError(msg)

        to_del = []
        for count, result in enumerate(search_results):
            if count <= self.max_to_keep:
                try:
                    f_name = os.path.join(pdb_db, '{0}', 'pdb{1}.ent.gz').format(result.pdb_code[1:3].lower(),
                                                                                 result.pdb_code.lower())
                    with gzip.open(f_name, 'rb') as f_in:
                        f_name_out = os.path.join(directory, '{0}.pdb'.format(result.pdb_code))
                        with open(f_name_out, 'w') as f_out:
                            f_out.write(f_in.read())
                except IOError:
                    logger.warning("Encountered problem copying PDB %s from %s - removing entry from list",
                                   result.pdb_code, pdb_db)
                    to_del.append(count)

        # Remove any errors for the search_results data
        for i in reversed(to_del):
            search_results.pop(i)
        self._search_results = search_results

    def download_results(self, directory=os.getcwd()):
        """Download the results directly from the PDB

        Parameters
        ----------
        directory : str
           The path to save results to [default: .]

        Raises
        ------
        ValueError
           No search results found/available
        ValueError
           Output directory does not exist
        RuntimeError
            Unable to download PDB
        """
        search_results = self.search_results
        if search_results is None:
            msg = "No search results found/available"
            raise ValueError(msg)

        if not os.path.isdir(directory):
            msg = "Output directory does not exist: {0}".format(directory)
            raise ValueError(msg)

        to_del = []
        # Deliberately take full results to bypass max_to_keep check
        for count, result in enumerate(search_results):
            if count < self.max_to_keep + len(to_del):
                try:
                    content = iotbx.pdb.fetch.fetch(result.pdb_code, data_type='pdb', format='pdb', mirror='pdbe')
                    download_state = content.msg
                    logger.debug("Downloading PDB %s from %s", result.pdb_code, content.url)
                except RuntimeError:
                    download_state = "FAIL"
                    
                if download_state == "OK":
                    with open(os.path.join(directory, result.pdb_code + '.pdb'), 'w') as f_out:
                        f_out.write(content.read())
                elif download_state == "FAIL":
                    logger.warning("Encountered problem downloading PDB %s - removing entry from list", result.pdb_code)
                    to_del.append(count)
                else:
                    logger.warning("Encountered problem downloading PDB %s from %s - removing entry from list", result.pdb_code, content.url)
                    to_del.append(count)
        
        # Remove any errors for the search_results data
        for i in reversed(to_del):
            search_results.pop(i)
        self._search_results = search_results

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
        with numpy.load(self.lattice_db_fname) as compressed:
            for entry in compressed['arr_0']:
                pdb_code = "".join(chr(c) for c in entry[:4].astype('uint8'))
                db_cell = entry[4:]

                if LatticeSearch.cell_within_tolerance(niggli_cell, db_cell, tol_niggli_cell):
                    total_pen, length_pen, angle_pen = LatticeSearch.calculate_penalty(niggli_cell, db_cell)
                    score = _LatticeParameterScore(pdb_code, db_cell, total_pen, length_pen, angle_pen)
                    results.append(score)

        self._search_results = results

    def summarize(self, csv_file=None):
        """Summarize the search results

        Parameters
        ----------
        csv_file : str
           The path for a backup CSV file

        Raises
        ------
        RuntimeError
           No results found - search was unsuccessful

        """
        search_results = self.search_results
        if search_results is None:
            msg = "No results found - search was unsuccessful"
            raise RuntimeError(msg)

        df = pandas.DataFrame(
            [r._as_dict() for r in search_results],
            index=[r.pdb_code for r in search_results],
            columns=['a', 'b', 'c', 'alpha', 'beta', 'gamma', 'length_penalty', 'angle_penalty', 'total_penalty'],
        )
        
        # Create a CSV for reading later
        if csv_file:
            df.to_csv(csv_file)

        if df.empty:
            logger.info("The lattice parameter search found no matching structures")
        else:
            # Display table in stdout
            summary_table = """
The lattice parameter search found the following structures:

%s
"""
            logger.info(summary_table, df.to_string())

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
            space_group=space_group,
            correct_rhombohedral_setting_if_necessary=True
        )
        niggli_cell = xs.change_basis(xs.change_of_basis_op_to_niggli_cell()).unit_cell()
        niggli_cell = numpy.array(ast.literal_eval(str(niggli_cell))).tolist()
        logger.info("Niggli cell calculated as: [%s]", ", ".join(map(str, niggli_cell)))
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
        float
           Total penalty
        float
           Length penalty
        float
           Angle penalty
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
    
def main():
    """SIMBAD lattice search function"""
    p = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    sg = p.add_argument_group('Lattice search specific options')
    sg.add_argument('-mtz', type=str, default=None,
                    help="Path to an input MTZ file")
    sg.add_argument('-cell', type=str, default=None,
                    help="Input unit cell instead of MTZ file, format 'a,b,c,alpha,beta,gamma'")
    sg.add_argument('-space_group', type=str, default=None,
                    help="Input space group instead of MTZ file")
    sg.add_argument('-latt_db', type=str, default=simbad.LATTICE_DB,
                    help='Path to local copy of the lattice database')
    sg.add_argument('-max_to_keep', type=int, default=10,
                    help='Number of results to display')
    sg.add_argument('-debug_lvl', type=str, default='info',
                   help='The console verbosity level < notset | info | debug | warning | error | critical > ')
    args = p.parse_args()
    
    # Logger setup
    logger = simbad.command_line.setup_logging(level=args.debug_lvl)
    
    # Print a fancy header
    simbad.command_line.print_header()

    # Start taking time
    stopwatch = StopWatch()
    stopwatch.start()
    
    # Run lattice parameter search
    logger.info("Running lattice parameter search")
    
    if args.cell and args.space_group:
        cell_parameters, space_group = args.cell, args.space_group
        logger.debug("Using input cell parameters: {0} and space_group: {1}".format(cell_parameters,
                                                                           space_group))
    elif args.mtz:
        space_group, _, cell_parameters = simbad.util.mtz_util.crystal_data(args.mtz)
        logger.debug("Using cell parameters: {0} and space_group: {1} from MTZ: {2}".format(cell_parameters,
                                                                                            space_group,
                                                                                            args.mtz))
    else:
        msg = "Cell parameters and space group not provided"
        logger.debug(msg)
        raise RuntimeError(msg)
        
    lattice_search = LatticeSearch(
        cell_parameters, space_group, args.latt_db, max_to_keep=args.max_to_keep
    )
    lattice_search.search()
    if lattice_search.search_results:
        lattice_search.summarize()
    
    # Calculate and display the runtime in hours
    stopwatch.stop()
    logger.info("All processing completed in %d days, %d hours, %d minutes, and %d seconds", *stopwatch.time_pretty)

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        msg = "Error running main SIMBAD program: {0}".format(e.message)
        simbad.exit.exit_error(*sys.exc_info())
    
    
    
