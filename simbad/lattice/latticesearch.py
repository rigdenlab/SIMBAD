"""Module to skim the PDB for similar unit cells"""

from __future__ import division

__author__ = "Adam Simpkin & Felix Simkovic"
__date__ = "30 Jun 2017"
__version__ = "1.0"

import ast
import cctbx.crystal
import cctbx.uctbx
import datetime
import logging
import numpy
import os

from simbad.lattice.latticescore import LatticeSearchResult
from simbad.util import pdb_edit

logger = logging.getLogger(__name__)


class LatticeSearch(object):
    """A class to do a search for PDB entries with similar unit cell dimensions

    """
    def __init__(self, lattice_db_fname):
        """Initialize a new Lattice Search class

        Parameters
        ----------
        lattice_db_fname : str
           The path to the lattice database [pickle format]
        max_to_keep : int, optional
           The maximum number of results to keep [default: 20]

        """
        self._lattice_db_fname = None
        self.lattice_db_fname = lattice_db_fname

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

    def search(self, space_group, unit_cell, tolerance=0.05, max_to_keep=50, max_penalty=12):
        """Search for similar Niggli cells

        Parameters
        ----------
        unit_cell : list, tuple
           The parameters of the unit cell
        space_group : str
           The space group
        tolerance : int, float, optional
           The tolerance applied for Niggli cell comparison [default: 0.05]
        max_to_keep : int, optional
           The top-N number of results to return [default: 50]
        penalty_cut_off : int, optional
           The total penalty score over which results are ignored [default: 12]

        Returns
        -------
        list
           The results ordered by the total penalty score

        """
        space_group = LatticeSearch.check_sg(space_group)
        niggli_cell = LatticeSearch.calculate_niggli_cell(unit_cell, space_group)
        niggli_cell = numpy.asarray(niggli_cell)

        tol_niggli_cell = niggli_cell * tolerance

        results = []
        with numpy.load(self.lattice_db_fname) as compressed:
            for entry in compressed["arr_0"]:
                pdb_code = "".join(chr(c) for c in entry[:4].astype('uint8'))
                alt_cell = chr(int(entry[4])) if entry[4] != 0.0 else ' '
                db_cell = entry[5:]

                if self.cell_within_tolerance(niggli_cell, db_cell, tol_niggli_cell):
                    total_pen, length_pen, angle_pen = self.calculate_penalty(niggli_cell, db_cell)
                    vol_diff = self.calculate_volume_difference(niggli_cell, db_cell)
                    if total_pen < max_penalty:
                        prob = self.calculate_probability(total_pen)
                        score = LatticeSearchResult(pdb_code, alt_cell, db_cell, vol_diff, total_pen, length_pen,
                                                    angle_pen, prob)
                        results.append(score)

        results_sorted = sorted(results, key=lambda x: float(x.total_penalty), reverse=False)
        return results_sorted[:max_to_keep]
    
    @classmethod
    def calculate_penalty(cls, query, reference):
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

    @classmethod
    def calculate_probability(cls, penalty_score):
        """Calculate the probability that a penalty score will give a solution

        Parameters
        ----------
        penalty_score : float
            The total penalty score calculate for a search result

        Returns
        -------
        float
            Probability score
        """
        # Calculate probability score using exp equation calculated from test set
        probability = numpy.exp(-0.41208106 * penalty_score)

        # Set to 3dp
        probability = float("{0:.3}".format(probability))
        return probability

    @classmethod
    def cell_within_tolerance(cls, query, reference, tolerance):
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
    
    @classmethod
    def calculate_volume_difference(cls, query, reference):
        """Calculate the difference in volume between the query unit cell and the reference unit cell
        
        Parameters
        ----------
        query : list, tuple
           The query cell parameters
        reference : list, tuple
           The reference cell parameters
           
        Returns
        -------
        float
            The absolute difference in cell volumes
        """
        
        cell_volume_1 = cctbx.uctbx.unit_cell(query).volume()
        cell_volume_2 = cctbx.uctbx.unit_cell(reference).volume()
        
        difference = abs(cell_volume_1 - cell_volume_2)
        
        return float("{0:.3}".format(difference))
        
 
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
        unit_cell = cctbx.uctbx.unit_cell(unit_cell)
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
    def check_sg(sg):
        """Check the space group for known anomalies"""
        sg_conversion = {
            'A1': 'P1', 'B2': 'B112', 'C1211': 'C2', 'F422': 'I422',
            'I21': 'I2', 'I1211': 'I2', 'P21212A': 'P212121',
            'R3': 'R3:R', 'C4212': 'P422',
        }
        return sg_conversion.get(sg, sg)
 
    @staticmethod
    def copy_results(results, source, destination):
        """Copy the results from a local copy of the PDB

        Parameters
        ----------
        results: list
           A list of SIMBAD search results
        source : str
           The path to copy results from
        destination : str
           The path to save results to 

        Raises
        ------
        ValueError
           No search results found/available
        ValueError
           Output directory does not exist
        IOError
            Search result not found in installed PDB

        """
        if not results:
            msg = "No search results found/available"
            raise ValueError(msg)

        if not os.path.isdir(destination):
            msg = "Output directory does not exist: {0}".format(destination)
            raise ValueError(msg)

        import gzip

        to_del = []
        for count, result in enumerate(results):
            f_name = os.path.join(source, '{0}', 'pdb{1}.ent.gz').format(result.pdb_code[1:3].lower(),
                                                                         result.pdb_code.lower())
            f_name_out = os.path.join(destination, '{0}.pdb'.format(result.pdb_code))
            try:
                with gzip.open(f_name, 'rb') as f_in, open(f_name_out, 'w') as f_out:
                    f_out.write(f_in.read())
                pdb_edit.to_single_chain(f_name_out, f_name_out)
            except IOError:
                logger.warning("Encountered problem copying PDB %s from %s - removing entry from list",
                               result.pdb_code, source)
                to_del.append(count)

        # Remove any errors for the results data
        for i in reversed(to_del):
            results.pop(i)

    @staticmethod
    def download_results(results, destination):
        """Download the results directly from the PDB

        Parameters
        ----------
        results: list
           A list of SIMBAD search results
        destination : str
           The path to save results to 

        Raises
        ------
        ValueError
           No search results found/available
        ValueError
           Output directory does not exist
        RuntimeError
            Unable to download PDB

        """
        if results is None:
            msg = "No search results found/available"
            raise ValueError(msg)

        if not os.path.isdir(destination):
            msg = "Output directory does not exist: {0}".format(destination)
            raise ValueError(msg)
        
        import iotbx.pdb.fetch

        to_del = []
        for count, result in enumerate(results):
            try:
                content = iotbx.pdb.fetch.fetch(result.pdb_code, data_type='pdb', format='pdb', mirror='pdbe')
                logger.debug("Downloading PDB %s from %s", result.pdb_code, content.url)
                download_state = content.msg
            except RuntimeError:
                download_state = "FAIL"
        
            if download_state == "OK":
                f_name_out = os.path.join(destination, result.pdb_code + '.pdb')
                with open(f_name_out, 'w') as f_out:
                    f_out.write(content.read())
                pdb_edit.to_single_chain(f_name_out, f_name_out)

            elif download_state == "FAIL":
                logger.warning("Encountered problem downloading PDB %s - removing entry from list", result.pdb_code)
                to_del.append(count)
            else:
                logger.warning("Encountered problem downloading PDB %s from %s - removing entry from list",
                               result.pdb_code, content.url)
                to_del.append(count)
        
        # Remove any errors for the results data
        for i in reversed(to_del):
            results.pop(i)

    @staticmethod
    def summarize(results, csvfile=None):
        """Summarize the search results

        Parameters
        ----------
        results : list
           The lattice search results
        csvfile : str, optional
           The path to an output CSV file

        """
        from simbad.util import summarize_result
        columns = ['alt', 'a', 'b', 'c', 'alpha', 'beta', 'gamma', 'length_penalty', 'angle_penalty', 'total_penalty',
                   'volume_difference', 'probability_score']
        summarize_result(results, csv_file=csvfile, columns=columns)
    
