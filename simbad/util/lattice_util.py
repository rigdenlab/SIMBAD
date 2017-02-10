"""
Class to skim the PDB for similar unit cells
@author: hlasimpk
"""

# Imports
import os
import numpy as np
from collections import namedtuple
import cctbx.crystal
import cctbx.uctbx
import cPickle
import gzip
import logging

from simbad.constants import SIMBAD_LATTICE_DB

logger = logging.getLogger(__name__)

class Lattice_search(object):
    def __init__(self, optd=None):

        self.unit_cell      = None
        self.space_group    = None
        self.database       = []
        self.niggli_cell    = None
        self.results        = []
        self.pickle_file    = SIMBAD_LATTICE_DB
        self.total_files    = 0

        if optd:
            self.init(optd)

        if self.unit_cell and self.space_group:
            self.calculate_niggli_cell()
            logger.info("Niggli cell calculated as: {0}".format(self.niggli_cell))

        if not self.database:
            self.load_database()

        if self.database:
            Lattice_parameter_scores = namedtuple("Lattice_parameter_scores", ['PDB_code',
                                                                               'unit_cell',
                                                                               'Penalty_score',
                                                                               'Length_penalty',
                                                                               'Angle_penalty'])

            item_ct = 0
            lattice_pars = self.database[1]
            for db_cell in lattice_pars:
                pdb_code = self.database[0][item_ct]
                try:
                    # Calculate the tolerances for the input niggli cell once
                    self.tolerance_a = ((self.niggli_cell[0] / 100) * 5)
                    self.tolerance_b = ((self.niggli_cell[1] / 100) * 5)
                    self.tolerance_c = ((self.niggli_cell[2] / 100) * 5)
                    self.tolerance_alpha = ((self.niggli_cell[3] / 100) * 5)
                    self.tolerance_beta = ((self.niggli_cell[4] / 100) * 5)
                    self.tolerance_gamma = ((self.niggli_cell[5] / 100) * 5)
                    within_tol = self.cell_tolerance(db_cell)
                except ValueError:
                    within_tol = False

                if within_tol:
                    total_penalty, length_penalty, angle_penalty = self.calculate_penalty(db_cell)
                    a, b, c, alpha, beta, gamma = db_cell
                    self.results.append(Lattice_parameter_scores(PDB_code=pdb_code,
                                                                 unit_cell='{0},{1},{2},{3},{4},{5}'.format(a, b, c,
                                                                                                            alpha, beta,
                                                                                                            gamma),
                                                                 Penalty_score="%08.5f" % float(abs(total_penalty)),
                                                                 Length_penalty="%08.5f" % float(abs(length_penalty)),
                                                                 Angle_penalty="%08.5f" % float(abs(angle_penalty))))
                item_ct += 1


        if self.results:
            self.rank_scores()
            self.return_results()
            if not optd['PDB']:
                logger.info("Downloading structures of the top 20 results")
                #self.download_pdbs()
            else:
                logger.info("Obtaining structures for the top 20 results from: {0}".format(optd['PDB']))

        return

    def init(self, optd):
        '''Set input arguments as class variables'''

        if 'cell_parameters' in optd.keys() and optd['cell_parameters']:
            self.unit_cell = optd['cell_parameters']
        else:
            msg = "Cell parameters not specified"
            logger.critical(msg)
            raise RuntimeError(msg)

        if 'space_group' in optd.keys() and optd['space_group']:
            self.space_group = optd['space_group']
        else:
            msg = "Space group not specified"
            logger.critical(msg)
            raise RuntimeError(msg)

        return

    def calculate_niggli_cell(self):
        '''Function to calculate the niggli cell'''

        # Check for known anomalies, may need to add to this
        if self.space_group == "B2":
            self.space_group = "B112"
        elif self.space_group == "C1211":
            self.space_group = "C2"
        elif self.space_group == "P21212A":
            self.space_group = "P212121"
        elif self.space_group == "R3":
            self.space_group = "R3:R"
        elif self.space_group == "C4212":
            self.space_group = "P422"

        unit_cell = cctbx.uctbx.unit_cell(self.unit_cell)

        xs = cctbx.crystal.symmetry(
            unit_cell=unit_cell,
            space_group=self.space_group
        )
        niggli_cell = xs.change_basis(xs.change_of_basis_op_to_niggli_cell()).unit_cell()
        self.niggli_cell = np.array(eval(str(niggli_cell)))
        return

    def load_database(self):
        '''Function to load the database'''

        self.database = cPickle.load(open(self.pickle_file))
        self.total_files = len(self.database[0])
        logger.info("  - database loaded from pickle file")

        return

    def calculate_penalty(self, ref_niggli_cell):
        '''Calculate the linear cell variation between unit cells'''
        # Try needed to account for ref niggli cells that didn't work
        try:
            length_penalty, angle_penalty = self.penalty(ref_niggli_cell)

            total_penalty = length_penalty + angle_penalty

            return total_penalty, length_penalty, angle_penalty
        except ValueError:
            return None

    def cell_tolerance(self, ref_niggli_cell):
        '''Function to remove cells outside certain tolerances from being considered'''

        # Calculate if the 6 lattice parameters are within tolerances
        a = np.allclose(self.niggli_cell[0], ref_niggli_cell[0], atol=self.tolerance_a)
        b = np.allclose(self.niggli_cell[1], ref_niggli_cell[1], atol=self.tolerance_b)
        c = np.allclose(self.niggli_cell[2], ref_niggli_cell[2], atol=self.tolerance_c)
        alpha = np.allclose(self.niggli_cell[3], ref_niggli_cell[3], atol=self.tolerance_alpha)
        beta = np.allclose(self.niggli_cell[4], ref_niggli_cell[4], atol=self.tolerance_beta)
        gamma = np.allclose(self.niggli_cell[5], ref_niggli_cell[5], atol=self.tolerance_gamma)

        within_tol = False
        if a and b and c and alpha and beta and gamma:
            within_tol = True

        return within_tol

    def penalty(self, input_cell):
        '''Function to add penalty score depending on difference between unit cells'''
        # Split the experimental niggli cell into into component parts
        a_e, b_e, c_e, alpha_e, beta_e, gamma_e = self.niggli_cell

        # Split the input niggli cell into into component parts
        a_r, b_r, c_r, alpha_r, beta_r, gamma_r = input_cell

        delta_a = abs(a_e - a_r)
        delta_b = abs(b_e - b_r)
        delta_c = abs(c_e - c_r)
        delta_alpha = abs(alpha_e - alpha_r)
        delta_beta = abs(beta_e - beta_r)
        delta_gamma = abs(gamma_e - gamma_r)

        lengths_penalty_score = delta_a + delta_b + delta_c
        angles_penalty_score = delta_alpha + delta_beta + delta_gamma
        return lengths_penalty_score, angles_penalty_score

    def rank_scores(self):
        '''Rank the scores in self.results by Penalty score'''
        self.results = sorted(self.results, key=lambda x: float(x.Penalty_score), reverse=False)
        return

    def return_results(self):
        '''Log the results from the lattice parameter search'''
        logger.info("The lattice parameter search found the following structures:")
        logger.info("PDB_CODE |   A   |   B   |   C   | ALPHA |  BETA | GAMMA | LENGTH_PENALTY | ANGLE_PENALTY | TOTAL_PENALTY ")
        for result in self.results:
            uc = result.unit_cell.split(',')
            logger.info('  {0}   | {1} | {2} | {3} | {4} | {5} | {6} |    {7}    |    {8}    |    {9}'.format(result.PDB_code,
                                                                                                              ("{0:.2f}".format(float(uc[0]))),
                                                                                                              ("{0:.2f}".format(float(uc[1]))),
                                                                                                              ("{0:.2f}".format(float(uc[2]))),
                                                                                                              ("{0:.2f}".format(float(uc[3]))),
                                                                                                              ("{0:.2f}".format(float(uc[4]))),
                                                                                                              ("{0:.2f}".format(float(uc[5]))),
                                                                                                              result.Length_penalty,
                                                                                                              result.Angle_penalty,
                                                                                                              result.Penalty_score,))

        return

    def download_pdbs(self):
        '''Copy across the top 20 PDB files identified by lattice parameter search'''
        with open('lattice_parameter_search.txt', 'r') as f:
            for line in f:
                PDB_code = line.split(',')[0].lower()
                with gzip.open('/data3/opt/database/pdb/{0}/pdb{1}.ent.gz'.format(PDB_code[1:3], PDB_code), 'rb') as f:
                    file_content = f.read()
                    with open(os.path.join(os.getcwd(), 'lattice_input_models/{0}.pdb'.format(PDB_code.upper())), 'w') as o:
                        o.write(file_content)
        return






