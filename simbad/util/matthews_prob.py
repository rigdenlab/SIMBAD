"""Module for running matthews probabilities to calculate solvent content"""

__author__ = "Adam Simpkin"
__date__ = "3 June 2020"
__version__ = "1.0"

import abc
import math
import numpy as np

from simbad.util.pdb_util import PdbStructure

ABC = abc.ABCMeta('ABC', (object,), {})


class _MatthewsCoefficient(ABC):
    def __init__(self, cell_volume):
        self.cell_volume = cell_volume

    @abc.abstractmethod
    def calculate_from_file(self, pdb):
        """ Abstract method to calculate Matthews Coefficient from input PDB"""
        pass

    @abc.abstractmethod
    def calculate_from_struct(self, struct):
        """ Abstract method to calculate Matthews Coefficient from PDB util :obj:"""
        pass

    def get_macromolecule_fraction(self, vm):
        """Calculate the macromolecule fraction"""
        return 1. / (6.02214e23 * 1e-24 * 1.35 * vm)


class SolventContent(_MatthewsCoefficient):
    def __init__(self, cell_volume):
        super(SolventContent, self).__init__(cell_volume)

    def calculate_from_file(self, pdb):
        struct = PdbStructure.from_file(pdb)
        return self.calculate_from_struct(struct)

    def calculate_from_struct(self, struct):
        return self._calculate(struct.molecular_weight)

    def _calculate(self, mw):
        if mw <= 0:
            raise ValueError("Incorrect Molecular Weight")
        vm = self.cell_volume / mw
        macromolecule_fraction = self.get_macromolecule_fraction(vm)
        solvent_fraction = 1.0 - macromolecule_fraction
        return solvent_fraction * 100


class MatthewsProbability(_MatthewsCoefficient):
    def __init__(self, cell_volume):
        super(MatthewsProbability, self).__init__(cell_volume)

    def calculate_from_file(self, pdb):
        struct = PdbStructure.from_file(pdb)
        return self.calculate_from_struct(struct)

    def calculate_from_struct(self, struct):
        return self._calculate(struct.molecular_weight)

    def _calculate(self, mw):
        n_copies = 0
        solvent_fraction = 1.0
        scores = []
        while solvent_fraction > 0:
            n_copies += 1
            vm = self.cell_volume / (mw * n_copies)
            macromolecule_fraction = self.get_macromolecule_fraction(vm)
            if macromolecule_fraction > 1:
                break
            solvent_fraction = 1.0 - macromolecule_fraction
            probability = self._calculate_solvent_probability(solvent_fraction)
            scores.append((n_copies, solvent_fraction, probability))
        return self._get_max_score(scores)

    def _calculate_solvent_probability(self, solvent):
        """Calculate solvent fraction probability"""
        assert 0 < solvent < 1

        coeffs = [-14.105436736742137,
                  -0.47015366358636385,
                  -2.9151681976244639,
                  -0.49308859741473005,
                  0.90132625209729045,
                  0.033529051311488103,
                  0.088901407582105796,
                  0.10749856607909694,
                  0.055000918494099861,
                  -0.052424473641668454,
                  -0.045698882840119227,
                  0.076048484096718036,
                  -0.097645159906868589,
                  0.03904454313991608,
                  -0.072186667173865071]

        chebyshev_poly = np.polynomial.Chebyshev(coeffs, domain=[0, 1])
        return math.exp(chebyshev_poly(solvent))

    def _get_max_score(self, scores):
        """Use the probability score to guess the number of copies and solvent content"""
        max = 0.0
        n_copies = 1
        solvent_fraction = 0.5
        for i in range(len(scores)):
            prob = scores[i][-1]
            if prob > max:
                max = prob
                n_copies = scores[i][0]
                solvent_fraction = scores[i][1]
        return solvent_fraction, n_copies

