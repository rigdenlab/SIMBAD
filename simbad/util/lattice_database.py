"""Script to update the database for the lattice parameter search"""

__author__ = "Adam Simpkin"
__date__ = "28 Apr 2017"
__version__ = "0.1"

import cctbx.crystal
import cctbx.uctbx
import cPickle
import numpy
import os
import urllib


class _LatticeParameters(object):
    """A basic lattice parameter sorting class"""

    __slots__ = ('pdb_code', 'A', 'B', 'C', 'ALPHA', 'BETA', 'GAMMA', 'space_group')

    def __init__(self, pdb_code, A, B, C, ALPHA, BETA, GAMMA, space_group):
        self.pdb_code = pdb_code
        self.A = A
        self.B = B
        self.C = C
        self.ALPHA = ALPHA
        self.BETA = BETA
        self.GAMMA = GAMMA
        self.space_group = space_group

    def __repr__(self):
        return "{0}(pdb_code={1} A={2} B={3} C={4} ALPHA={5}, BETA={6}, GAMMA={7}, space_group={8}".format(
            self.__class__.__name__, self.pdb_code, self.A, self.B, self.C,
            self.ALPHA, self.BETA, self.GAMMA, self.space_group)

    def _as_dict(self):
        """Convert the :obj:`_LatticeParameterScore <simbad.lattice.search._LatticeParameterScore>`
        object to a dictionary"""
        dictionary = {}
        for k in self.__slots__:
            dictionary[k] = getattr(self, k)
        return dictionary


def custom_report():
    link = 'http://www.rcsb.org/pdb/rest/customReport.csv?pdbids=*&customReportColumns=lengthOfUnitCellLatticeA,' \
           'lengthOfUnitCellLatticeB,lengthOfUnitCellLatticeC,unitCellAngleAlpha,unitCellAngleBeta,unitCellAngleGamma,' \
           'spaceGroup&format=csv&service=wsfile'

    urllib.urlretrieve(link, os.path.join(os.getcwd(), 'results.csv'))

    results = []
    with open(os.path.join(os.getcwd(), 'results.csv'), 'r') as f:
        for line in f:
            if line.startswith('structureId'):
                continue
            else:
                try:
                    tmp = line.split('","')
                    pdb_code = tmp[0][1:]
                    A = float(tmp[1])
                    B = float(tmp[2])
                    C = float(tmp[3])
                    ALPHA = float(tmp[4])
                    BETA = float(tmp[5])
                    GAMMA = float(tmp[6])
                    space_group = tmp[7][:-2]

                    data = _LatticeParameters(pdb_code, A, B, C, ALPHA, BETA, GAMMA, space_group)
                    results.append(data)
                except ValueError:
                    pass

    os.unlink(os.path.join(os.getcwd(), 'results.csv'))

    return results


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
        space_group=space_group
    )
    niggli_cell = xs.change_basis(xs.change_of_basis_op_to_niggli_cell()).unit_cell()
    niggli_cell = numpy.array(eval(str(niggli_cell))).tolist()
    return niggli_cell


def create_niggli_cell_db(crystal_data):
    pickle_file = 'niggli_database.cpk'

    names = []
    lattice_pars = []

    for i, c in enumerate(crystal_data):
        unit_cell = ' '.join(str(p) for p in [c.A, c.B, c.C, c.ALPHA, c.BETA, c.GAMMA])

        if c.space_group.replace(' ', '') == "A1":
            space_group = "P1"
        elif c.space_group.replace(' ', '') == "B2":
            space_group = "B112"
        elif c.space_group.replace(' ', '') == "C1211":
            space_group = "C2"
        elif c.space_group.replace(' ', '') == "F422":
            space_group = "I422"
        elif c.space_group.replace(' ', '') == "I21":
            space_group = "I2"
        elif c.space_group.replace(' ', '') == "I1211":
            space_group = "I2"
        elif c.space_group.replace(' ', '') == "P21212A":
            space_group = "P212121"
        elif c.space_group.replace(' ', '') == "R3":
            space_group = "R3:R"
        elif c.space_group.replace(' ', '') == "C4212":
            space_group = "P422"
        else:
            space_group = c.space_group

        try:
            niggli_cell = calculate_niggli_cell(unit_cell, space_group)
        except AssertionError:
            pass
        except ValueError:
            pass
        except RuntimeError:
            pass

        if niggli_cell:
            vals = numpy.array(niggli_cell)
            if len(vals) == 6:
                lattice_pars.append(vals)
                names.append(c.pdb_code)

        print "Total files loaded = {0}".format(i)

    database = [names, lattice_pars]

    cPickle.dump(database, open(pickle_file, "w"))

    return


if __name__ == "__main__":
    results = custom_report()
    create_niggli_cell_db(results)
