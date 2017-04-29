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

import simbad.constants 


class _LatticeParameters(object):
    """A basic lattice parameter sorting class"""

    __slots__ = ('pdb_code', 'a', 'b', 'c', 'alpha', 'beta', 'gamma', 'space_group')

    def __init__(self, pdb_code, a, b, c, alpha, beta, gamma, space_group):
        self.pdb_code = pdb_code
        self.a = a
        self.b = b
        self.c = c
        self.alpha = alpha
        self.beta = beta
        self.gamma = gamma
        self.space_group = space_group

    def __repr__(self):
        return "{0}(pdb_code={1} a={2} b={3} c={4} alpha={5}, beta={6}, gamma={7}, space_group={8}".format(
            self.__class__.__name__, self.pdb_code, self.a, self.b, self.c,
            self.alpha, self.beta, self.gamma, self.space_group)

    def _as_dict(self):
        """Convert the :obj:`_LatticeParameterScore <simbad.lattice.search._LatticeParameterScore>`
        object to a dictionary"""
        dictionary = {}
        for k in self.__slots__:
            dictionary[k] = getattr(self, k)
        return dictionary


def custom_report():
    """Create a custom report from the PDB and returns a :obj: _LatticeParameters containing 
    the pdb code, lattice parameters and space group for every structure in the PDBs"""
    
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
                    a = float(tmp[1])
                    b = float(tmp[2])
                    c = float(tmp[3])
                    alpha = float(tmp[4])
                    beta = float(tmp[5])
                    gamma = float(tmp[6])
                    space_group = tmp[7][:-2]

                    data = _LatticeParameters(pdb_code, a, b, c, alpha, beta, gamma, space_group)
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
    """Generate a new database for the lattice parameter search, overwrites old file"""
    
#     pickle_file = 'niggli_database.cpk'
    pickle_file = simbad.constants.SIMBAD_LATTICE_DB
    names = []
    lattice_pars = []

    for i, c in enumerate(crystal_data):
        unit_cell = ' '.join(str(p) for p in [c.a, c.b, c.c, c.alpha, c.beta, c.gamma])

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
