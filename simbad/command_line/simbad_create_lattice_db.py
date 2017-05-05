"""Script to update the database for the lattice parameter search"""

__author__ = "Felix Simkovic & Adam Simpkin"
__date__ = "28 Apr 2017"
__version__ = "1.0"

import argparse
import cctbx.crystal
import numpy as np
import urllib2

import simbad.constants 
import simbad.command_line

logger = simbad.command_line.get_logger('create_lattice_db', level='info')


def rcsb_custom_report():
    """Create a custom report from the PDB

    Returns 
    -------
    list
        A list of results

    """
    base = 'http://www.rcsb.org/pdb/rest/'
    report = 'customReport.csv?pdbids=*&customReportColumns=lengthOfUnitCellLatticeA,'\
             + 'lengthOfUnitCellLatticeB,lengthOfUnitCellLatticeC,unitCellAngleAlpha,'\
             + 'unitCellAngleBeta,unitCellAngleGamma,spaceGroup,experimentalTechnique'\
             + '&service=wsfile&format=csv'
    
    # Space group conversions
    sg_conversion = {
        'A1': 'P1', 'B2': 'B112', 'C1211': 'C2', 'F422': 'I422',
        'I21': 'I2', 'I1211': 'I2', 'P21212A' : 'P212121',
        'R3': 'R3:R', 'C4212': 'P422',
    }

    # Obtain the data by querying the RCSB server 
    logger.info('Querying the RCSB Protein DataBank')
    results, errors = [], []
    for line in urllib2.urlopen(base + report):
        
        # Ignore the header
        if line.startswith('structureId'):
            continue
        
        pdb_code, rest = line.replace('"', "").split(',', 1)
        unit_cell, space_group, exp_tech = rest.rsplit(',', 2)
        
        # Ignore non-xtal structures
        if exp_tech.strip().upper() != "X-RAY DIFFRACTION":
            continue

        # Some entries do not have stored unit cell parameters
        try:
            unit_cell = map(float, unit_cell.split(','))
        except ValueError as e:
            logger.debug('Skipping pdb entry %s\t-\t%s', pdb_code, e)
            errors.append(pdb_code)
            continue

        space_group = space_group.replace(" ", "").strip()
        space_group = sg_conversion.get(space_group, space_group)   

        try:
            symmetry = cctbx.crystal.symmetry(unit_cell=unit_cell, space_group=space_group)
        except Exception as e:
            logger.debug('Skipping pdb entry %s\t-\t%s', pdb_code, e)
            errors.append(pdb_code)
            continue

        results.append((pdb_code, symmetry))
    
    logger.info('\tError with %d pdb entries', len(errors))
    return results


def create_niggli_cell_data(crystal_data):
    """Generate a new data for the lattice parameter search
    
    Parameters
    ----------
    crystal_data : list, tuple, generator
       A list of lists for PDB entries with the pdb code as 
       first and the crystal info as second argument.

    """
    logger.info('Calculating the Niggli cells')
    data_ascii_encoded = np.zeros((0, 10))
    for i, xtal_data in enumerate(crystal_data):
        d = np.zeros(10)
        d[:4] = np.fromstring(xtal_data[0], dtype='uint8')
        d[4:] = np.asarray(xtal_data[1].niggli_cell().unit_cell().parameters())
        data_ascii_encoded = np.concatenate((data_ascii_encoded, d[np.newaxis, :]), axis=0)
    logger.info("Total Niggli cells loaded: %d", i + 1)
    return data_ascii_encoded


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument('-db', '--database', help='File to store the database in. Default overwrites the general SIMBAD database')
    args = p.parse_args()
    
    if args.database:
        logger.info('Creating a new database file: %s', args.database)
        lattice_db = args.database
    else:
        logger.info('Overwriting the default SIMBAD Niggli database %s', simbad.constants.SIMBAD_LATTICE_DB)
        lattice_db = simbad.constants.SIMBAD_LATTICE_DB

    niggli_data = create_niggli_cell_data(rcsb_custom_report())
    np.savez_compressed(lattice_db, niggli_data)

