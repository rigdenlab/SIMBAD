"""Script to download or update SIMBAD-related databases"""

__author__ = "Felix Simkovic & Adam Simpkin"
__date__ = "17 May 2017"
__version__ = "0.2"

import argparse
import numpy as np
import sys
import time
import urllib2

import cctbx.crystal
import simbad.constants 
import simbad.command_line
import simbad.exit

logger = None

# The space groups in the list below cannot be recognized by CCTBX, so we convert them
# to similar ones understandle by the library
CCTBX_ERROR_SG = {
    'A1': 'P1', 'B2': 'B112', 'C1211': 'C2', 'F422': 'I422', 'I21': 'I2', 'I1211': 'I2', 
    'P21212A' : 'P212121', 'R3': 'R3:R', 'C4212': 'P422',
}

def create_lattice_db(database):
    """Create a lattice search database

    Parameters
    ----------
    database : str
       The path to the database file

    """
    logger.info('Querying the RCSB Protein DataBank')

    url = 'http://www.rcsb.org/pdb/rest/customReport.csv?pdbids=*&customReportColumns=lengthOfUnitCellLatticeA,'\
          + 'lengthOfUnitCellLatticeB,lengthOfUnitCellLatticeC,unitCellAngleAlpha,unitCellAngleBeta,unitCellAngleGamma,'\
          + 'spaceGroup,experimentalTechnique&service=wsfile&format=csv'

    crystal_data, errors = [], []
    for line in urllib2.urlopen(url):
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
            logger.debug('Skipping pdb entry %s\t%s', pdb_code, e)
            errors.append(pdb_code)
            continue
        space_group = space_group.replace(" ", "").strip()
        space_group = CCTBX_ERROR_SG.get(space_group, space_group)   
        try:
            symmetry = cctbx.crystal.symmetry(unit_cell=unit_cell, space_group=space_group)
        except Exception as e:
            logger.debug('Skipping pdb entry %s\t%s', pdb_code, e)
            errors.append(pdb_code)
            continue
        crystal_data.append((pdb_code, symmetry))
    logger.debug('Error processing %d pdb entries', len(errors))

    logger.info('Calculating the Niggli cells')
    niggli_data = np.zeros((len(crystal_data), 10))
    for i, xtal_data in enumerate(crystal_data):
        niggli_data[i][:4] = np.fromstring(xtal_data[0], dtype='uint8').astype(np.float64)
        niggli_data[i][4:] = np.asarray(xtal_data[1].niggli_cell().unit_cell().parameters())
    logger.info("Total Niggli cells loaded: %d", len(crystal_data))

    if not database.endswith('.npz'):
        database += ".npz"
    logger.info('Storing database in file: %s', database)
    np.savez_compressed(database, niggli_data)


def create_db_argparse():
    """Argparse function database creationg"""
    p = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    sp = p.add_subparsers(help='Database-specific commands')

    pa = sp.add_parser('lattice', help='lattice database')
    pa.set_defaults(which="lattice")
    simbad.command_line._argparse_lattice_options(pa)

    return p


def main():
    """SIMBAD database creation function"""
    args = create_db_argparse().parse_args()
    # Logger setup
    global logger
    args.debug_lvl = 'debug'
    logger = simbad.command_line.setup_logging(level=args.debug_lvl)
    # Print a fancy header
    simbad.command_line.print_header()
    # Take a time snapshot
    time_start = time.time()
    # Create the requested database
    if args.which == "lattice":
        create_lattice_db(args.latt_db)
    # Calculate and display the runtime 
    days, hours, mins, secs = simbad.command_line.calculate_runtime(time_start, time.time())
    logger.info("Database creation completed in %d days, %d hours, %d minutes, and %d seconds", days, hours, mins, secs)


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        msg = "Error running main SIMBAD program: {0}".format(e.message)
        simbad.exit.exit_error(*sys.exc_info())
