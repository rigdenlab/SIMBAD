"""Script to update the database for the lattice parameter search"""

__author__ = "Felix Simkovic & Adam Simpkin"
__date__ = "28 Apr 2017"
__version__ = "0.1"

import argparse
import numpy as np
import os
import sys
import time
import urllib2

import cctbx.crystal
import simbad.constants 
import simbad.command_line
import simbad.util.exit_util

logger = None


def create_db_argparse():
    """Argparse function database creationg"""
    p = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    simbad.command_line._argparse_lattice_options(p)
    # In simbad.command_line._argparse_core_options but don't want other options here to create for now
    p.add_argument('-debug_lvl', type=str, default='info',
                    help='The console verbosity level < notset | info | debug | warning | error | critical > ')
    return p.parse_args()


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
            logger.debug('Skipping pdb entry %s\t%s', pdb_code, e)
            errors.append(pdb_code)
            continue

        space_group = space_group.replace(" ", "").strip()
        space_group = sg_conversion.get(space_group, space_group)   

        try:
            symmetry = cctbx.crystal.symmetry(unit_cell=unit_cell, space_group=space_group)
        except Exception as e:
            logger.debug('Skipping pdb entry %s\t%s', pdb_code, e)
            errors.append(pdb_code)
            continue

        results.append((pdb_code, symmetry))
    
    logger.debug('Error processing %d pdb entries', len(errors))
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
    data_ascii_encoded = np.zeros((len(crystal_data), 10))
    for i, xtal_data in enumerate(crystal_data):
        data_ascii_encoded[i][:4] = np.fromstring(xtal_data[0], dtype='uint8').astype(np.float64)
        data_ascii_encoded[i][4:] = np.asarray(xtal_data[1].niggli_cell().unit_cell().parameters())
    logger.info("Total Niggli cells loaded: %d", len(crystal_data))
    return data_ascii_encoded


def main():
    """SIMBAD database creation function"""
    args = create_db_argparse()

    # Logger setup
    global logger
    logger = simbad.command_line.setup_logging(level=args.debug_lvl)

    # Print some fancy header
    simbad.command_line.print_header()
    
    # Take a time snapshot
    time_start = time.time()
    
    # Generate the Niggli data
    niggli_data = create_niggli_cell_data(rcsb_custom_report())

    # Create the library
    lattice_db = args.latt_db
    if not lattice_db.endswith('.npz'):
        lattice_db += ".npz"
    logger.info('Storing database in file: %s', lattice_db)
    np.savez_compressed(lattice_db, niggli_data)
    
    # Calculate and display the runtime 
    days, hours, mins, secs = simbad.command_line.calculate_runtime(time_start, time.time())
    logger.info("Database creation completed in %d days, %d hours, %d minutes, and %d seconds", days, hours, mins, secs)


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        msg = "Error running main SIMBAD program: {0}".format(e.message)
        simbad.util.exit_util.exit_error(msg, sys.exc_info()[2])
