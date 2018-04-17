"""Various miscellaneous functions"""

__author__ = "Adam Simpkin, Felix Simkovic & Jens Thomas"
__date__ = "05 May 2017"
__version__ = "1.0"

import glob
import logging
import os
import pandas as pd
import shutil

# Constants that need to be accessed externally (e.g. by CCP4I2)
SIMBAD_DIRNAME = 'SIMBAD'
SIMBAD_PYRVAPI_SHAREDIR = 'jsrview'
EXPORT = "SET" if os.name == "nt" else "export"
CMD_PREFIX = "call" if os.name == "nt" else ""

logger = logging.getLogger(__name__)


def output_files(run_dir, result, output_pdb, output_mtz):
    """Return output pdb/mtz from best result in result obj"""
    pdb_code = result[0]
    stem = os.path.join(run_dir, 'mr_search', pdb_code, 'mr', '*', 'refine')
    input_pdb = glob.glob(os.path.join(stem, '{0}_refinement_output.pdb'.format(pdb_code)))[0]
    input_mtz = glob.glob(os.path.join(stem, '{0}_refinement_output.mtz'.format(pdb_code)))[0]
    shutil.copyfile(input_pdb, output_pdb)
    shutil.copyfile(input_mtz, output_mtz)


def result_by_score_from_csv(f, score, ascending=True):
    """Return result with the best defined score"""
    df = pd.read_csv(f)
    df.sort_values(score, ascending=ascending, inplace=True)
    return df.loc[0, ["pdb_code", score]].tolist()


def summarize_result(results, csv_file=None, columns=None):
    """Summarize the search results"""
    kwargs = {}
    if columns:
        kwargs["columns"] = ["pdb_code"] + columns
    df = pd.DataFrame([r._as_dict() for r in results], **kwargs)
    df.set_index("pdb_code", inplace=True)

    if csv_file:
        logger.debug("Storing results in file: %s", csv_file)
        df.to_csv(csv_file)

    if df.empty:
        logger.info("No results found")
    else:
        summary_table = "The results for this search are:\n\n%s\n"
        logger.info(summary_table, df.to_string())
