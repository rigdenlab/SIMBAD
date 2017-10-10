"""Various miscellaneous functions"""

__author__ = "Adam Simpkin, Felix Simkovic & Jens Thomas"
__date__ = "05 May 2017"
__version__ = "1.0"

import logging
import pandas as pd

logger = logging.getLogger(__name__)


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
