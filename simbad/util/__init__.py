"""Various miscellaneous functions"""

__author__ = "Adam Simpkin, Felix Simkovic & Jens Thomas"
__date__ = "05 May 2017"
__version__ = "1.0"

import glob
import logging
import math
import os
import pandas as pd
import shutil
import tempfile

from simbad.db import convert_pdb_to_dat
from simbad.util.pdb_util import PdbStructure

from pyjob.factory import TaskFactory

# Constants that need to be accessed externally (e.g. by CCP4I2)
SIMBAD_DIRNAME = "SIMBAD"
SIMBAD_PYRVAPI_SHAREDIR = "jsrview"
EXPORT = {"nt": "SET"}.get(os.name, "export")
CMD_PREFIX = {"nt": "call"}.get(os.name, "")

logger = logging.getLogger(__name__)


def get_sequence(input_f, output_s):
    """Output sequence file from input pdb file"""
    ps = PdbStructure.from_file(input_file=input_f)
    seq_info = ps.get_sequence_info

    content = "\n".join(">{}\n{}".format(i, seq_info[i]) for i in seq_info)
    with open(output_s, "w") as fh:
        fh.write(content)


def get_mrbump_ensemble(mrbump_dir, final):
    """Output ensemble from mrbump directory to a dat file"""
    if os.path.isdir(mrbump_dir):
        ensemble = glob.iglob(os.path.join(mrbump_dir, "models", "domain_*", "ensembles", "gesamtEnsTrunc_*_100.0_SideCbeta.pdb"))[0]
        convert_pdb_to_dat(ensemble, final)
    else:
        logger.critical("Directory missing: %s", mrbump_dir)


def output_files(run_dir, result, output_pdb, output_mtz):
    """Return output pdb/mtz from best result in result obj"""
    pdb_code = result[0]
    stem = os.path.join(run_dir, "output_files", pdb_code)
    input_pdb = os.path.join(stem, "{}_refinement_output.pdb".format(pdb_code))
    input_mtz = os.path.join(stem, "{}_refinement_output.mtz".format(pdb_code))
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
    df = pd.DataFrame([r._asdict() for r in results], **kwargs)
    df.set_index("pdb_code", inplace=True)

    if csv_file:
        logger.debug("Storing results in file: %s", csv_file)
        df.to_csv(csv_file)

    if df.empty:
        logger.info("No results found")
    else:
        logger.info("The results for this search are:\n\n%s\n", df.to_string())


def tmp_dir(directory=None, prefix="tmp", suffix=""):
    """Return a filename for a temporary directory

    Parameters
    ----------
    directory : str, optional
       Path to a directory to write the files to.
    prefix : str, optional
       A prefix to the temporary filename
    suffix : str, optional
       A suffix to the temporary filename

    """
    return tempfile.mkdtemp(dir=directory, prefix=prefix, suffix=suffix)


def tmp_file(delete=False, directory=None, prefix="tmp", stem=None, suffix=""):
    """Return a filename for a temporary file
    The naming convention of scripts will be ``prefix`` + ``stem`` + ``suffix``.

    Parameters
    ----------
    delete : bool, optional
       Delete the file, thus return name only [default: True]
    directory : str, optional
       Path to a directory to write the files to
    prefix : str, optional
       A prefix to the temporary filename
    stem : str, optional
       The steam part of the script name
    suffix : str, optional
       A suffix to the temporary filename
    """

    if directory is None:
        directory = tempfile.gettempdir()
    if stem is None:
        tmpf = tempfile.NamedTemporaryFile(delete=delete, dir=directory, prefix=prefix, suffix=suffix)
        tmpf.close()
        return tmpf.name
    else:
        tmpf = os.path.join(directory, "".join([prefix, stem, suffix]))
        if not delete:
            open(tmpf, "w").close()
        return tmpf


def source_ccp4():
    """Function to return bash command to source CCP4"""
    if os.name == "nt":
        return
    return "source {}".format(os.path.join(os.environ["CCP4"], "bin", "ccp4.setup-sh"))


def submit_chunk(collector, run_dir, nproc, job_name, submit_qtype, submit_queue, permit_nonzero, monitor, success_func):
    """Submit jobs in small chunks to avoid using too much disk space

    Parameters
    ----------
    collector : list
        :obj:`~pyjob.script.ScriptCollector` containing run scripts
    nproc : int, optional
        The number of processors to run the job on
    job_name : str
        The name of the job to submit
    submit_qtype : str
        The cluster submission queue type - currently support SGE and LSF
    submit_queue : str
        The queue to submit to on the cluster
    permit_nonzero : bool
        Permit non-zero return codes from TaskFactory
    success_func : func
        function to check for success

    """

    if submit_qtype == "local":
        processes = nproc
        array_size = None
    else:
        processes = None
        array_size = nproc

    with TaskFactory(
        submit_qtype,
        collector,
        cwd=run_dir,
        name=job_name,
        processes=processes,
        max_array_size=array_size,
        queue=submit_queue,
        permit_nonzero=permit_nonzero,
        shell="/bin/bash",
        priority=-10,
        cleanup=True,
    ) as task:
        task.run()
        interval = int(math.log(len(collector.scripts)) / 3)
        interval_in_seconds = interval if interval >= 5 else 5
        task.wait(interval=interval_in_seconds, monitor_f=monitor, success_f=success_func)
