"""Module to run the rotation search"""

__author__ = "Adam Simpkin & Felix Simkovic"
__date__ = "12 April 2018"
__version__ = "0.3"

import math
import pyjob


def rotation_search_factory(method):
    if method == "amore":
        from simbad.rotsearch.amore_search import AmoreRotationSearch
        return AmoreRotationSearch
    elif method == "phaser":
        from simbad.rotsearch.phaser_search import PhaserRotationSearch
        return PhaserRotationSearch
    else:
        raise ValueError("Unrecognised program entered to perform the rotation search: %s", method)


def mr_job_succeeded(r_fact, r_free):
    """Check values for job success"""
    return r_fact < 0.45 and r_free < 0.45


def get_chunk_size(total, size):
    return total if size == 0 else size


def get_total_chunk_cycles(total, step):
    total_chunk_cycles, remainder = divmod(total, step)
    if remainder > 0:
        return total_chunk_cycles + 1
    else:
        return total_chunk_cycles


def submit_chunk(chunk_scripts, run_dir, nproc, job_name, submit_qtype, submit_queue, monitor, success_func):
    """Submit jobs in small chunks to avoid using too much disk space

    Parameters
    ----------
    chunk_scripts : list
        List of scripts for each chunk
    nproc : int, optional
        The number of processors to run the job on
    job_name : str
        The name of the job to submit
    submit_qtype : str
        The cluster submission queue type - currently support SGE and LSF
    submit_queue : str
        The queue to submit to on the cluster
    success_func : func
        function to check for success

    """
    j = pyjob.Job(submit_qtype)
    j.submit(chunk_scripts, directory=run_dir, name=job_name, nproc=nproc,
             max_array_jobs=nproc, queue=submit_queue, permit_nonzero=True, priority=-10)
    interval = int(math.log(len(chunk_scripts)) / 3)
    interval_in_seconds = interval if interval >= 5 else 5
    j.wait(interval=interval_in_seconds, monitor=monitor, check_success=success_func)

