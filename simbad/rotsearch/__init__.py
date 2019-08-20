"""Module to run the rotation search"""

__author__ = "Adam Simpkin & Felix Simkovic"
__date__ = "12 April 2018"
__version__ = "0.3"


def rotation_search_factory(method):
    if method == "amore":
        from simbad.rotsearch.amore_search import AmoreRotationSearch

        return AmoreRotationSearch
    elif method == "phaser":
        from simbad.rotsearch.phaser_search import PhaserRotationSearch

        return PhaserRotationSearch
    else:
        raise ValueError("Unrecognised program entered to perform the rotation search: %s", method)


def get_chunk_size(total, size):
    return total if size == 0 else size


def get_total_chunk_cycles(total, step):
    total_chunk_cycles, remainder = divmod(total, step)
    if remainder > 0:
        return total_chunk_cycles + 1
    else:
        return total_chunk_cycles
