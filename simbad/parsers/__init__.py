"""Parser classes for SIMBAD"""

__author__ = "Felix Simkovic"
__date__ = "04 May 2017"
__version__ = "0.1"


class _Parser(object):
    def __init__(self, fname):
        self.logfile = fname
