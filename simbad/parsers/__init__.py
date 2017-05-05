"""Parser classes for SIMBAD"""

__author__ = "Felix Simkovic"
__date__ = "04 May 2017"
__version__ = "0.1"


class _Parser(object):
    """Daddy Parser Yo
    
    Do no instantiate directly

    """
    def __init__(self, f):
        """Create a new parser"""
        self._logfile = None

        self.logfile = f

    @property
    def logfile(self):
        """Return logfile"""
        return self._logfile

    @logfile.setter
    def logfile(self, f):
        """Define logfile"""
        self._logfile = f


