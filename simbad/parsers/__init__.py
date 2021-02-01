"""Parser classes for SIMBAD"""

__author__ = "Felix Simkovic & Adam Simpkin & Filomeno Sanchez Rodriguez"
__date__ = "27 May 2020"
__version__ = "0.1"


import abc
import os
import logging

ABC = abc.ABCMeta('ABC', (object,), {})


class _Parser(ABC):
    """Parser abstract class
    This class contains general methods and data structures to extract information from the output created by the
    classes at :py:obj`~simbad.parsers`

    Attributes
    ----------
    fname : str
        Name of the file to be parsed (default: None)
    stdout : str
        The stdout to be parsed (default: None)
    logger : :py:obj:`logger`
        Loggig interface for the parser (default: None)

    Raises
    ------
    BoolError
        If an error has occured while parsing the file of stdout
    """

    def __init__(self, fname=None, stdout=None, logger=None):
        self.fname = fname
        self.error = False
        self.stdout = stdout
        self.inputfile_contents = None
        if logger is None:
            self.logger = logging.getLogger(__name__)
        else:
            self.logger = logger

    @abc.abstractmethod
    def parse(self):
        """ Abstract method to run the parser"""
        pass
