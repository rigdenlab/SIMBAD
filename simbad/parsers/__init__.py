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
        self.check_input()

    # ------------------ Abstract methods and properties ------------------

    @abc.abstractmethod
    def parse(self):
        """ Abstract method to run the parser"""
        pass

    @property
    @abc.abstractmethod
    def summary(self):
        """Abstract property to store a summary of the parsed figures of merit"""
        pass

    # ------------------ Some general methods ------------------

    def check_input(self):
        """Check if :py:attr:`~swamp.parsers.parser.fname` exists"""

        if self.fname is not None and not os.path.isfile(self.fname):
            self.error = True
            self.logger.error('Cannot find input file, please make sure it exists: %s' % self.fname)