__author__ = "Adam Simpkin"
__date__ = "16 April 2018"
__version__ = "0.1"

import abc


class ScoreBase(object):
    """Abstract class for storing scores"""
    __metaclass__ = abc.ABCMeta

    @abc.abstractmethod
    def _as_dict(self):
        """Convert the object to a dictionary"""
        return
