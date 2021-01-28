from enum import Enum
import json
import logging.config
import os

from simbad import LOGGING_CONFIG


class LogColors(Enum):
    """Color container for log messages"""

    CRITICAL = 31
    DEBUG = 34
    DEFAULT = 0
    ERROR = 31
    WARNING = 33


class LogColorFormatter(logging.Formatter):
    """Formatter for log messages"""

    def format(self, record):
        if record.levelname in LogColors.__members__:
            prefix = '\033[1;{}m'.format(LogColors[record.levelname].value)
            postfix = '\033[{}m'.format(LogColors["DEFAULT"].value)
            record.msg = os.linesep.join([prefix + msg + postfix for msg in str(record.msg).splitlines()])
        return logging.Formatter.format(self, record)


def setup_logging(level, logfile=None, debugfile=None):
    """Read JSON config for logger and return root logger"""
    if not os.path.isfile(LOGGING_CONFIG):
        raise RuntimeError("Cannot find SIMBAD logging config file: {}".format(LOGGING_CONFIG))
    with open(LOGGING_CONFIG, 'rt') as f:
        config = json.load(f)

    # Reset some of the defaults
    config["handlers"]["console_handler"]["level"] = level.upper()
    if logfile:
        config["handlers"]["file_handler"]["filename"] = logfile
    if debugfile:
        config["handlers"]["debug_file_handler"]["filename"] = debugfile

    logging.config.dictConfig(config)
    return logging.getLogger()



