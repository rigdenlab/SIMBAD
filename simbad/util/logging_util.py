import logging
import sys


def setup_console_logging(level=logging.INFO,
                          formatstr='%(message)s\n' # Always add a blank line after every print
                          ):
    """
    Set up logging to the console for the root logger.

    :param level: int - Sets the threshold for the console output to level.
    :param formatstr: str - The string used to formate the log messages
    :return:
    """

    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)

    # First create console logger for outputting stuff
    # create file handler and set level to debug
    # Seems they changed the api in python 2.6->2.7
    try:
        cl = logging.StreamHandler(stream=sys.stdout)
    except TypeError:
        cl = logging.StreamHandler(stream=sys.stdout)
    cl.setLevel(level)
    formatter = logging.Formatter(formatstr)
    cl.setFormatter(formatter)
    logger.addHandler(cl)
    return logger

def setup_file_logging(logfile,
                       level=logging.DEBUG,
                       formatstr='%(asctime)s - %(name)s [%(lineno)d] - %(levelname)s - %(message)s'):
    """
    The path to the logfile that output will be written to

    :param logfile: The path to the logfile that output will be written to
    :param level: int - Sets the threshold for the console output to level.
    :param formatstr: str - The string used to format the log messages
    :return:
    """

    logger = logging.getLogger()
    fl = logging.FileHandler(logfile)
    fl.setLevel(level)
    formatter = logging.Formatter(formatstr)
    fl.setFormatter(formatter)
    logger.addHandler(fl)

    return logger
