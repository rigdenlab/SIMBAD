"""Main route into Simbad

The code here needs to be protected by if __name__ == '__main__'
to prevent multiprocessing on windows rerunning the main ample code when it
imports the main module again.

A bug in multiprocessing means that we can't use the __main__.py file and -m switch
to ccp4-python under windows.
"""

if __name__ == '__main__':
    import os, sys

    from simbad import main
    from simbad.util import exit_util

    try:
        main.SIMBAD().main()
    except Exception as e:
        msg = "Error running main SIMBAD program: {0}".format(e.message)
        exit_util.exit_error(msg, sys.exc_info()[2])

