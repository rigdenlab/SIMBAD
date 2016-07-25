'''
Various miscellaneous functions
Largely stolen from ample
'''

import logging
import os
import shlex
import subprocess
import tempfile

_logger = logging.getLogger(__name__)

def ccp4_version():
    """Return the CCP4 version as a tuple"""
    global CCP4_VERSION
    if CCP4_VERSION is None:
        # Currently there seems no sensible way of doing this other then running a program and grepping the output
        cmd=['pdbcur']
        logf = tempfile.NamedTemporaryFile(delete=False)
        run_command(cmd, stdin="", logfile=logf.name)
        logf.seek(0) # rewind logfile
        tversion=None
        for i, line in enumerate(logf):
            if i > 20:break
            if line.startswith(' ### CCP4'):
                tversion=line.split()[2].rstrip(':')
                break
        
        logf.close()
        if not tversion: raise RuntimeError,"Cannot determine CCP4 version"
        vsplit = tversion.split('.')
        if len(vsplit) == 2:
            major = int(vsplit[0])
            minor =  int(vsplit[1])
            rev = '-1'
        elif len(vsplit) == 3:
            major = int(vsplit[0])
            minor = int(vsplit[1])
            rev = int(vsplit[2])
        else: raise RuntimeError,"Cannot split CCP4 version: {0}".format(tversion)
    os.unlink(logf.name)
    return (major,minor,rev)

def filename_append(filename=None, astr=None,directory=None, separator="_"):
    """Append astr to filename, before the suffix, and return the new filename."""
    dirname, fname = os.path.split( filename )
    name, suffix = os.path.splitext( fname )
    name  =  name + separator + astr + suffix
    if directory is None:
        directory = dirname
    return os.path.join( directory, name )

def is_exe(fpath):
    return fpath and os.path.exists(fpath) and os.access(fpath, os.X_OK)

def run_command(cmd, logfile=None, directory=None, dolog=True, stdin=None, check=False, env=None):
    """Execute a command and return the exit code.

    We take care of outputting stuff to the logs and opening/closing logfiles

    Args:
    cmd - command to run as a list
    stdin - a string to use as stdin for the command
    logfile (optional) - the path to the logfile
    directory (optional) - the directory to run the job in (cwd assumed)
    dolog: bool - whether to output info to the system log
    """
    assert type(cmd) is list, "run_command needs a list!"
    if check:
        if not is_exe(cmd[0]): raise RuntimeError,"run_command cannot find executable: {0}".format(cmd[0])

    if not directory:  directory = os.getcwd()
    if dolog: _logger.debug("In directory {0}\nRunning command: {1}".format(directory, " ".join(cmd)))
    file_handle=False
    if logfile:
        if type(logfile)==file:
            file_handle=True
            logf=logfile
            logfile=os.path.abspath(logf.name)
        else:
            logfile = os.path.abspath(logfile)
            logf = open(logfile, "w")
        if dolog: _logger.debug("Logfile is: {0}".format(logfile))
    else:
        logf = tempfile.TemporaryFile()
        
    if stdin != None:
        stdinstr = stdin
        stdin = subprocess.PIPE

    # Windows needs some special treatment
    kwargs = {}
    if os.name == "nt":
        kwargs = { 'bufsize': 0, 'shell' : "False" }
    p = subprocess.Popen(cmd, stdin=stdin, stdout=logf, stderr=subprocess.STDOUT, cwd=directory, env=env, **kwargs)

    if stdin != None:
        p.stdin.write( stdinstr )
        p.stdin.close()
        if dolog: _logger.debug("stdin for cmd was: {0}".format( stdinstr ) )

    p.wait()
    if not file_handle: logf.close()
    return p.returncode

def run_job(command_line, logfile, key=""):
    """ Generic job runner """

    if os.name == "nt":
        process_args = shlex.split(command_line, posix=False)
        p = subprocess.Popen(process_args, bufsize=0, shell="False", stdin = subprocess.PIPE,
                             stdout = subprocess.PIPE, stderr = subprocess.STDOUT)
    else:
        process_args = shlex.split(command_line)
        p = subprocess.Popen(process_args, stdin = subprocess.PIPE,
                             stdout = subprocess.PIPE, stderr = subprocess.STDOUT)


    (child_stdin, child_stdout_and_stderr) = (p.stdin, p.stdout)

    # Write the keyword input
    child_stdin.write(key)
    child_stdin.close()

    # Watch the output for successful termination
    out=child_stdout_and_stderr.readline()

    log=open(logfile, "w")

    while out:
        out=child_stdout_and_stderr.readline()
        log.write(out)

    child_stdout_and_stderr.close()
    log.close()
