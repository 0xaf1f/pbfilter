import os
import sys
import logging
import functools
import time

from pbcore.util.Process import backticks

log = logging.getLogger(__name__)


def _validate_resource(func, resource):
    """Validate the existence of a file/dir"""
    if func(resource):
        return os.path.abspath(resource)
    else:
        raise IOError("Unable to find {f}".format(f=resource))

validate_file = functools.partial(_validate_resource, os.path.isfile)
validate_dir = functools.partial(_validate_resource, os.path.isdir)
validate_output_dir = functools.partial(_validate_resource, os.path.isdir)


def _nfs_exists_check(ff):
    """Return whether a file or a dir ff exists or not.
    Call ls instead of python os.path.exists to eliminate NFS errors.
    """
    # this is taken from Yuan
    cmd = "ls %s" % ff
    output, errCode, errMsg = backticks(cmd)
    return errCode == 0


def validate_fofn(fofn):
    """Validate existence of FOFN and files within the FOFN.

    :param fofn: (str) Path to File of file names.
    :raises: IOError if any file is not found.
    :return: (str) abspath of the input fofn
    """
    if os.path.isfile(fofn):
        file_names = fofn_to_files(os.path.abspath(fofn))
        log.debug("Found {n} files in FOFN {f}.".format(n=len(file_names), f=fofn))
        return os.path.abspath(fofn)
    else:
        raise IOError("Unable to find {f}".format(f=fofn))


def fofn_to_files(fofn):
    """Util func to convert a bas/bax fofn file to a list of bas/bax files."""
    if os.path.exists(fofn):
        with open(fofn, 'r') as f:
            bas_files = [line.strip() for line in f.readlines()]

        for bas_file in bas_files:
            if not os.path.isfile(bas_file):
                # try one more time to find the file by
                # performing an NFS refresh
                found = _nfs_exists_check(bas_file)
                if not found:
                    raise IOError("Unable to find bas/bax file '{f}'".format(f=bas_file))

        if len(set(bas_files)) != len(bas_files):
            raise IOError("Detected duplicate files in fofn. {n} unique, {m} total files in {f}".format(f=fofn, m=len(bas_files), n=len(set(bas_files))))

        return bas_files
    else:
        raise IOError("Unable to find FOFN {f}".format(f=fofn))


def setup_log(alog, level=logging.INFO, file_name=None, log_filter=None,
              str_formatter='[%(levelname)s] %(asctime)-15s [%(name)s %(funcName)s %(lineno)d] %(message)s'):
    """Core Util to setup log handler

    :param alog: a log instance
    :param level: (int) Level of logging debug
    :param file_name: (str, None) if None, stdout is used, str write to file
    :param log_filter: (LogFilter, None)
    :param str_formatter: (str) log formatting str
    """
    alog.setLevel(logging.DEBUG)

    if file_name is None:
        handler = logging.StreamHandler(sys.stdout)
    else:
        handler = logging.FileHandler(file_name)

    formatter = logging.Formatter(str_formatter)
    handler.setFormatter(formatter)
    handler.setLevel(level)

    if log_filter:
        handler.addFilter(log_filter)

    alog.addHandler(handler)


def log_timing(func):
    """Simple deco to log the runtime of func"""
    started_at = time.time()

    def wrapper(*args, **kw):
        return func(*args, **kw)

    run_time = time.time() - started_at
    name = func.__name__
    log.info("Func {f} took {s:.2f} sec ({m:.2f} min)".format(f=name, s=run_time, m=run_time / 60.0))

    return wrapper