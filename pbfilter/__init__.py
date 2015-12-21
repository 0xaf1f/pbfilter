VERSION = (1, 2)

# $Date: 2014/08/24 $
# $Revision: #1 $

_changelist = "$Change: 139483 $"


def get_version():
    """Return the version as a string. "O.7"

    This uses a major.minor

    Each python module of the system (e.g, butler, detective, siv_butler.py)
    will use this version +  individual changelist. This allows top level
    versioning, and sub-component to be versioned based on a p4 changelist.

    .. note:: This should be improved to be compliant with PEP 386.
    """
    return ".".join([str(i) for i in VERSION])


def _get_changelist(perforce_str):
    import re
    rx = re.compile(r'Change: (\d+)')
    match = rx.search(perforce_str)
    if match is None:
        v = 'UnknownChangelist'
    else:
        try:
            v = int(match.group(1))
        except (TypeError, IndexError):
            v = "UnknownChangelist"
    return v


def get_changelist():
    return _get_changelist(_changelist)
