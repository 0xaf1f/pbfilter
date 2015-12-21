#!/usr/bin/env python
# EASY-INSTALL-ENTRY-SCRIPT: 'pbfilter==1.2','console_scripts','filter_artifacts.py'
__requires__ = 'pbfilter==1.2'
import sys
from pkg_resources import load_entry_point

sys.exit(
   load_entry_point('pbfilter==1.2', 'console_scripts', 'filter_artifacts.py')()
)
