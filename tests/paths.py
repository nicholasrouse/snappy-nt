"""
This module provides functions for manipulating the file paths so that the pytest files
can find the relevant data in the package. Currently, it just works by manipulating
the __file__ variable. However, eventually, it should use importlib.resources. The issue
with using it now is that it requires Python 3.9. The goal is for the other testing
modules to be agnostic about the implementation, though we'll see how doable this is
when we update to importlib.resources.
"""

import os

import snappynt


def convert_rel_to_abs(path):
    """
    This converts a path relative to the snappynt directory to an absolute one usable
    by the filesystem.
    """
    return os.path.join(snappynt.__dict__["__path__"][0], path)
