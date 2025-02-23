"""
This type stub file was generated by pyright.
"""

import os
import sys
from .version import version as __version__
from . import config as _config
from .utils.state import ScienceState
from .tests.runner import TestRunner
from .logger import _init_log, _teardown_log
from .utils.misc import find_api_page
from types import ModuleType as __module_type__

"""
Astropy is a package intended to contain core functionality and some
common tools needed for performing astronomy and astrophysics research with
Python. It also provides an index for other astronomy packages and tools for
managing them.
"""
if 'dev' in __version__:
    online_docs_root = ...
else:
    online_docs_root = ...
class Conf(_config.ConfigNamespace):
    """
    Configuration parameters for `astropy`.
    """
    unicode_output = ...
    use_color = ...
    max_lines = ...
    max_width = ...


conf = ...
class base_constants_version(ScienceState):
    """
    Base class for the real version-setters below
    """
    _value = ...
    _versions = ...
    @classmethod
    def validate(cls, value): # -> str:
        ...
    
    @classmethod
    def set(cls, value): # -> _ScienceStateContext:
        """
        Set the current constants value.
        """
        ...
    


class physical_constants(base_constants_version):
    """
    The version of physical constants to use
    """
    _value = ...
    _versions = ...


class astronomical_constants(base_constants_version):
    """
    The version of astronomical constants to use
    """
    _value = ...
    _versions = ...


test = ...
__citation__ = ...
log = ...
def online_help(query): # -> None:
    """
    Search the online Astropy documentation for the given query.
    Opens the results in the default web browser.  Requires an active
    Internet connection.

    Parameters
    ----------
    query : str
        The search query.
    """
    ...

__dir_inc__ = ...
