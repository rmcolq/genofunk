from pkg_resources import get_distribution

try:
    __version__ = get_distribution("genofunk").version
except:
    __version__ = "local"

__all__ = ["annotate", "merge", "apply", "translate", "subcommands", "editfile", "refparser"]

from genofunk import *
