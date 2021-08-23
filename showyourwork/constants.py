from pathlib import Path
import subprocess
import os


__all__ = [
    "TEMP",
    "ROOT",
    "USER",
    "FIGURE_EXTENSIONS",
    "UNKNOWN_SCRIPT",
    "ScriptDoesNotExist",
    "ScriptNotVersionControlled",
    "ScriptHasUncommittedChanges",
    "ScriptUpToDate",
]


# Useful paths
ROOT = Path(__file__).parents[0].absolute()
USER = Path(os.getcwd()).absolute()
TEMP = USER / ".showyourwork"

# Recognized figure extensions
FIGURE_EXTENSIONS = ["pdf", "png", "eps", "jpg", "jpeg", "gif", "svg", "tiff"]

#
UNKNOWN_SCRIPT = "unknown-script"


# Make temporary paths if needed
(TEMP / "tex").mkdir(exist_ok=True, parents=True)
(TEMP / "figures").mkdir(exist_ok=True, parents=True)
(TEMP / "data").mkdir(exist_ok=True, parents=True)
(TEMP / "tree").mkdir(exist_ok=True, parents=True)


# Error codes
class _ShowYourWorkError:
    def __str__(self):
        return "syw{}".format(self.__class__.__name__[1:])

    def __eq__(self, other):
        return self.code == other.code

    def __gt__(self, other):
        return self.code > other.code

    def __ge__(self, other):
        return self.code >= other.code

    def __lt__(self, other):
        return self.code < other.code

    def __le__(self, other):
        return self.code <= other.code


class _ScriptDoesNotExist(_ShowYourWorkError):
    code = 3


class _ScriptNotVersionControlled(_ShowYourWorkError):
    code = 2


class _ScriptHasUncommittedChanges(_ShowYourWorkError):
    code = 1


class _ScriptUpToDate(_ShowYourWorkError):
    code = 0


ScriptDoesNotExist = _ScriptDoesNotExist()
ScriptNotVersionControlled = _ScriptNotVersionControlled()
ScriptHasUncommittedChanges = _ScriptHasUncommittedChanges()
ScriptUpToDate = _ScriptUpToDate()
