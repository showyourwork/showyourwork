from pathlib import Path
import subprocess
import os


__all__ = [
    "TEMP",
    "ROOT",
    "USER",
    "PROJECT",
    "FIGURE_EXTENSIONS",
    "ScriptDoesNotExist",
    "ScriptNotVersionControlled",
    "ScriptHasUncommittedChanges",
    "ScriptUpToDate",
]


def get_project_name():
    """
    Attempt to infer the name of the current git repo.

    """
    try:
        url = (
            subprocess.check_output(
                ["git", "config", "--get", "remote.origin.url"],
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
            )
            .decode()
            .replace("\n", "")
        )
        name = url.split("/")[-1]
    except:
        try:
            name = (
                Path(
                    subprocess.check_output(
                        ["git", "rev-parse", "--show-toplevel"],
                        stdout=subprocess.DEVNULL,
                        stderr=subprocess.DEVNULL,
                    )
                    .decode()
                    .replace("\n", "")
                )
                .absolute()
                .name
            )
        except:
            name = Path(os.getcwd()).absolute().name
    return name


# Useful paths
TEMP = Path.home() / ".showyourwork"
ROOT = Path(__file__).parents[0].absolute()
USER = Path(os.getcwd()).absolute()


# Name of the user's project
PROJECT = get_project_name()


# Recognized figure extensions
FIGURE_EXTENSIONS = ["pdf", "png", "eps", "jpg", "jpeg", "gif", "svg", "tiff"]


# Make temporary paths if needed
(TEMP / PROJECT / "tex").mkdir(exist_ok=True, parents=True)
(TEMP / PROJECT / "figures").mkdir(exist_ok=True, parents=True)
(TEMP / PROJECT / "tree").mkdir(exist_ok=True, parents=True)
(TEMP / "user").mkdir(exist_ok=True, parents=True)


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
