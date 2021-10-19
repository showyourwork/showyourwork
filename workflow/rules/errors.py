"""
Figure "error" (status) codes.

"""


class _ShowYourWorkError:
    """
    Showyourwork ``error'' codes for figure GitHub links.

    """

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


#: The figure script does not exist
ScriptDoesNotExist = _ScriptDoesNotExist()

#: The figure script is not version controlled
ScriptNotVersionControlled = _ScriptNotVersionControlled()

#: The figure script has uncommitted changes
ScriptHasUncommittedChanges = _ScriptHasUncommittedChanges()

#: The figure script is version controlled and up to date
ScriptUpToDate = _ScriptUpToDate()
