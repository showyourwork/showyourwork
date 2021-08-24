from . import entry_points
from . import constants
from . import utils
from . import meta
from . import cache
from .showyourwork_version import __version__


# Path to the main Snakefile
Snakefile = str(
    constants.ROOT.parents[0].absolute()
    / "showyourwork-workflow"
    / "Snakefile"
)
