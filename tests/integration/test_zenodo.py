import pytest
from helpers import (
    ShowyourworkRepositoryActions,
    TemporaryShowyourworkRepository,
)

from showyourwork.config import edit_yaml

pytestmark = pytest.mark.remote

test_figure_script = """
import matplotlib.pyplot as plt
import numpy as np
import json
import paths

# Set up the figure
fig = plt.figure(figsize=(7, 6))
ax = fig.subplot_mosaic(
    '''
    AB
    CC
    '''
)
fig.subplots_adjust(wspace=0.1, hspace=0.15)

# Images
for n, dataset in zip(["A", "B"], ["S06", "S07"]):
    with open(paths.data / "TOI640" / f"{dataset}.json", "r") as f:
        image_data = json.load(f)
    aperture = np.array(image_data["aperture"])
    image = np.array(image_data["image"]) / 1e3
    im = ax[n].imshow(image, origin="lower", vmin=0, vmax=4)
    ax[n].set_xticks([])
    ax[n].set_yticks([])
    ax[n].set_title(dataset)
    ax[n].set_aspect("auto")
fig.colorbar(im, ax=[ax["A"], ax["B"]], shrink=1, label="flux [ke$^-$/s]")

# Dummy colorbar to get the spacing right
cbar = fig.colorbar(im, ax=ax["C"], shrink=1)
cbar.ax.set_visible(False)

# Load planet data
with open(paths.data / "TOI640" / "planet.json", "r") as f:
    planet_data = json.load(f)
phase = planet_data["phase"]
period = planet_data["period"]

# Photometry
for dataset in ["S06", "S07"]:

    # Load the data
    time, flux, *_ = np.loadtxt(paths.data / "TOI640" / f"{dataset}.txt").T

    # Divide out the rolling 1-day median baseline
    dt = np.median(np.diff(time))
    w = int(1 / dt)
    baseline = [np.median(flux[a : a + w]) for a in range(len(flux))]
    flux /= baseline

    # Fold the data
    time = (time % period) / period + phase

    # Plot it
    ax["C"].plot(time, flux, ".", label=dataset)

# Format the figure
ax["C"].legend()
ax["C"].set_xlim(-0.1, 0.1)
ax["C"].set_xlabel("orbital phase")
ax["C"].set_ylabel("normalized flux")
fig.savefig(paths.figures / "test_figure.pdf", bbox_inches="tight", dpi=300)
"""


class TestZenodo(
    TemporaryShowyourworkRepository, ShowyourworkRepositoryActions
):
    """Test a repo that downloads data from Zenodo."""

    def customize(self):
        """Add files and edit settings to customize the repo."""
        # Add the figure script
        with open(self.cwd / "src" / "scripts" / "test_figure.py", "w") as f:
            print(test_figure_script, file=f)

        # Add the figure environment to the tex file
        self.add_figure_environment()

        # Edit the config file to specify deps & zenodo metadata
        with edit_yaml(self.cwd / "showyourwork.yml") as config:
            config["dependencies"] = {
                "src/scripts/test_figure.py": [
                    "src/data/TOI640/planet.json",
                    "src/data/TOI640/S06.txt",
                    "src/data/TOI640/S06.json",
                    "src/data/TOI640/S07.txt",
                    "src/data/TOI640/S07.json",
                ]
            }
            config["datasets"] = {
                "10.5281/zenodo.6515864": {
                    "destination": "src/data/TOI640",
                    "contents": {
                        "README.md": None,
                        "TOI640b.json": "src/data/TOI640/planet.json",
                        "images.tar.gz": {
                            "README.md": None,
                            "S06": {"image.json": "src/data/TOI640/S06.json"},
                            "S07": {"image.json": "src/data/TOI640/S07.json"},
                        },
                        "lightcurves.zip": {
                            "lightcurves": {
                                "README.md": None,
                                "S06": {"lc.txt": "src/data/TOI640/S06.txt"},
                                "S07": {"lc.txt": "src/data/TOI640/S07.txt"},
                            }
                        },
                    },
                }
            }
