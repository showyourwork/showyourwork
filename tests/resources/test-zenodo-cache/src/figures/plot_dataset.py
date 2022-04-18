import numpy as np
import matplotlib.pyplot as plt
import corner
from pathlib import Path


# Resolve paths
DATA = Path(__file__).resolve().parents[1] / "data"
FIGURES = Path(__file__).resolve().parents[1] / "tex" / "figures"


# Load the results
data = np.load(DATA / "dataset.npz")
samples = data["samples"]


# Plot the corner plot
fig = corner.corner(samples)
fig.savefig(FIGURES / "dataset.pdf", bbox_inches="tight")