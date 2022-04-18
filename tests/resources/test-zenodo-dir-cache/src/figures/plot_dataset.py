import numpy as np
import matplotlib.pyplot as plt
import corner
from pathlib import Path


# Resolve paths
DATA = Path(__file__).resolve().parents[1] / "data"
FIGURES = Path(__file__).resolve().parents[1] / "tex" / "figures"


# Load the light curves
files = list((DATA / "dataset").glob("dataset*.txt"))
nfiles = len(files)
data = np.empty((nfiles, 100))
for n, file in enumerate(files):
    data[n] = np.loadtxt(file)

# Plot
fig = plt.figure(figsize=(8, 5))
plt.plot(data)
fig.savefig(FIGURES / "dataset.pdf", bbox_inches="tight")