import numpy as np
import matplotlib.pyplot as plt
import corner
import paths


# Load the results
data = np.load(paths.data / "dataset.npz")
samples = data["samples"]


# Plot the corner plot
fig = corner.corner(samples)
fig.savefig(paths.figures / "dataset.pdf", bbox_inches="tight")