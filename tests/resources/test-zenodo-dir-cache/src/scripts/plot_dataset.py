import numpy as np
import matplotlib.pyplot as plt
import corner
import paths


# Load the light curves
files = list((paths.data / "dataset").glob("dataset*.txt"))
nfiles = len(files)
data = np.empty((nfiles, 100))
for n, file in enumerate(files):
    data[n] = np.loadtxt(file)

# Plot
fig = plt.figure(figsize=(8, 5))
plt.plot(data)
fig.savefig(paths.figures / "dataset.pdf", bbox_inches="tight")