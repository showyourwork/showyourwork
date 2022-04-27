import numpy as np
import paths
import os


# Check
if os.getenv("CI", "false") == "true":
    raise Exception(
        "This script should never run on GitHub Actions, "
        "as its output should always be cached on Zenodo."
    )


# Generate 50 fake datasets and store them one per file
(paths.data / "dataset").mkdir(exist_ok=True)
for n in range(50):
    np.random.seed(n)
    data = np.random.randn(100)
    np.savetxt(paths.data / "dataset" / f"dataset{n:02d}.txt", data)
