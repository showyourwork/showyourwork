import numpy as np
from pathlib import Path
import os


# Check
if os.getenv("CI", "false") == "true":
    raise Exception(
        "This script should never run on GitHub Actions, "
        "as its output should always be cached on Zenodo."
    )


# Resolve paths
DATA = Path(__file__).resolve().parents[1] / "data"


# Generate 50 fake datasets and store them one per file
(DATA / "dataset").mkdir(exist_ok=True)
for n in range(50):
    np.random.seed(n)
    data = np.random.randn(100)
    np.savetxt(DATA / "dataset" / f"dataset{n:02d}.txt", data)
