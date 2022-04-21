import numpy as np
import paths
import os


# Check
if os.getenv("CI", "false") == "true":
    raise Exception(
        "This script should never run on GitHub Actions, "
        "as its output should always be cached on Zenodo."
    )


# Generate samples
np.random.seed(0)
N = 3
mu = np.random.randn(N)
L = np.tril(np.random.randn(N))
cov = L @ L.T
samples = np.random.multivariate_normal(mu, cov, 100000)


# Store the results
np.savez(paths.data / "dataset.npz", samples=samples)
