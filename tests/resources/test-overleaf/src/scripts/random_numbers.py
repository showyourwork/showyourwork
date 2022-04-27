import matplotlib.pyplot as plt
import numpy as np
import paths

# Generate some data
random_numbers = np.random.randn(100, 10)

# Plot and save
fig = plt.figure(figsize=(7, 6))
plt.plot(random_numbers)
fig.savefig(paths.figures / "random_numbers.pdf", bbox_inches="tight", dpi=300)