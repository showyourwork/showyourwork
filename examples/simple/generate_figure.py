import sys

import matplotlib.pyplot as plt
import numpy as np

x = np.linspace(0, 4 * np.pi, 400)
plt.plot(x, np.sin(x))
plt.savefig(sys.argv[1], bbox_inches="tight")
