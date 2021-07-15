import numpy as np
import matplotlib.pyplot as plt

# The Heart Curve. See https://mathworld.wolfram.com/HeartCurve.html
fig = plt.figure(figsize=(7, 5))
t = np.linspace(-np.pi, np.pi, 1000)
x = 16 * np.sin(t) ** 3
y = 13 * np.cos(t) - 5 * np.cos(2 * t) - 2 * np.cos(3 * t) - np.cos(4 * t)
plt.plot(x, y, "C0-", lw=2)
plt.gca().axis("off")
fig.savefig("heart.pdf", bbox_inches="tight")
