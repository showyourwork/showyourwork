"""Plot a pretty fractal.

This script plots the Mandelbrot set (mandelbrot.pdf).
The code below was adapted from

https://scipy-lectures.org/intro/numpy/auto_examples/plot_mandelbrot.html

"""
import numpy as np
import matplotlib.pyplot as plt
from numpy import newaxis
import copy


def compute_mandelbrot(N_max, some_threshold, nx, ny):
    x = np.linspace(-2, 1, nx)
    y = np.linspace(-1.5, 1.5, ny)
    c = x[:, newaxis] + 1j * y[newaxis, :]
    z = c
    with np.warnings.catch_warnings():
        np.warnings.simplefilter("ignore")
        for j in range(N_max):
            z = z ** 2 + c
        mandelbrot_set = abs(z) < some_threshold
    return mandelbrot_set


fig = plt.figure(figsize=(8, 8))
mandelbrot_set = np.round(1 - compute_mandelbrot(500, 50.0, 601, 401))
tab10 = copy.copy(plt.get_cmap("tab10"))
tab10.set_over("w")
plt.imshow(
    mandelbrot_set.T,
    extent=[-2, 1, -1.5, 1.5],
    interpolation="nearest",
    cmap=tab10,
    vmin=0,
    vmax=0.9,
)
plt.gca().axis("off")
fig.savefig("mandelbrot.pdf", bbox_inches="tight")
