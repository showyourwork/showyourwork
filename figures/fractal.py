import numpy as np
import matplotlib.pyplot as plt
from numpy import newaxis


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


def rotate(p, d):
    a = np.radians(d)
    m = np.array([[np.cos(a), np.sin(a)], [-np.sin(a), np.cos(a)]])
    return np.dot(p, m)


def koch_curve(p, q):
    p, q = np.array(p), np.array(q)
    u = p + (q - p) / 3
    v = q - (q - p) / 3
    w = rotate(v - u, 60) + u
    return u.tolist(), v.tolist(), w.tolist()


def snow(triangle, k):
    for i in range(k):
        result = list()
        t_len = len(triangle)
        for j in range(t_len):
            p = triangle[j]
            q = triangle[(j + 1) % t_len]
            u, v, w = koch_curve(p, q)
            result.extend([p, u, w, v])
        triangle = result.copy()
    triangle.append(triangle[0])
    return triangle


# Panel A: The Mandelbrot set
# https://scipy-lectures.org/intro/numpy/auto_examples/plot_mandelbrot.html
fig = plt.figure(figsize=(8, 8))
mandelbrot_set = np.round(1 - compute_mandelbrot(500, 50.0, 601, 401))
tab10 = plt.get_cmap("tab10").copy()
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
fig.savefig("fractal_mandelbrot.pdf", bbox_inches="tight")

# Panel B: A Koch snowflake
# https://programmer.group/use-matplotlib-to-draw-snowflakes-and-snow-scenes.html
fig = plt.figure(figsize=(8, 8))
data = np.array(snow([(0, 0), (0.5, 0.8660254), (1, 0)], 5))
x, y = np.split(data, 2, axis=1)
plt.plot(x, y)
plt.gca().margins(0.2, 0.2)
plt.gca().axis("off")
fig.savefig("fractal_koch.pdf", bbox_inches="tight")
