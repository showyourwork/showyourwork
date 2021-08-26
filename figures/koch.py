"""Plot some pretty fractals.

This script plots the Koch snowflake (koch.pdf). The code below was adapted from

https://programmer.group/use-matplotlib-to-draw-snowflakes-and-snow-scenes.html

"""
import numpy as np
import matplotlib.pyplot as plt


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


# Panel 1
fig = plt.figure(figsize=(8, 8))
data = np.array(snow([(0, 0), (0.5, 0.8660254), (1, 0)], 3))
x, y = np.split(data, 2, axis=1)
plt.plot(x, y)
plt.gca().margins(0.2, 0.2)
plt.gca().axis("off")
fig.savefig("koch1.pdf", bbox_inches="tight")

# Panel 2
fig = plt.figure(figsize=(8, 8))
data = np.array(snow([(0, 0), (0.5, 0.8660254), (1, 0)], 5))
x, y = np.split(data, 2, axis=1)
plt.plot(x, y)
plt.gca().margins(0.2, 0.2)
plt.gca().axis("off")
fig.savefig("koch2.pdf", bbox_inches="tight")
