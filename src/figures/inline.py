import matplotlib.pyplot as plt

fig, ax = plt.subplots(1)
ax.add_patch(plt.Polygon([[0, 0], [0.5, 0], [0.25, 0.5]], color="C0"))
ax.add_patch(plt.Polygon([[0.5, 0], [1, 0], [0.75, 0.5]], color="C0"))
ax.add_patch(plt.Polygon([[0.25, 0.5], [0.75, 0.5], [0.5, 1]], color="C0"))
ax.axis("off")
fig.savefig("inline.pdf", bbox_inches="tight")
