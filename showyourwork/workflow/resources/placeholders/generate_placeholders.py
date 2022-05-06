import matplotlib.pyplot as plt

fig, ax = plt.subplots(1, figsize=(6, 4))
ax.text(0.5, 0.5, "placeholder", fontsize=48, color="red", ha="center", va="center")
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.set_xticks([])
ax.set_yticks([])

figexts = plt.gcf().canvas.get_supported_filetypes().keys()
for ext in figexts:
    fig.savefig(f"placeholder.{ext}", bbox_inches="tight")