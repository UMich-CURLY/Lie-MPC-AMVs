from domatplotlib import np, plt, loadmat

# Load data.
data = loadmat("vdposcillator.mat")
controllers = ["LMPC", "NMPC"]
colors = ["red", "green"]

# Make figure.
[fig, ax] = plt.subplots(figsize=(6, 3), nrows=3, sharex="col")
[x1ax, x2ax, uax] = ax

# Add plots.
for a in ax:
    a.axhline(color="gray", label="Setpoint")
for (controller, color) in zip(controllers, colors):
    x = data[controller]["x"]
    u = data[controller]["u"]
    u = np.concatenate((u, u[-1:]))
    t = data[controller]["t"]
    
    x1ax.plot(t, x[0,:], color=color, label=controller)
    x2ax.plot(t, x[1,:], color=color)
    uax.plot(t, u, color=color, drawstyle="steps-post")

# Clean up.
for [a, label] in zip(ax, ["$x_1$", "$x_2$", "$u$"]):
    a.margins(x=0, y=0.1)
    a.set_ylabel(label, rotation=0, labelpad=10)
uax.set_xlabel("Time")
leg = x1ax.legend(loc="lower center", ncol=3, bbox_to_anchor=(0.5, 1.01))
fig.tight_layout(pad=0.5, rect=(0, 0, 1, 0.9))

# Save figure.
fig.savefig("vdposcillator.pdf")
