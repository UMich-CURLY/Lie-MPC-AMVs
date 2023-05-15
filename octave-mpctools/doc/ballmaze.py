from domatplotlib import loadmat, plt

# Load data.
data = loadmat("ballmaze.mat")
xmax = data["xmax"]
cushion = 0.05*xmax
circles = data["circles"]
holes = [tuple(circles[i,:]) for i in range(circles.shape[0])]

# Make plot.
[fig, ax] = plt.subplots(figsize=(3, 3))

for (x, y, r) in holes:
    circ = plt.Circle((x, y), r, edgecolor="red", facecolor=(1,0,0,.5))
    ax.add_artist(circ)

x = data["x"]
ax.plot(x[0,:], x[1,:], color="black", marker="o", markersize=1.5,
        linewidth=0.25)
ax.set_xlabel("$x_1$")
ax.set_ylabel("$x_2$", rotation=0, labelpad=10)
ax.set_xlim((-cushion, xmax + cushion))
ax.set_ylim((-cushion, xmax + cushion))
    
fig.tight_layout(pad=0.5)
fig.savefig("ballmaze.pdf")
