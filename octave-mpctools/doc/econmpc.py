from domatplotlib import np, plt, loadmat
from matplotlib import cm

colormap = getattr(cm, "viridis", cm.Set1)

def plot(trajectories, Delta=1):
    """Makes a phase plot and cost plot."""
    Ntraj = len(trajectories)
    colors = [colormap(c) for c in np.linspace(0, 1, Ntraj)]
    
    [fig, ax] = plt.subplots(figsize=(3, 3))
    for (traj, color) in zip(trajectories, colors):
        # Plot phase trajectory with collocation points.
        x = traj["x"]
        X = traj["X"]
        ax.plot(X[0,:], X[1,:], linestyle="-", marker="o", markersize=2,
                color=color, markerfacecolor="none", markeredgecolor=color)
        ax.plot(x[0,:], x[1,:], linestyle="", marker="o", color=color,
                markerfacecolor=color, markeredgecolor=color, markersize=4)
    
    # Clean up.
    ax.set_xlabel("$c_A$ (mol/L)")
    ax.set_ylabel("$c_B$ (mol/L)")
    ax.axis("equal")
    fig.tight_layout(pad=.5)
    
    # Return.
    return [fig, ax]

data = loadmat("econmpc.mat")
[fig, ax] = plot(data["Regularized"])
fig.savefig("econmpc.pdf")
