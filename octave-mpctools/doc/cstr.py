import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
import plottools
plottools._init(family="sans-serif", sansmath=True)

def cstrplot(data, title=None, Delta=1):
    # Read data.
    x = data["x"]
    u = data["u"]
    ysp = data.get("ysp") # This one is optional.
    
    # Get sizes.
    (Nx, Nt) = x.shape
    (Nu, _) = u.shape
    
    # Grab vectors.
    u = dupfinalcol(u)
    t = np.arange(0, Nt)*Delta
    ylabelsx = ["$c$ (mol/L)", "$T$ (K)", "$h$ (m)"]
    ylabelsu = ["$T_c$ (K)", "$F$ (kL/min)"]
    if ysp is None:
        ysp = np.NaN*np.ones((Nx, Nt))
    else:
        ysp = dupfinalcol(ysp)
    
    # Make plots.
    gs = gridspec.GridSpec(Nx*Nu, 2)    
    fig = plt.figure(figsize=(6, 3.5))
    for i in range(Nx):
        ax = fig.add_subplot(gs[i*Nu:(i + 1)*Nu,0])
        ax.plot(t, x[i,:], "-ok", markersize=2.5)
        ax.step(t, ysp[i,:], "-r", where="post")
        ax.set_ylabel(ylabelsx[i])
        ax.margins(x=0, y=0.1)
    ax.set_xlabel("Time (min)")
    for i in range(Nu):
        ax = fig.add_subplot(gs[i*Nx:(i + 1)*Nx,1])
        ax.step(t, u[i,:], "-k", where="post")
        ax.set_ylabel(ylabelsu[i])
        ax.margins(x=0, y=0.1)
    ax.set_xlabel("Time (min)")
    fig.tight_layout(pad=.5)
    if title is not None:
        fig.canvas.set_window_title(title)
    return fig

def dupfinalcol(arr):
    """Duplicates the final column of an array for easy stairstep plotting."""
    return np.concatenate((arr, arr[:,-1:]), axis=1)

# Load data and plot.
data = plottools.loadmat("cstr.mat")
fig = cstrplot(data)
fig.savefig("cstr.pdf")
