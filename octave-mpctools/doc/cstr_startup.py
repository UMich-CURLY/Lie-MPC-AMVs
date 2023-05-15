# Linear and nonlinear control of startup of a CSTR.
from domatplotlib import np, plt, loadmat
import matplotlib.gridspec as gridspec

# Load data.
data = loadmat("cstr_startup.mat")
x = data["xcl"]
u = data["ucl"]
xsp = data["xsp"]
contvars = data["contvars"] - 1 # Need zero-based indices.

# Define plotting function.
def cstrplot(x, u, xsp=None, contvars=[], title=None, colors={}, labels={},
             markers={}, keys=None, bounds=None, ilegend=0, Delta=1):
    if keys is None:
        keys = list(x.keys())
    for k in keys:    
        u[k] = np.concatenate((u[k], np.NaN*np.ones_like(u[k][:,-1:])), axis=1)
    xsp = np.concatenate((xsp, np.NaN*np.ones_like(xsp[:,-1:])), axis=1)
    k = keys[0]
    Nx = x[k].shape[0]
    Nu = u[k].shape[0]
    ylabelsx = ["$c$ (mol/L)", "$T$ (K)", "$h$ (m)"]
    ylabelsu = ["$T_c$ (K)", "$F$ (kL/min)"]
    
    gs = gridspec.GridSpec(Nx*Nu,2)    
    fig = plt.figure(figsize=(6, 4))
    leglines = []
    leglabels = []
    for i in range(Nx):
        ax = fig.add_subplot(gs[i*Nu:(i+1)*Nu,0])
        for k in keys:
            t = np.arange(0, x[k].shape[1])*Delta
            args = {"color": colors.get(k, "black"), "label": labels.get(k, k),
                    "marker": markers.get(k, "")}
            [line] = ax.plot(t, x[k][i,:], markeredgecolor="none", **args)
            if i == ilegend:
                leglines.append(line)
                leglabels.append(args["label"])
        if i in contvars and xsp is not None:
            ax.step(t, xsp[i,:], linestyle="--", color="black", where="post")
        ax.set_ylabel(ylabelsx[i])
        ax.margins(x=0, y=0.1)
    ax.set_xlabel("Time (min)")
    for i in range(Nu):
        ax = fig.add_subplot(gs[i*Nx:(i+1)*Nx,1])
        for k in keys:
            t = np.arange(0, u[k].shape[1])*Delta
            args = {"color":colors.get(k, "black"), "label": labels.get(k, k)}
            ax.step(t, u[k][i,:], where="post", **args)
        if bounds is not None:
            for b in set(["uub", "ulb"]).intersection(bounds.keys()):
                ax.plot(np.array([t[0], t[-1]]), np.ones((2,))*bounds[b][i],
                        '--k')
        ax.set_ylabel(ylabelsu[i])
        ax.margins(x=0, y=0.1)
    ax.set_xlabel("Time (min)")
    fig.legend(leglines,leglabels,loc="lower center",ncol=len(keys))
    fig.tight_layout(pad=.5, rect=(0,.075,1,1))
    if title is not None:
        fig.canvas.set_window_title(title)
    return fig

# Make plots.
keys = ["uncont", "lmpc", "nmpc"]
colors = {"lmpc": "blue", "nmpc": "green", "uncont": "red"}
labels = {"lmpc": "LMPC", "nmpc": "NMPC", "uncont": "Uncontrolled"}
markers = {"lmpc": "s", "nmpc": "o", "uncont": "^"}
plotbounds = {k : data[k] for k in ["ulb","uub"]}
fig = cstrplot(x, u, xsp, colors=colors, contvars=contvars, labels=labels,
               keys=keys, markers={}, bounds=plotbounds, ilegend=2)
fig.savefig("cstr_startup.pdf")
