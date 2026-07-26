import pandas as pd
df_properties = pd.read_csv("/home/nicolas/Documents/Research/PhD/JWST-Data/Medium_Resolution/Physical_Properties/physical_properties.tsv", sep='\t')

Te_JADES = df_properties['Te JADES']
Te_MAST = df_properties['Te MAST']


import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

#%matplotlib widget
#================= Plot Settings =================#
# Latex rendering
plt.rc('font', **{'family': 'serif', 'serif': ['palatino']}) #Helvetica
plt.rc('text', usetex=True)
plt.rc('font', weight='bold')

# --- Cosmology palette (normalized RGB) ---


COSMO_COLORS = {
        "color2": (73/255, 153/255, 51/255),
        "color1": (67/255, 55/255, 63/255),
        "deep_blue": (30/255, 60/255, 120/255),
        "warm_gold": (200/255, 140/255, 0/255),
        "dark_teal": (0/255, 130/255, 90/255),
        "soft_crimson": (150/255, 40/255, 40/255),
        "cool_gray": (120/255, 120/255, 120/255),
        }


sns.set_theme(
    style="ticks",
    context='talk',
    rc={
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.top": True,
        "ytick.right": True,
        "axes.grid": True,
    })

sns.set_palette([
        COSMO_COLORS['color1'],
        COSMO_COLORS['color2'],
        COSMO_COLORS["deep_blue"],
        COSMO_COLORS["warm_gold"],
        COSMO_COLORS["cool_gray"],
        COSMO_COLORS["dark_teal"],
        COSMO_COLORS["soft_crimson"],
        "black",
        ])

#General Settings
plt.rc('axes', titlesize=16)     # fontsize of the axes title
plt.rc('axes', labelsize=16)     # fontsize of the x and y labels
plt.rc('axes', linewidth=1.65 )   # width of the frame
plt.rc('xtick', labelsize=14)    # fontsize of the tick labels
plt.rc('ytick', labelsize=14)    # fontsize of the tick labels
plt.rc('legend', fontsize=12)    # legend fontsize
plt.rc('font', size=16)          # controls default text sizes




gray_red_cmap = LinearSegmentedColormap.from_list(
    'gray_red',
    [COSMO_COLORS['color2'], COSMO_COLORS['warm_gold'],  'red']
)

def parse_asymmetric_error(err):
    if isinstance(err, str):
        try:
            lo, hi = err.strip("()").split(",")
            return float(lo), float(hi)
        except ValueError:
            return np.nan, np.nan
    return np.nan, np.nan

errors_jades = df_properties["Te err JADES"].map(parse_asymmetric_error)
errors_mast = df_properties["Te err MAST"].map(parse_asymmetric_error)

df_plot = df_properties.copy()
df_plot["Te_JADES_err_lo"] = errors_jades.map(lambda x: x[0])
df_plot["Te_JADES_err_hi"] = errors_jades.map(lambda x: x[1])
df_plot["Te_MAST_err_lo"] = errors_mast.map(lambda x: x[0])
df_plot["Te_MAST_err_hi"] = errors_mast.map(lambda x: x[1])

mask = df_plot["Te JADES"].notna() & df_plot["Te MAST"].notna()
df_plot = df_plot[mask]

x = df_plot["Te JADES"]
y = df_plot["Te MAST"]
xerr = np.vstack([
    df_plot["Te_JADES_err_lo"].fillna(0.0),
    df_plot["Te_JADES_err_hi"].fillna(0.0),
])
yerr = np.vstack([
    df_plot["Te_MAST_err_lo"].fillna(0.0),
    df_plot["Te_MAST_err_hi"].fillna(0.0),
])

fig, ax = plt.subplots(figsize=(8, 6))
scatter = ax.scatter(
    x,
    y,
    c=df_plot["redshift"],
    cmap=gray_red_cmap,
    s=80,
    edgecolor="k",
    linewidth=0.7,
)

ax.errorbar(
    x,
    y,
    xerr=xerr,
    yerr=yerr,
    fmt="none",
    ecolor="gray",
    alpha=0.7,
    capsize=3,
    zorder=0,
)

for _, row in df_plot.iterrows():
    ax.text(
        row["Te JADES"],
        row["Te MAST"],
        str(int(row["ID"])),
        fontsize=8,
        ha="left",
        va="bottom",
    )

cbar = fig.colorbar(scatter, ax=ax)
cbar.set_label("redshift")

ax.set_xlabel("Te JADES [K]")
ax.set_ylabel("Te MAST [K]")
ax.set_title("Te JADES vs Te MAST with asymmetric uncertainties")
ax.grid(True, linestyle=":", alpha=0.6)
fig.tight_layout()
lims = [min(x.min(), y.min()), max(x.max(), y.max())]
ax.plot(lims, lims, color="gray", linestyle="--", linewidth=1, zorder=-1)
ax.set_title("")
ax.grid(False)
plt.savefig('temperatures_mast_jades.pdf')
plt.show()
