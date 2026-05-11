import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


# Latex rendering
plt.rc('font', **{'family': 'serif', 'serif': ['palatino']}) #Helvetica
plt.rc('text', usetex=True)
plt.rc('font', weight='bold')


sns.set_theme(
    style="ticks",
    context='paper',
    rc={
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.top": True,
        "ytick.right": True,
        "axes.grid": False,
    })


#General Settings
plt.rc('axes', titlesize=16)     # fontsize of the axes title
plt.rc('axes', labelsize=16)     # fontsize of the x and y labels
plt.rc('axes', linewidth=1.65 )   # width of the frame
plt.rc('xtick', labelsize=14)    # fontsize of the tick labels
plt.rc('ytick', labelsize=14)    # fontsize of the tick labels
plt.rc('legend', fontsize=12)    # legend fontsize
plt.rc('font', size=16)          # controls default text sizes

# Read the TSV file
df_g235m = pd.read_csv('./Output_data/Line_fluxes_g235m_v5_optimal.tsv', sep='\t')
df_g395m = pd.read_csv('./Output_data/Line_fluxes_g395m_v5_optimal.tsv', sep='\t')

# Create a sub dataframe with ID, redshift, and all columns containing "MAST"
mast_cols_g235m = [col for col in df_g235m.columns if 'MAST' in col]
mast_df_g235m = df_g235m[['ID', 'redshift'] + mast_cols_g235m]

mast_cols_g395m = [col for col in df_g395m.columns if 'MAST' in col]
mast_df_g395m = df_g395m[['ID', 'redshift'] + mast_cols_g395m]

# Create a sub dataframe with ID, redshift, and all columns containing "JADES"
#jades_cols = [col for col in df_g395m.columns if 'JADES' in col]
#jades_df = df_g235m[['ID', 'redshift'] + jades_cols]

# Extract the specific columns of interest
flux_o3_g235m = mast_df_g235m['Flux ([OIII]5007) MAST']
flux_hb_g235m = mast_df_g235m['Flux (Hb 4861) MAST']
flux_n2_g235m = mast_df_g235m['Flux ([NII]6583) MAST']
flux_ha_g235m = mast_df_g235m['Flux (Ha 6563) MAST']

flux_o3_g395m = mast_df_g395m['Flux ([OIII]5007) MAST']
flux_hb_g395m = mast_df_g395m['Flux (Hb 4861) MAST']
flux_n2_g395m = mast_df_g395m['Flux ([NII]6583) MAST']
flux_ha_g395m = mast_df_g395m['Flux (Ha 6563) MAST']

# Calculate log ratios (x and y for scatter plot)

y235 = np.log10(flux_o3_g235m / flux_hb_g235m)
x235 = np.log10(flux_n2_g235m / flux_ha_g235m)

y395 = np.log10(flux_o3_g395m / flux_hb_g395m)
x395 = np.log10(flux_n2_g395m / flux_ha_g395m)

# Kewley maximum starburst line-
x_kewly = np.linspace(-2.0, 0.4, 500)
y_kewly = 0.61 / (x_kewly - 0.47) + 1.19

# Plot
plt.figure(figsize=(12, 7.4165))
plt.scatter(x235, y235, alpha=0.6, color='purple',label='MAST G235M')
plt.scatter(x395, y395, alpha=0.6, color='grey', label='MAST G395M')
plt.plot(x_kewly, y_kewly, color='k', label='Kewley et al. 2001')
plt.minorticks_on()

plt.xlim(-1.8, 0.5)
plt.ylim(-1.3, 1.5)
plt.xlabel(r'$\log([NII]6583/H\alpha)$')
plt.ylabel(r'$\log([OIII]5007/H\beta)$') 

for x, y, obj_id in zip(x235, y235, mast_df_g235m['ID']):
    plt.annotate(
        str(obj_id),
        (x, y),
        textcoords="offset points",
        xytext=(3, 3),
        ha="left",
        fontsize=8,
        color='purple'
    )

for x, y, obj_id in zip(x395, y395, mast_df_g395m['ID']):
    plt.annotate(
        str(obj_id),
        (x, y),
        textcoords="offset points",
        xytext=(3, 3),
        ha="left",
        fontsize=8,
        color='grey'
    )

plt.legend(loc='lower left')
plt.savefig('./Output_data/BPT_diagram.pdf', bbox_inches='tight')
plt.show()
