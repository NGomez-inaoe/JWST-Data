import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Read the TSV file
df = pd.read_csv('./Output_data/Line_fluxes_combined_v4.tsv', sep='\t')

# Create a sub dataframe with ID, redshift, and all columns containing "MAST"
mast_cols = [col for col in df.columns if 'MAST' in col]
sub_df = df[['ID', 'redshift'] + mast_cols]

# Extract the specific columns of interest
flux_oiii = sub_df['Flux ([OIII]5007) MAST']
flux_hb = sub_df['Flux (Hb 4861) MAST']
flux_n2 = sub_df['Flux ([NII]6583) MAST']
flux_ha = sub_df['Flux (Ha 6563) MAST']

# Calculate log ratios (x and y for scatter plot)
# y: log10([OIII]/Hb)
# x: log10([NII]/Ha)
y = np.log10(flux_oiii / flux_hb)
x = np.log10(flux_n2 / flux_ha)

# Kewley et al. (2001) maximum starburst line-
x_kewly = np.linspace(-2.0, 0.4, 500)
y_kewly = 0.61 / (x_kewly - 0.47) + 1.19

# Plot
plt.figure(figsize=(7, 6))
plt.plot(x_kewly, y_kewly, color='k', label='Kewley+2001 (AGN boundary)')
plt.xlim(-2.0, 0.5)
plt.ylim(-1.5, 1.5)
# Create scatter plot
plt.scatter(x, y, alpha=0.6)

plt.xlabel(r'$\log([NII]6583/H\alpha)$')
plt.ylabel(r'$\log([OIII]5007/H\beta)$') 

plt.show()
