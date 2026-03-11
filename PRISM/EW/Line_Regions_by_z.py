import pandas as pd
import numpy as np

# Read the file (tsv = tab separated)
df = pd.read_csv("EW-Data-Output/EW_output_v6.tsv", sep="\t")

# Define columns
z_array = df.iloc[:, 1]

Ha_ini = df.iloc[:, 14]
Ha_end = df.iloc[:, 15]
Hb_ini = df.iloc[:, 16]
Hb_end = df.iloc[:, 17]


# Create masks for the different z ranges
groups = {
    "z < 3": z_array < 3,
    "3 <= z < 5": (z_array >= 3) & (z_array < 5),
    "5 <= z < 7": (z_array >= 5) & (z_array < 7),
    "z >= 7": z_array >= 7
}

results = []




# Compute averages
for label, mask in groups.items():
    avg_Ha_ini = Ha_ini[mask].mean()
    avg_Ha_end = Ha_end[mask].mean()
    avg_Hb_ini = Hb_ini[mask].mean()
    avg_Hb_end = Hb_end[mask].mean()

    results.append([label, avg_Ha_ini, avg_Ha_end, avg_Hb_ini, avg_Hb_end])


# Create output dataframe
output_df = pd.DataFrame(
    results,
    columns=["z_range", "Ha_ini", "Ha_end", "Hb_ini", "Hb_end"]
)

# Save to file
output_df.to_csv("Lines_Regions_by_z.tsv", sep="\t", index=False)


print("EW-Data-Output/Output written to EW_output_grouped.tsv")
