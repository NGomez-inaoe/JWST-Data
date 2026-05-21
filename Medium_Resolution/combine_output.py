import numpy as np
import pandas as pd
from pathlib import Path


# Cargar archivos

df1 = pd.read_csv("/home/nicolas/Documents/Research/PhD/JWST-Data/Medium_Resolution/Output_data/Line_fluxes_g395m_v4.tsv", sep='\t')
df2 = pd.read_csv("/home/nicolas/Documents/Research/PhD/JWST-Data/Medium_Resolution/Output_data/Line_fluxes_g235m_v4.tsv", sep='\t')
#df3 = pd.read_csv("/home/nicolas/Documents/Research/PhD/JWST-Data/Medium_Resolution/Output_data/Line_fluxes_g140m_v4.tsv", sep='\t')


# Columnas comunes que identifican cada fila
fix_columns = ["ID", "redshift"]
# Columnas que quieres combinar
combine_cols = ["Flux (Ha 6563) MAST",
        "F_err(Ha 6563) MAST",
        "Flux (Ha 6563) JADES",
        "F_err(Ha 6563) JADES",
        "Flux (Hb 4861) MAST",
        "F_err(Hb 4861) MAST",
        "Flux (Hb 4861) JADES",
        "F_err(Hb 4861) JADES",
        "Flux ([OIII]4959) MAST",
        "F_err([OIII]4959) MAST",
        "Flux ([OIII]4959) JADES",
        "F_err([OIII]4959) JADES",
        "Flux ([OIII]5007) MAST",
        "F_err([OIII]5007) MAST",
        "Flux ([OIII]5007) JADES",
        "F_err([OIII]5007) JADES",
        "Flux ([NII]6583) MAST",
        "F_err([NII]6583) MAST"]

# Unir los dos dataframes por las columnas clave
df_merge = df1.merge(
    df2,
    on=fix_columns,
    how="outer",
    suffixes=("_1", "_2")
)

# Función para combinar dos valores
def combinar_valores(a, b):
    if pd.notna(a) and pd.notna(b):
        return b
    elif pd.notna(a):
        return a
    elif pd.notna(b):
        return b
    else:
        return np.nan

# Crear dataframe final con las columnas clave
df_final = df_merge[fix_columns].copy()

# Combinar cada una de las columnas
for col in combine_cols:
    df_final[col] = df_merge.apply(
        lambda fila: combinar_valores(fila[f"{col}_1"], fila[f"{col}_2"]),
        axis=1
    )

# Guardar resultado
df_final.to_csv("Output_data/Line_fluxes_combined_v4.tsv", sep='\t', index=False)

print(df_final)