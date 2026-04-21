import numpy as np
import pandas as pd
from pathlib import Path


# Cargar archivos

df2 = pd.read_csv("/home/nicolas/Documents/Research/PhD/JWST-Data/Medium_Resolution/Output_data/Line_fluxes_g395m_v1.tsv", sep='\t')
df1 = pd.read_csv("/home/nicolas/Documents/Research/PhD/JWST-Data/Medium_Resolution/Output_data/Line_fluxes_g235m_v1.tsv", sep='\t')
df1 = pd.read_csv("/home/nicolas/Documents/Research/PhD/JWST-Data/Medium_Resolution/Output_data/Line_fluxes_g140m_v1.tsv", sep='\t')


# Columnas comunes que identifican cada fila
fix_columns = ["ID", "redshift"]
# Columnas que quieres combinar
combine_cols = ["LF(Ha) MAST", "LF err(Ha) MAST", "LF(Ha) JADES", "LF err(Ha) JADES", "LF(Hb) MAST","LF err(Hb) MAST",
                "LF(Hb) JADES","LF err(Hb) JADES","LF([OIII]4959) MAST","LF err([OIII]4959) MAST","LF([OIII]4959) JADES",
                "LF err([OIII]4959) JADES","LF([OIII]5007) MAST","LF err([OIII]5007) MAST","LF([OIII]5007) JADES","LF err([OIII]5007) JADES"]

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
        return (a + b) / 2
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
df_final.to_csv("Output_data/Line_fluxes_combined.tsv", sep='\t', index=False)

print(df_final)