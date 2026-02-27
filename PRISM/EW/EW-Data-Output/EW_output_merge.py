import pandas as pd
import glob
import os


folder='./'
files=glob.glob(os.path.join(folder, "*.tsv"))

df_merged=pd.concat((pd.read_csv(f) for f in files), ignore_index=True)

df_merged.to_csv("EW_output_v4.tsv", sep='\t')

