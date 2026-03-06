import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Leer archivo (ajusta el nombre del archivo)
data = pd.read_csv("/home/nicolas/Documents/Research/PhD/JWST-Data/PRISM/EW/EW-Data-Output/EW_output_v5_2.tsv", sep='\t')

# Quitar " Angstrom" de todas las columnas que lo tengan
data = data.replace(" Angstrom", "", regex=True)


# Convertir todo lo posible a número
#data = data.apply(pd.to_numeric, errors="ignore")



# Guardar columnas en arrays
ID_array = data.iloc[:,0].to_numpy()
z_array = data.iloc[:,1].to_numpy()

#Halpha data
mast_Ha = data.iloc[:,2].to_numpy()
mast_Ha_err = data.iloc[:,3].to_numpy()
jades_Ha = data.iloc[:,4].to_numpy()
jades_Ha_err = data.iloc[:,5].to_numpy()


#Hbeta data
mast_Hb = data.iloc[:,6].to_numpy()
mast_Hb_err = data.iloc[:,7].to_numpy()
jades_Hb = data.iloc[:,8].to_numpy()
jades_Hb_err = data.iloc[:,9].to_numpy()


# Definir variables para la gráfica
mast_dHa = [float(mast_Ha_err[i]) for i in range(len(mast_Ha_err))] 
jades_dHa = [float(jades_Ha_err[i]) for i in range(len(jades_Ha_err))] 
mast_dHb = [float(mast_Hb_err[i]) for i in range(len(mast_Ha_err))] 
jades_dHb = [float(jades_Hb_err[i]) for i in range(len(jades_Ha_err))] 



#General Settings
plt.rc('axes', titlesize=18)     # fontsize of the axes title
plt.rc('axes', labelsize=18)     # fontsize of the x and y labels
plt.rc('axes', linewidth=1.65 )   # width of the frame
plt.rc('xtick', labelsize=16)    # fontsize of the tick labels
plt.rc('ytick', labelsize=16)    # fontsize of the tick labels
plt.rc('legend', fontsize=12)    # legend fontsize
plt.rc('font', size=18)          # controls default text sizes


# Scatter plot con barras de error
plt.figure(figsize=(8,8))

plt.errorbar(
    mast_Ha, jades_Ha,
    xerr=mast_dHa,
    yerr=jades_dHa,
    fmt='o',
    capsize=3,
    markersize=5
)
plt.title(r"MAST vs JADES ($H\alpha$)")
plt.plot(mast_Ha, mast_Ha, color='black', linestyle='--')
plt.xlabel("MAST")
plt.ylabel("JADES")

#plt.grid(True)
plt.tight_layout()
plt.savefig('EW_comparisons_Ha_2.pdf')
plt.show()
plt.close()


#///////// H beta //////////
plt.figure(figsize=(8,8))

plt.errorbar(
    mast_Hb, jades_Hb,
    xerr=mast_dHb,
    yerr=jades_dHb,
    fmt='o',
    capsize=3,
    markersize=5
)
plt.title(r"MAST  vs JADES ($H\beta$)")
plt.plot(mast_Hb, mast_Hb, color='black', linestyle='--')
plt.xlabel("MAST")
plt.ylabel("JADES")

#plt.grid(True)
plt.tight_layout()
plt.savefig('EW_comparisons_Hb_2.pdf')
plt.show()

