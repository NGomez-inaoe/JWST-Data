from pathlib import Path
from astropy.io import fits
from astropy.table import Table
import shutil
import re
import numpy as np
import matplotlib.pyplot as plt


source=Path("/home/nicolas/Documents/Research/PhD/JWST-Data/JWST Data/PRISM/Spectra1D/")

mast_folder=source / "MAST"
jades_folder=source / "JADES"



#Extract IDs and plugs them into an array
IDs_array = []


for file in jades_folder.iterdir():
    match=re.search(r"deephst-100(\d+)", file.name)
    if match:
        obj_ID = match.group(1)
        IDs_array.append(obj_ID)     
        
#Definte the function that receives the files and the ID and make the plot
def plot_spectrum(MAST_file, JADES_file, ID):
  
    
    # ~ fits.open('file.fits')
    t_MAST=Table.read(MAST_file, hdu="EXTRACT1D")
    w_MAST=t_MAST['WAVELENGTH']
    flux_MAST=t_MAST['FLUX'] #flux in Jy
    #z=1
    #rest_w=w/(1+z) #in micras
    #rest_w_Angstrom=rest_w*1e-6/1e-10 # rest wavelength in Angstrom
    w_Angstrom=w_MAST*1e-6/1e-10 # observed wavelength in Angstrom
    flux_ergs = 2.99792458E-05 *flux_MAST / w_Angstrom**2. #flux in erg s^-1 cm^-2 A^-1



    #archivo de JADES
    t_JADES=Table.read(JADES_file, hdu="EXTRACT1D")
    w_JADES=t_JADES['WAVELENGTH']
    flux_JADES=t_JADES['FLUX'] #flux in Jy

    
    plt.rcParams["figure.figsize"] = (8, 6)
    plt.title(f'ID {ID}, z=', fontsize=18.)
    
    plt.step(w_MAST/1e-4,flux_ergs, label='MAST Data-Archive', color='Blue')
    plt.step(w_JADES/1e-4,flux_JADES, label='JADES Survey', color='Green')
    plt.xlabel(r'$\lambda_{obs}$ [$\AA$]', fontsize=22)
    plt.ylabel(r'FLUX [erg s$^{-1}$ cm$^{-2}$ $\AA^{-1}$]', fontsize=22)
    legend=plt.legend(loc=4,labelspacing=0.1)
    legend=plt.legend(loc='best',labelspacing=0.1)
    plt.setp(legend.get_texts(),fontsize='14')
    plt.savefig(f"PRISM_Plots/{ID}.pdf")
    plt.close()
    
for ID in IDs_array:
    mast_file = './Spectra1D/MAST/jw01210-o001_s0000'+ID+'_nirspec_clear-prism_x1d.fits'
    jades_file = './Spectra1D/JADES/hlsp_jades_jwst_nirspec_goods-s-deephst-100'+ID+'_clear-prism_v1.0_x1d.fits'

    plot_spectrum(mast_file, jades_file, ID)