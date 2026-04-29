import warnings
import numpy as np
import pandas as pd
from pathlib import Path
from astropy.io import fits
from astropy import units as u
import matplotlib.pyplot as plt
from specutils.analysis import line_flux
from astropy.nddata import StdDevUncertainty
from specutils import Spectrum, SpectralRegion


#Define source path where the files are
source = Path('/home/nicolas/Documents/Research/PhD/JWST-Data/Medium_Resolution/Spectra1D')
jades_folder = source / 'JADES'
mast_folder = source / 'MAST'


#Extract data with Pandas
df = pd.read_csv("list_of_candidates.tsv", sep='\t')
ID_array = df["NIRSpec_ID"]
z_array = df["redshift"]
JADES_files = df["JADES_FILENAME_F170LP-G235M"]
MAST_files = df["MAST_FILENAME_F170LP-G235M"]
Ha_ini_array = df["Ha_ini"]
Ha_end_array = df["Ha_end"]
Hb_ini_array = df["Hb_ini"]
Hb_end_array = df["Hb_end"]
O3_ini_array = df["[OIII4959]_ini"]
O3_end_array = df["[OIII4959]_end"]
o3_ini_array = df["[OIII5007]_ini"]
o3_end_array = df["[OIII5007]_end"]

#Create arrays to save the data with save_EW()
ID_data =[]
z_data = []
mast_Ha_data = []
mast_Ha_err_data =[]
mast_Hb_data = []
mast_Hb_err_data =[]
mast_O3_data = []
mast_O3_err_data =[]
mast_o3_data = []
mast_o3_err_data =[]

jades_Ha_data = []
jades_Ha_err_data =[]
jades_Hb_data = []
jades_Hb_err_data= []
jades_O3_data =[]
jades_O3_err_data=[]
jades_o3_data =[]
jades_o3_err_data=[]



def main():
    
    
    for index in range(62):
        save_line_fluxes(index, plot=False)


    #Save data in data frame
    LF_data = {
        "ID": ID_data,
        "redshift": z_data,
        "LF(Ha) MAST": mast_Ha_data,
        "LF err(Ha) MAST": mast_Ha_err_data,
        "LF(Ha) JADES": jades_Ha_data,
        "LF err(Ha) JADES": jades_Ha_err_data,
        "LF(Hb) MAST": mast_Hb_data,
        "LF err(Hb) MAST": mast_Hb_err_data,
        "LF(Hb) JADES": jades_Hb_data,
        "LF err(Hb) JADES": jades_Hb_err_data,
        "LF([OIII]4959) MAST": mast_O3_data,
        "LF err([OIII]4959) MAST": mast_O3_err_data,
        "LF([OIII]4959) JADES": jades_O3_data,
        "LF err([OIII]4959) JADES": jades_O3_err_data,
        "LF([OIII]5007) MAST": mast_o3_data,
        "LF err([OIII]5007) MAST": mast_o3_err_data,
        "LF([OIII]5007) JADES": jades_o3_data,
        "LF err([OIII]5007) JADES": jades_o3_err_data
    }
    
    LineFlux_df=pd.DataFrame( LF_data)
    LineFlux_df.to_csv('./Output_data/Line_fluxes_g235m_v3.tsv', sep='\t', index=False)
    
        


#==================================================
#/// /////         Functions                 //////
#==================================================
def flux_stdDev(lambda_array, flux_array, regionLeft, regionRight):

    #Compute the standard deviation of the flux around a given line
    lambda_array = np.array(lambda_array)
    flux_array = np.array(flux_array)

    # Crear máscara booleana para el rango deseado
    lambda_ini = regionLeft[0]
    lambda_end = regionLeft[1]
    mask_left = (lambda_array >= lambda_ini) & (lambda_array <= lambda_end)

    lambda_ini = regionRight[0]
    lambda_end = regionRight[1]
    mask_right = (lambda_array >= lambda_ini) & (lambda_array <= lambda_end)

    # Extraer los valores correspondientes
    flux_left = flux_array[mask_left]
    flux_right = flux_array[mask_right]

    stdDev_left = np.std(flux_left)
    stdDev_right = np.std(flux_right)
    if np.isnan(stdDev_left):
        stdDev_left = 0
    if np.isnan(stdDev_right):
        stdDev_right = 0

    return (stdDev_left + stdDev_right)/2
#------------------------------------------------
#
#------------------------------------------------
def compute_line_flux(lamb, flux, region, z):
    
    r_ini = region[0] 
    r_end = region[1] 
    vlight = 2.99792458e-5 
    
    with warnings.catch_warnings(): #ignore warnings
        warnings.simplefilter('ignore')


        if np.isfinite(r_ini) and np.isfinite(r_end):

            #Flux Uncerainty
            if r_end > 5050 or r_ini < 4900:
                flux_err = flux_stdDev(lamb.value, flux.value, [r_ini - 50, r_ini], [r_end, r_end + 50]) 
            else:
                flux_err = flux_stdDev(lamb.value, flux.value, [4950 - 50, 4950], [5015, 5015 + 50]) 

            flux_err_Jy = flux_err * np.ones(len(flux)) * u.Jy
            flux_uncertainty = StdDevUncertainty(flux_err_Jy)
            flux_err_Jy = flux_err * np.ones(len(flux)) * u.Jy
            flux_uncertainty = StdDevUncertainty(flux_err_Jy)

            spectrum = Spectrum(spectral_axis=lamb, flux=flux, uncertainty=flux_uncertainty)

    
            try:
                l_flux = line_flux(spectrum, regions=SpectralRegion(r_ini* u.AA, r_end* u.AA))
                l_flux_unc = l_flux.uncertainty
                #Convert to ergs s-1 cm-2
                l_flux = l_flux.value * vlight / ((1+float(z))*np.mean([r_ini, r_end]))**2
                l_flux_unc = l_flux_unc.value * vlight / ((1+float(z))*np.mean([r_ini, r_end]))**2

        
            except IndexError:
                l_flux = np.nan
                l_flux_unc = np.nan

        else:
            l_flux = np.nan
            l_flux_unc = np.nan


    return l_flux*1e18, l_flux_unc*1e18
#------------------------------------------------
#
#------------------------------------------------        
def get_line_fluxes(lambda_rest, flux, indx):

    #Arrays for the lines
    z = z_array[indx]
    Ha_region = [Ha_ini_array[indx], Ha_end_array[indx] ] 
    Hb_region = [Hb_ini_array[indx], Hb_end_array[indx] ] 
    O3_region = [O3_ini_array[indx], O3_end_array[indx] ]
    o3_region = [o3_ini_array[indx], o3_end_array[indx] ]
    
    flux_Jy = flux * u.Jy
    lambda_rest_AA = lambda_rest *u.AA
    
    
    Ha_flux, Ha_flux_unc = compute_line_flux(lambda_rest_AA, flux_Jy, Ha_region, z) 
    Hb_flux, Hb_flux_unc = compute_line_flux(lambda_rest_AA, flux_Jy, Hb_region, z)
    O3_flux, O3_flux_unc = compute_line_flux(lambda_rest_AA, flux_Jy, O3_region, z)
    o3_flux, o3_flux_unc = compute_line_flux(lambda_rest_AA, flux_Jy, o3_region, z)


    return Ha_flux, Ha_flux_unc, Hb_flux, Hb_flux_unc, O3_flux, O3_flux_unc, o3_flux, o3_flux_unc
#------------------------------------------------
#
#------------------------------------------------            
def save_line_fluxes(index, plot=False, plotRegion=False):
    
    
    #Extract file to use as input in compute_EW()
    jades_file = jades_folder / JADES_files[index]
    mast_file = mast_folder / MAST_files[index]
    z = z_array[index]

    if plot:
        plot_spectrum_full(mast_file, jades_file, index)

    if plotRegion:
        region = plotRegion
        plot_spectrum_region(mast_file, jades_file, index, region)
    
    
    
    # ===== JADES =====
    with fits.open(jades_file) as f:
        specdata = f[1].data
        
        #Flux vector
        flux = specdata['FLUX']  #flux in erg cm-2 s-1 AA-1'
        mask = np.isfinite(flux)
        flux_erg = flux[mask]
        
        #Compute rest wavelength
        lambda_obs = specdata['WAVELENGTH'] * 1e-6 / 1e-10 #convert from um to AA.
        lambda_rest = lambda_obs[mask] / (1.+float(z))
    
        flux = flux_erg * lambda_obs[mask]**2 / 2.99792458e-5
        

        jades_fluxes = get_line_fluxes(lambda_rest, flux, index) 
        
        
    # ===== MAST =====
    try:
        with fits.open(mast_file) as f:
            specdata = f[1].data
            
            #Flux vector
            flux = specdata['FLUX']  #flux in Jy already
            mask = np.isfinite(flux)
            flux = flux[mask]
            
            #Compute rest wavelength
            lambda_obs = specdata['WAVELENGTH'] * 1e-6 / 1e-10 #convert from um to AA.
            lambda_rest = lambda_obs[mask] / (1.+float(z))
        
            
            mast_fluxes = get_line_fluxes(lambda_rest, flux, index)

    except FileNotFoundError:
        print(f'MAST file for Object {ID_array[index]} not found')
        mast_fluxes = np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan
        pass


    # ==== Save the data ====
        
    
    #Save output appending to arrays
    ID_data.append(ID_array[index])
    z_data.append(z_array[index])
    
    
    mast_Ha_data.append( mast_fluxes[0] )
    mast_Ha_err_data.append(mast_fluxes[1] )
    mast_Hb_data.append( mast_fluxes[2] )
    mast_Hb_err_data.append(mast_fluxes[3] )
    mast_O3_data.append(  mast_fluxes[4] )
    mast_O3_err_data.append( mast_fluxes[5] )
    mast_o3_data.append(  mast_fluxes[6] )
    mast_o3_err_data.append( mast_fluxes[7] )


    jades_Ha_data.append( jades_fluxes[0] )
    jades_Ha_err_data.append(jades_fluxes[1] )
    jades_Hb_data.append( jades_fluxes[2] )
    jades_Hb_err_data.append(jades_fluxes[3] )
    jades_O3_data.append(  jades_fluxes[4] )
    jades_O3_err_data.append( jades_fluxes[5] )
    jades_o3_data.append(  jades_fluxes[6] )
    jades_o3_err_data.append( jades_fluxes[7] )

    
    print(f'Saved data for obj with ID={ID_array[index]}')
#------------------------------------------------
#
#------------------------------------------------        
def plot_spectrum_full(mast_file, jades_file, index):

    ID = ID_array[index]
    z = z_array[index]

    plt.figure()
    plt.title(f'ID={ID}, z={z:.2f}')
    plt.ylabel(r' Flux [$\mathrm{erg~s^{-1} cm^{-2} \AA^{-1} }$]')
    plt.xlabel(r' $\lambda_{rest}$ [$\AA$]')
    
    
    file = jades_file
    with fits.open(file) as f:
        specdata = f[1].data
        
        #Flux vector
        flux = specdata['FLUX']  #flux in erg cm-2 s-1 AA-1'
        mask = np.isfinite(flux)
        flux_erg = flux[mask]
        
        #Compute rest wavelength
        lambda_obs = specdata['WAVELENGTH'] * 1e-6 / 1e-10 #convert from um to AA.
        lambda_rest = lambda_obs[mask] / (1.+float(z))
        
        
        plt.step(lambda_rest, flux_erg, label='JADES')
        

    file = mast_file
    with fits.open(file) as f:
        specdata = f[1].data
        
        #Flux vector
        flux = specdata['FLUX']  #mast flux is in Jy
        mask = np.isfinite(flux)
        
        
        #Compute rest wavelength
        lambda_obs = specdata['WAVELENGTH'] * 1e-6 / 1e-10 #convert from um to AA.
        lambda_rest = lambda_obs[mask] / (1.+float(z))
        
        flux_erg = flux[mask] *2.99792458e-5 / lambda_obs[mask]**2
        
    
        plt.step(lambda_rest, flux_erg, label='MAST')
        
    plt.legend(loc='best',)
    plt.show()
#------------------------------------------------
#
#------------------------------------------------
def plot_spectrum_region(mast_file, jades_file, index, region):

    ID = ID_array[index]
    z = z_array[index]

    plt.figure()
    plt.title(f'ID={ID}, z={z:.2f}')
    plt.ylabel(r' Flux [$\mathrm{erg~s^{-1} cm^{-2} \AA^{-1} }$]')
    plt.xlabel(r' $\lambda_{rest}$ [$\AA$]')
    plt.xlim(region[0], region[1])
    
    file = jades_file
    with fits.open(file) as f:
        specdata = f[1].data
        
        #Flux vector
        flux = specdata['FLUX']  #flux in erg cm-2 s-1 AA-1'
        mask = np.isfinite(flux)
        flux_erg = flux[mask]
        
        #Compute rest wavelength
        lambda_obs = specdata['WAVELENGTH'] * 1e-6 / 1e-10 #convert from um to AA.
        lambda_rest = lambda_obs[mask] / (1.+float(z))
        
        
        plt.step(lambda_rest, flux_erg)
        
    file = mast_file
    
    with fits.open(file) as f:
        specdata = f[1].data
        
        #Flux vector
        flux = specdata['FLUX']  #mast flux is in Jy
        mask = np.isfinite(flux)
        
        
        #Compute rest wavelength
        lambda_obs = specdata['WAVELENGTH'] * 1e-6 / 1e-10 #convert from um to AA.
        lambda_rest = lambda_obs[mask] / (1.+float(z))
        
        flux_erg = flux[mask] *2.99792458e-5 / lambda_obs[mask]**2
        
    
        plt.step(lambda_rest, flux_erg)
        
    plt.show()
#------------------------------------------------
#
#------------------------------------------------


main()