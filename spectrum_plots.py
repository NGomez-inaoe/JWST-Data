import warnings
from isort import file
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
path = Path('/home/nicolas/Documents/Research/PhD/JWST-Data/')
prism_folder = path / 'PRISM/Spectra1D'
medium_folder = path / 'Medium_Resolution/Spectra1D'
high_folder = path / 'High_Resolution/Spectra1D'

#Extract data with Pandas
df = pd.read_csv("list_of_candidates.tsv", sep='\t')
ID_array = df["NIRSpec_ID"]
z_array = df["redshift"]


def main():

    index = 62
    filter = 'PRISM'

    plt.figure()

    plot_spectrum(source='MAST', index=index, filter=filter, frame='rest')
    plot_spectrum(source='JADES', index=index, filter=filter, frame='rest')
    
    plt.show()

#================================================================
#====================== Plotting functions ======================
#================================================================
def plot_spectrum(source, index, filter='PRISM', frame='rest'):

    
    if filter == 'PRISM':
        print('prism')
        folder = prism_folder
    elif filter == 'G235M' or filter == 'G395M':
        print('medium')
        folder = medium_folder
    elif filter == 'G235H' or filter == 'G395H':
        folder = high_folder
    else:
        raise ValueError("Unknown Filter ")
    

    files_list = df[f"{source}_FILENAME_{filter}"]
        
    file = folder / source / files_list[index]  
    ID = ID_array[index]
    z = z_array[index]


    with fits.open(file) as f:
        #Extract data from the file
        specdata=f[1].data

        lambda_obs = specdata['WAVELENGTH'] * 1e-6 / 1e-10 #convert from um to AA.
        flux = specdata['FLUX'] * 1e19
        if source == 'MAST':
            flux = flux * 2.99792458E-05 /lambda_obs**2

        if frame == 'rest': 
            
            lambda_rest = lambda_obs / (1.+float(z))
            wl = lambda_rest
        else:
            wl = lambda_obs
        
    
    plt.title(f'ID={ID}, z={z:.3f}')
    plt.ylabel(r' Flux [$10^{-19} ~ \mathrm{erg~s^{-1} cm^{-2} \AA^{-1} }$]')
    plt.xlabel(r' $\lambda_{rest}$ [$\AA$]')
    plt.tick_params(axis='x',labelsize=10)
    plt.tick_params(axis='y',labelsize=10)
    
    plt.step(wl, flux, label=source)
    


"""
def plot_medium_resolution(source, index, filter='G235M',frame='rest'):

    
    files_list = df[f"{source}_FILENAME_{filter}"]
    
    file = medium_folder / source / files_list[index]  
    ID = ID_array[index]
    z = z_array[index]

    

    with fits.open(file) as f:
        #Extract data from the file
        specdata=f[1].data

        lambda_obs = specdata['WAVELENGTH'] * 1e-6 / 1e-10 #convert from um to AA.
        flux = specdata['FLUX'] * 1e19
        if source == 'MAST':
            flux = flux * 2.99792458E-05 /lambda_obs**2

        if frame == 'rest': 
            
            lambda_rest = lambda_obs / (1.+float(z))
            wl = lambda_rest
        else:
            wl = lambda_obs
        
    
    plt.title(f'ID={ID}, z={z:.3f}')
    plt.ylabel(r' Flux [$10^{-19} \, \mathrm{erg~s^{-1} cm^{-2} \AA^{-1} }$]')
    plt.xlabel(r' $\lambda_{rest}$ [$\AA$]')
    plt.step(wl, flux, label=source)




def plot_spectrum(mast_file, jades_file, z, ID):


    with fits.open(mast_file) as m_f:
        #Extract data from the file
        specdata=m_f[1].data

        #Compute rest wavelength
        lambda_obs = specdata['WAVELENGTH'] * 1e-6 / 1e-10 #convert from um to AA.
        lambda_rest = lambda_obs / (1.+float(z))
        mast_lambda_rest_Angstrom = lambda_rest * u.AA #u.AA just includes the units

        #Extract and convert flux to working units
        flux = specdata['FLUX'] * 2.99792458E-05 /lambda_obs**2
        flux_Units = flux * u.Unit('erg cm-2 s-1 AA-1')


        #Define the spectrum over which the EW will be calculated, use class Spectrum1D from specutils
        mast_spectrum = Spectrum(spectral_axis=mast_lambda_rest_Angstrom, flux=flux_Units)

        with warnings.catch_warnings(): #ignore warnings
            warnings.simplefilter('ignore')

            # Normalize the spectrum by its continuum
                #Regions to exclude
            lamb = mast_lambda_rest_Angstrom
            First_region = SpectralRegion(lamb[0], lamb_ini)
            Last_region = SpectralRegion(lamb[-1] - lamb_end, lamb[-1])
            All_regions = [First_region, A_region, B_region, C_region, D_region, E_region, Last_region]
            mast_spec_continuum_fitted = fit_generic_continuum(mast_spectrum, exclude_regions=All_regions)(mast_spectrum.spectral_axis)



    with fits.open(jades_file) as j_f:
        #Extract data from the file
        specdata=j_f[1].data

        #Compute rest wavelength
        lambda_obs = specdata['WAVELENGTH'] * 1e-6 / 1e-10 #convert from um to AA.
        lambda_rest = lambda_obs / (1.+float(z))
        jades_lambda_rest_Angstrom = lambda_rest * u.AA #u.AA just includes the units

        #Extract and convert flux to working units
        flux = specdata['FLUX']
        flux_Units = flux * u.Unit('erg cm-2 s-1 AA-1')

        #Define the spectrum over which the EW will be calculated, use class Spectrum1D from specutils
        jades_spectrum = Spectrum1D(spectral_axis=jades_lambda_rest_Angstrom, flux=flux_Units)

        with warnings.catch_warnings(): #ignore warnings
            warnings.simplefilter('ignore')

            #Regions to exclude
            lamb=jades_lambda_rest_Angstrom
            First_region = SpectralRegion(lamb[0], lamb_ini)
            Last_region = SpectralRegion(lamb[-1] - lamb_end, lamb[-1])
            All_regions = [First_region, A_region, B_region, C_region, D_region, E_region, Last_region]
            jades_spec_continuum_fitted = fit_generic_continuum(jades_spectrum, exclude_regions=All_regions)(jades_spectrum.spectral_axis)



    ## //// Plot the spectrum /////
    plt.figure()
    output_folder=objects_folder / f'{ID}'

            #Only JADES
    plt.rcParams["figure.figsize"] = (12, 9)

    plt.step(jades_spectrum.spectral_axis, jades_spectrum.flux, label='JADES', color='g')
    plt.step(jades_lambda_rest_Angstrom, jades_spec_continuum_fitted, label='JADES NC',  color='black')


    #plt.step(normalized_continuum_spec.wavelength, normalized_continuum_spec.flux)
    plt.ylabel(r' Flux [ergs s$^{-1}$ cm$^{-2}$ $\AA^{-1}$]', fontsize=15)
    plt.xlabel(r' $\lambda_{rest}$ [$\AA$]', fontsize=15)
    plt.title(f'ID={ID}, z={z}',fontsize=14)
    legend=plt.legend(loc='best',labelspacing=0.1)
    plt.setp(legend.get_texts(),fontsize='14')
    plt.tick_params(axis='x',labelsize=10)
    plt.tick_params(axis='y',labelsize=10)
    plt.savefig(f'{output_folder}/EW-{ID}-J.pdf')
    plt.close()


        #Only MAST
    plt.rcParams["figure.figsize"] = (12, 9)

    plt.step(mast_spectrum.spectral_axis, mast_spectrum.flux, label='MAST', color='b')
    plt.step(mast_lambda_rest_Angstrom, mast_spec_continuum_fitted, label='MAST NC', color='black')


    plt.ylabel(r' Flux [ergs s$^{-1}$ cm$^{-2}$ $\AA^{-1}$]', fontsize=15)
    plt.xlabel(r' $\lambda_{rest}$ [$\AA$]', fontsize=15)
    plt.title(f'ID={ID}, z={z}',fontsize=14)
    legend=plt.legend(loc='best',labelspacing=0.1)
    plt.setp(legend.get_texts(),fontsize='14')
    plt.tick_params(axis='x',labelsize=10)
    plt.tick_params(axis='y',labelsize=10)
    
    


       
    #plt.savefig(f'{output_folder}/EW-{ID}-B.pdf')
    plt.show()
    



def plot_spectrum(source='MAST', filter='PRISM', ID='3184'):

    if source == 'MAST':
        folder = mast_folder
    elif source == 'JADES':
        folder = jades_folder
    else:
        raise ValueError("Source must be either 'MAST' or 'JADES'")
    



    jades_folder = source / 'JADES'
    mast_folder = source / 'MAST'
"""
#================================================================
#====================== _End of functions_ ======================
#================================================================
main()