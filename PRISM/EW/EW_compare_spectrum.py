from specutils.fitting import fit_generic_continuum
from specutils import SpectralRegion
from specutils import Spectrum1D
import matplotlib.pyplot as plt
from astropy import units as u
from astropy.io import fits
from pathlib import Path
import pandas as pd
import warnings




#Define source path where the files are
source=Path("/home/nicolas/Documents/Research/PhD/JWST-Data/PRISM/Spectra1D/")
mast_folder=source / "MAST"
jades_folder=source / "JADES"

df = pd.read_csv("short-table.tsv", sep='\t')

ID_array = df["NIRSpec_ID"]
JADES_files = df["JADES_FILENAME"]
MAST_files = df["MAST_FILENAME"]
z_array = df["redshift"]



def compare_spectrum(mast_file, jades_file, z, ID, whichSpectrum):

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
        mast_spectrum = Spectrum1D(spectral_axis=mast_lambda_rest_Angstrom, flux=flux_Units)

        with warnings.catch_warnings(): #ignore warnings
            warnings.simplefilter('ignore')

            # Normalize the spectrum by its continuum
            mast_spec_continuum_fitted = fit_generic_continuum(mast_spectrum)(mast_spectrum.spectral_axis)
            mast_normalized_continuum_spec = mast_spectrum / mast_spec_continuum_fitted



    with fits.open(jades_file) as j_f:
        #Extract data from the file
        specdata=j_f[1].data

        #Compute rest wavelength
        lambda_obs = specdata['WAVELENGTH'] * 1e-6 / 1e-10 #convert from um to AA.
        lambda_rest = lambda_obs / (1.+float(z))
        lambda_rest_Angstrom = lambda_rest * u.AA #u.AA just includes the units

        #Extract and convert flux to working units
        flux = specdata['FLUX']
        flux_Units = flux * u.Unit('erg cm-2 s-1 AA-1')

        #Define the spectrum over which the EW will be calculated, use class Spectrum1D from specutils
        jades_spectrum = Spectrum1D(spectral_axis=lambda_rest_Angstrom, flux=flux_Units)

        with warnings.catch_warnings(): #ignore warnings
            warnings.simplefilter('ignore')

            # Normalize the spectrum by its continuum
            jades_spec_continuum_fitted = fit_generic_continuum(jades_spectrum)(jades_spectrum.spectral_axis)
            jades_normalized_continuum_spec = jades_spectrum / jades_spec_continuum_fitted



    #Plot the spectrum
    plt.figure()


        #Original
    if whichSpectrum == "Ori":
        plotSpectrum_1 = mast_spectrum
        plotSpectrum_2 = jades_spectrum

        #Normalized Continuum
    elif whichSpectrum == "NC":
        plotSpectrum_1 = mast_normalized_continuum_spec
        plotSpectrum_2 = jades_normalized_continuum_spec

        #H alpha
    elif whichSpectrum == "Ha":
        plotSpectrum_1 = mast_normalized_continuum_spec
        plotSpectrum_2 = jades_normalized_continuum_spec
        plt.xlim(6450 , 6650)

        #Hbeta
    elif whichSpectrum == "Hb":
        plotSpectrum_1 = mast_normalized_continuum_spec
        plotSpectrum_2 = jades_normalized_continuum_spec
        plt.xlim(4840 , 4880)


    ##plots

    plt.rcParams["figure.figsize"] = (12, 9)
    plt.step(plotSpectrum_1.spectral_axis, plotSpectrum_1.flux, label='MAST', color='b', linestyle='--')
    plt.step(plotSpectrum_2.spectral_axis, plotSpectrum_2.flux, label='JADES', color='g')

    plt.step(mast_lambda_rest_Angstrom, mast_spec_continuum_fitted, label='MAST NC', color='purple')
    plt.step(lambda_rest_Angstrom, jades_spec_continuum_fitted, label='JADES NC', linestyle='--', color='k')


    #plt.step(normalized_continuum_spec.wavelength, normalized_continuum_spec.flux)
    plt.ylabel(r' Flux [ergs s$^{-1}$ cm$^{-2}$ $\AA^{-1}$]', fontsize=15)
    plt.xlabel(r' $\lambda_{rest}$ [$\AA$]', fontsize=15)
    plt.title(f'ID={ID}, z={z}',fontsize=14)
    legend=plt.legend(loc='best',labelspacing=0.1)
    plt.setp(legend.get_texts(),fontsize='14')
    plt.tick_params(axis='x',labelsize=10)
    plt.tick_params(axis='y',labelsize=10)
    plt.show()
    #plt.close()

def show_spectrum(ID):

    for i in range(len(ID_array)):
        if ID == ID_array[i]:

            jades_file = jades_folder / JADES_files[i]
            mast_file = mast_folder / MAST_files[i]
            z = z_array[i]

            print(f'ID={ID}')
            compare_spectrum(mast_file, jades_file, z, ID, 'Ori')
            compare_spectrum(mast_file, jades_file, z, ID, 'NC')
            compare_spectrum(mast_file, jades_file, z, ID, 'Ha')
            compare_spectrum(mast_file, jades_file, z, ID, 'Hb')


def save_spectrum(ID):
    
    for i in range(len(ID_array)):
        if ID == ID_array[i]:

            jades_file = jades_folder / JADES_files[i]
            mast_file = mast_folder / MAST_files[i]
            z = z_array[i]
            compare_spectrum(mast_file, jades_file, z, ID, 'Ori')
            plt.savefig(f'EW-plots/EW-{ID}')

save_spectrum(7507)
