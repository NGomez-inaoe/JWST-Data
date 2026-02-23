from specutils.fitting import fit_generic_continuum
from specutils.analysis import equivalent_width
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



def main():

    show_spectrum(7507)






# ======================================
# Functions
# ======================================

def plot_spectra(mast_file, jades_file, z, ID):

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



    fig, axs = plt.subplots(2, 2, figsize=(16, 12))
    axs = axs.flatten()  # para indexar fácil

    cases = ["Ori", "NC", "Ha", "Hb"]

    for ax, whichSpectrum in zip(axs, cases):

        # Selección de espectro
        if whichSpectrum == "Ori":
            plotSpectrum_1 = mast_spectrum
            plotSpectrum_2 = jades_spectrum

        elif whichSpectrum == "NC":
            plotSpectrum_1 = mast_normalized_continuum_spec
            plotSpectrum_2 = jades_normalized_continuum_spec

        elif whichSpectrum == "Ha":
            plotSpectrum_1 = mast_normalized_continuum_spec
            plotSpectrum_2 = jades_normalized_continuum_spec
            ax.set_xlim(6450, 6650)

        elif whichSpectrum == "Hb":
            plotSpectrum_1 = mast_normalized_continuum_spec
            plotSpectrum_2 = jades_normalized_continuum_spec
            ax.set_xlim(4840, 4880)

        # Plots
        ax.step(plotSpectrum_1.spectral_axis,
                plotSpectrum_1.flux,
                label='MAST',
                linestyle='--',
                color='b')

        ax.step(plotSpectrum_2.spectral_axis,
                plotSpectrum_2.flux,
                label='JADES',
                color='g')

        ax.step(mast_lambda_rest_Angstrom,
                mast_spec_continuum_fitted,
                label='MAST NC',
                linestyle='--',
                color='cyan')

        ax.step(lambda_rest_Angstrom,
                jades_spec_continuum_fitted,
                label='JADES NC',
                color='lightgreen')

        ax.set_title(whichSpectrum)
        ax.set_xlabel(r'$\lambda_{rest}$ [$\AA$]')
        ax.set_ylabel(r'Flux')

        ax.tick_params(axis='both', labelsize=9)

    # Leyenda solo en el primer panel para no saturar
    axs[0].legend(fontsize=10)

    fig.suptitle(f'ID={ID}, z={z}', fontsize=16)
    fig.tight_layout()
    plt.show()



def show_spectrum(ID):

    for i in range(len(ID_array)):
        if ID == ID_array[i]:

            jades_file = jades_folder / JADES_files[i]
            mast_file = mast_folder / MAST_files[i]
            z = z_array[i]

            print(f'ID={ID}')
            plot_spectra(mast_file, jades_file, z, ID)


main()
