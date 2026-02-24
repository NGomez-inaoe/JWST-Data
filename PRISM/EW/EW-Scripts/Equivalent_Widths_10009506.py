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
source=Path("/home/nicolas/Documents/Research/PhD/JWST-Data/PRISM/")
mast_folder=source / "Spectra1D/MAST"
jades_folder=source / "Spectra1D/JADES"
output_folder=source / "EW/EW-Data-Output"
plots_folder=source / "EW/EW-Plots"

df = pd.read_csv("short-table.tsv", sep='\t')

ID_array = df["NIRSpec_ID"]
JADES_files = df["JADES_FILENAME"]
MAST_files = df["MAST_FILENAME"]
z_array = df["redshift"]



#Create arrays to save the data with save_EW()
ID_data = []
z_data = []
mast_EW_Ha_data = []
jades_EW_Ha_data = []
mast_EW_Hb_data = []
jades_EW_Hb_data = []
mast_EW_O_data = []
jades_EW_O_data = []

lamb_ini = 2500 * u.AA
#For this object we also exclude final regions:
lamb_end = 1 * u.AA

def main():


    #Execute the functions

      

    #Save EW estimate
    n = 49
    save_EW(n)

    
    EW_data = {
        "ID": ID_data,
        "redshift": z_data,
        "EW(Ha) MAST": mast_EW_Ha_data,
        "EW(Ha) JADES": jades_EW_Ha_data,
        "EW(Hb) MAST": mast_EW_Hb_data,
        "EW(Hb) JADES": jades_EW_Hb_data,
        "EW([OIII]) MAST": mast_EW_O_data,
        "EW([OIII]) JADES": jades_EW_O_data
        }

    ID = ID_data[0]
    EWD = pd.DataFrame(EW_data)
    
    EWD.to_csv(f'{output_folder}/EW_output_{ID}.tsv', sep="\t", index=False)
    
    #Save spectrum with fitted continuum
    save_spectrum(ID)



#///////////////////////////////////////
#       Functions                    ///
#///////////////////////////////////////
def compare_EW(mast_file, jades_file, z, printValue=False):


    with fits.open(mast_file) as m_f:
        #Extract data from the file
        specdata=m_f[1].data

        #Compute rest wavelength
        lambda_obs = specdata['WAVELENGTH'] * 1e-6 / 1e-10 #convert from um to AA.
        lambda_rest = lambda_obs / (1.+float(z))
        lambda_rest_Angstrom = lambda_rest * u.AA #u.AA just includes the units

        #Extract and convert flux to working units
        flux = specdata['FLUX'] * 2.99792458E-05 /lambda_obs**2
        flux_Units = flux * u.Unit('erg cm-2 s-1 AA-1')


        #Define the spectrum over which the EW will be calculated, use class Spectrum1D from specutils
        mast_spectrum = Spectrum1D(spectral_axis=lambda_rest_Angstrom, flux=flux_Units)


        with warnings.catch_warnings(): #ignore warnings
            warnings.simplefilter('ignore')

                #Initial region to exlude
            lamb = lambda_rest_Angstrom
            Lower_lamb_region = SpectralRegion(lamb[0], lamb_ini)
            Last_lamb_region = SpectralRegion(lamb[-1] - lamb_end, lamb[-1])
                #Compute continuum by fitting
            spec_continuum_fitted = fit_generic_continuum(mast_spectrum, exclude_regions=[Lower_lamb_region, Last_lamb_region])(mast_spectrum.spectral_axis)
                # Normalize the spectrum by its continuum
            mast_normalized_continuum_spec = mast_spectrum / spec_continuum_fitted

        #Define Spectral Regions for the emission lines
        region_Ha = SpectralRegion(6530 * u.AA, 6660 * u.AA) #This is in rest frame since lambda=lambda_rest
        region_Hb = SpectralRegion(4815 * u.AA, 4906 * u.AA)
        region_OIII = SpectralRegion(4906 * u.AA, 5058 * u.AA)


        # Compute EW for Ha, Hb and [OIII]
        mast_EW_Ha=equivalent_width(mast_normalized_continuum_spec, regions=region_Ha)
        mast_EW_Hb=equivalent_width(mast_normalized_continuum_spec, regions=region_Hb)
        mast_EW_OIII=equivalent_width(mast_normalized_continuum_spec, continuum=1, regions=region_OIII)


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

            lamb = lambda_rest_Angstrom
            Lower_lamb_region = SpectralRegion(lamb[0], lamb_ini)
            Last_lamb_region = SpectralRegion(lamb[-1] - lamb_end, lamb[-1])
                #Compute continuum by fitting
            spec_continuum_fitted = fit_generic_continuum(jades_spectrum, exclude_regions=[Lower_lamb_region, Last_lamb_region])(jades_spectrum.spectral_axis)
            jades_normalized_continuum_spec = jades_spectrum / spec_continuum_fitted

        #Define Spectral Regions for the emission lines
        region_Ha = SpectralRegion(6530 * u.AA, 6660 * u.AA) #This is in rest frame since lambda=lambda_rest
        region_Hb = SpectralRegion(4815 * u.AA, 4906 * u.AA)
        region_OIII = SpectralRegion(4906 * u.AA, 5058 * u.AA)


        # Compute EW for Ha, Hb and [OIII]
        jades_EW_Ha=equivalent_width(jades_normalized_continuum_spec, regions=region_Ha)
        jades_EW_Hb=equivalent_width(jades_normalized_continuum_spec, regions=region_Hb)
        jades_EW_OIII=equivalent_width(jades_normalized_continuum_spec, continuum=1, regions=region_OIII)


    #Print to the screen the computed values
    if printValue:
        print("EW from mast:\n",
              "EW(Ha) =", mast_EW_Ha.value, "\n"
              "EW(Hb) =", mast_EW_Hb.value, "\n"
              "EW([OIII]) =", mast_EW_OIII.value, "\n")

        print('\n')

        print("EW from jades file:\n",
              "EW(Ha) =", jades_EW_Ha.value, "\n"
              "EW(Hb) =", jades_EW_Hb.value, "\n"
              "EW([OIII]) =", jades_EW_OIII.value, "\n")


    return mast_EW_Ha.value, mast_EW_Hb.value, mast_EW_OIII.value, jades_EW_Ha.value, jades_EW_Hb.value, jades_EW_OIII.value


def get_info(ID):

    for i in range(len(ID_array)):
        if ID == ID_array[i]:

            jades_file = jades_folder / JADES_files[i]
            mast_file = mast_folder / MAST_files[i]
            z = z_array[i]

            print(f'ID={ID}')
            compare_EW(mast_file, jades_file, z, printValue=True)



def save_EW(i):

    jades_file = jades_folder / JADES_files[i]
    mast_file = mast_folder / MAST_files[i]
    z = z_array[i]

    ma_Ha, ma_Hb, ma_O, jd_Ha, jd_Hb, jd_O = compare_EW(mast_file, jades_file, z)

    ID_data.append(ID_array[i])
    z_data.append(z_array[i])
    mast_EW_Ha_data.append(ma_Ha)
    jades_EW_Ha_data.append(jd_Ha)
    mast_EW_Hb_data.append(ma_Hb)
    jades_EW_Hb_data.append(jd_Hb)
    mast_EW_O_data.append(ma_O)
    jades_EW_O_data.append(jd_O)




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
            lamb=mast_lambda_rest_Angstrom
            Lower_lamb_region = SpectralRegion(lamb[0], lamb_ini)
            Last_lamb_region = SpectralRegion(lamb[-1] - lamb_end, lamb[-1])
            mast_spec_continuum_fitted = fit_generic_continuum(mast_spectrum, exclude_regions=[Lower_lamb_region, Last_lamb_region])(mast_spectrum.spectral_axis)
            mast_normalized_continuum_spec = mast_spectrum / mast_spec_continuum_fitted



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

            # Normalize the spectrum by its continuum
            lamb=jades_lambda_rest_Angstrom
            Lower_lamb_region = SpectralRegion(lamb[0], lamb_ini)
            Last_lamb_region = SpectralRegion(lamb[-1] - lamb_end, lamb[-1])
            jades_spec_continuum_fitted = fit_generic_continuum(jades_spectrum, exclude_regions=[Lower_lamb_region, Last_lamb_region])(jades_spectrum.spectral_axis)
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
    plt.step(plotSpectrum_1.spectral_axis, plotSpectrum_1.flux, label='MAST', color='b')
    plt.step(plotSpectrum_2.spectral_axis, plotSpectrum_2.flux, label='JADES', color='g')

    plt.step(mast_lambda_rest_Angstrom, mast_spec_continuum_fitted, label='MAST NC', color='cyan')
    plt.step(jades_lambda_rest_Angstrom, jades_spec_continuum_fitted, label='JADES NC',  color='lightgreen')


    #plt.step(normalized_continuum_spec.wavelength, normalized_continuum_spec.flux)
    plt.ylabel(r' Flux [ergs s$^{-1}$ cm$^{-2}$ $\AA^{-1}$]', fontsize=15)
    plt.xlabel(r' $\lambda_{rest}$ [$\AA$]', fontsize=15)
    plt.title(f'ID={ID}, z={z}',fontsize=14)
    legend=plt.legend(loc='best',labelspacing=0.1)
    plt.setp(legend.get_texts(),fontsize='14')
    plt.tick_params(axis='x',labelsize=10)
    plt.tick_params(axis='y',labelsize=10)
    
    #plt.close()



def save_spectrum(ID):
    
    for i in range(len(ID_array)):
        if ID == ID_array[i]:

            jades_file = jades_folder / JADES_files[i]
            mast_file = mast_folder / MAST_files[i]
            z = z_array[i]
            compare_spectrum(mast_file, jades_file, z, ID, 'Ori')
            plt.savefig(f'{plots_folder}/EW-{ID}.pdf')
            plt.show()

main()

