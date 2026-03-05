"""
in this new version: Equivalent widths are calculated
taking only the continuum around the lines of interes intead to the whole
spectrum
"""
from specutils.fitting import fit_generic_continuum
from specutils.analysis import equivalent_width
from astropy.nddata import StdDevUncertainty
from specutils import SpectralRegion
from specutils import Spectrum
import matplotlib.pyplot as plt
from astropy import units as u
from astropy.io import fits
from pathlib import Path
import seaborn as sns
import pandas as pd
import numpy as np
import warnings



#Define source path where the files are
source=Path("/home/nicolas/Documents/Research/PhD/JWST-Data/PRISM/")
mast_folder=source / "Spectra1D/MAST"
jades_folder=source / "Spectra1D/JADES"
output_folder=source / "EW/EW-Data-Output"
plots_folder=source / "EW/EW-Plots"
objects_folder=source / "EW/EW-Objects"

#Extract data with Pandas
df = pd.read_csv("short-table.tsv", sep='\t')
ID_array = df["NIRSpec_ID"]
JADES_files = df["JADES_FILENAME"]
MAST_files = df["MAST_FILENAME"]
z_array = df["redshift"]

#Create arrays to save the data with save_EW()
ID_data =[]
z_data = []
mast_EW_Ha_data = []
mast_EW_dHa_data =[]
mast_EW_Hb_data = []
mast_EW_dHb_data =[]
mast_EW_O_data = []
mast_EW_dO_data =[]

jades_EW_Ha_data = []
jades_EW_dHa_data =[]
jades_EW_Hb_data = []
jades_EW_dHb_data= []
jades_EW_O_data =[]
jades_EW_dO_data=[]



#Define Spectral Regions for the emission lines: this is the region
#we edit
#//////////////////////////////////////
indx=19

Ha_ini, Ha_end = 6510 *u.AA, 6608 *u.AA
Hb_ini, Hb_end = 4821 *u.AA, 4901 *u.AA
O_ini, O_end = 4918 *u.AA, 5060 *u.AA
#//////////////////////////////////////



def main():

    #Save EW estimate
    save_EW(indx)
    save_spectra(indx)

    #Save data in data frame
    EW_data = {
        "ID": ID_data,
        "redshift": z_data,
        "EW(Ha) MAST": mast_EW_Ha_data,
        "EW(Ha) unc MAST": mast_EW_dHa_data,
        "EW(Ha) JADES": jades_EW_Ha_data,
        "EW(Ha) unc JADES": jades_EW_dHa_data,
        "EW(Hb) MAST": mast_EW_Hb_data,
        "EW(Hb) unc MAST": mast_EW_dHb_data,
        "EW(Hb) JADES": jades_EW_Hb_data,
        "EW(Hb) unc JADES": jades_EW_dHb_data,
        "EW([OIII]) MAST": mast_EW_O_data,
        "EW([OIII]) unc MAST": mast_EW_dO_data,
        "EW([OIII]) JADES": jades_EW_O_data,
        "EW([OIII]) unc JADES": jades_EW_dO_data
    }

    #Save the EW estimations in external file
    ID = ID_array[indx]
    EWD = pd.DataFrame(EW_data)
    folder=objects_folder / f'{ID}'
    EWD.to_csv(f'{folder}/EW_output_{ID}.tsv', sep="\t", index=False)
    EWD.to_csv(f'{output_folder}/EW_output_v5.tsv', sep="\t", index=False, mode='a', header=False)


    print(f'Equitalent Widths for obj {ID}, saved!')




#//////////////////////////////////////////
#///       Functions                    ///
#//////////////////////////////////////////
def compute_EW(mast_file, jades_file, z):

    #///////// Extract MAST DATA /////////

    with fits.open(mast_file) as m_f:
        #Extract data from the file
        specdata=m_f[1].data

        #Compute rest wavelength
        lambda_obs = specdata['WAVELENGTH'] * 1e-6 / 1e-10 #convert from um to AA.
        lambda_rest = lambda_obs / (1.+float(z))
        lambda_rest_Angstrom = lambda_rest * u.AA #u.AA just includes the units

        #Extract flux in Jy from .fits
        flux = specdata['FLUX']
        flux_Jy = flux * u.Jy

        #Uncertainty of the flux
        flux_err = specdata['FLUX_ERROR']
        flux_err_Jy = flux_err * u.Jy
        flux_uncertainty = StdDevUncertainty(flux_err_Jy)

        #Define the spectrum over which the EW will be calculated, use class Spectrum1D from specutils
        spectrum = Spectrum(spectral_axis=lambda_rest_Angstrom,
                                 flux=flux_Jy,
                                 uncertainty=flux_uncertainty)


        #Remove NaN values
        mask = np.isfinite(spectrum.flux.value)
        clean_spectrum = Spectrum(spectral_axis=lambda_rest_Angstrom[mask],
                              flux=flux_Jy[mask],
                              uncertainty=flux_uncertainty[mask])


        with warnings.catch_warnings(): #ignore warnings
            warnings.simplefilter('ignore')

            #///////// Regions to exlude /////////

            lamb = lambda_rest_Angstrom
            #For Ha
            Ha_left = SpectralRegion(lamb[0], Ha_ini - 500 *u.AA)
            Ha_region = SpectralRegion(Ha_ini, Ha_end)
            Ha_region_excl = SpectralRegion(Ha_ini, Ha_end + 150 *u.AA )
            Ha_right = SpectralRegion(Ha_end + 350 *u.AA, lamb[-1])
            Ha_exclusion_regions = [Ha_left, Ha_region_excl, Ha_right]

            #For Hb and [OIII]
            Hb_left = SpectralRegion(lamb[0], Hb_ini - 500 *u.AA)
            Hb_region = SpectralRegion(Hb_ini, Hb_end)
            OIII_region = SpectralRegion(O_ini, O_end)
            OIII_right = SpectralRegion(O_end + 500 *u.AA, lamb[-1])
            Hb_exclusion_regions = [Hb_left, Hb_region, OIII_region, OIII_right]

            #//////// Compute continuum by fitting /////////

            Ha_continuum_fit = fit_generic_continuum(clean_spectrum, exclude_regions=Ha_exclusion_regions)(clean_spectrum.spectral_axis)
            Hb_continuum_fit = fit_generic_continuum(clean_spectrum, exclude_regions=Hb_exclusion_regions)(clean_spectrum.spectral_axis)


            #////// Normalize the spectrum by its continuum ///////////

            Ha_normalized_continuum_spec = clean_spectrum / Ha_continuum_fit
            Hb_normalized_continuum_spec = clean_spectrum / Hb_continuum_fit

            #////// Compute Equivalent Width ///////////

            Ha_EW = equivalent_width(Ha_normalized_continuum_spec, continuum=1, regions=Ha_region)
            Hb_EW = equivalent_width(Hb_normalized_continuum_spec, continuum=1, regions=Hb_region)
            OIII_EW = equivalent_width(Hb_normalized_continuum_spec, continuum=1, regions=OIII_region)

            # Save the results
            mast_output = [Ha_EW.value, Ha_EW.uncertainty, Hb_EW.value, Hb_EW.uncertainty, OIII_EW.value, OIII_EW.uncertainty]


        #///////// Extract JADES DATA /////////

    with fits.open(jades_file) as j_f:
        #Extract data from the file
        specdata=j_f[1].data

        #Compute rest wavelength
        lambda_obs = specdata['WAVELENGTH'] * 1e-6 / 1e-10 #convert from um to AA.
        lambda_rest = lambda_obs / (1.+float(z))
        lambda_rest_Angstrom = lambda_rest * u.AA #u.AA just includes the units

        #Extract and convert flux from working units to Jy
        flux = specdata['FLUX'] #flux in erg cm-2 s-1 AA-1'
        flux_Jy = flux * lambda_obs**2 / 2.99792458e-5 *u.Jy

        #Uncertainty of the flux
        flux_err = specdata['FLUX_ERR'] * lambda_obs**2 / 2.99792458e-5
        flux_err_Jy = flux_err * u.Jy
        flux_uncertainty = StdDevUncertainty(flux_err_Jy)

        #Define the spectrum over which the EW will be calculated, use class Spectrum1D from specutils
        spectrum = Spectrum(spectral_axis=lambda_rest_Angstrom,
                                 flux=flux_Jy,
                                 uncertainty=flux_uncertainty)


        #Remove NaN values
        mask = np.isfinite(spectrum.flux.value)
        clean_spectrum = Spectrum(spectral_axis=lambda_rest_Angstrom[mask],
                              flux=flux_Jy[mask],
                              uncertainty=flux_uncertainty[mask])


        with warnings.catch_warnings(): #ignore warnings
            warnings.simplefilter('ignore')

            #///////// Regions to exlude /////////

            lamb = lambda_rest_Angstrom
            #For Ha
            Ha_left = SpectralRegion(lamb[0], Ha_ini - 500 *u.AA)
            Ha_region = SpectralRegion(Ha_ini, Ha_end)
            Ha_region_excl = SpectralRegion(Ha_ini, Ha_end + 150 *u.AA )
            Ha_right = SpectralRegion(Ha_end + 350 *u.AA, lamb[-1])
            Ha_exclusion_regions = [Ha_left, Ha_region_excl, Ha_right]

            #For Hb and [OIII]
            Hb_left = SpectralRegion(lamb[0], Hb_ini - 500 *u.AA)
            Hb_region = SpectralRegion(Hb_ini, Hb_end)
            OIII_region = SpectralRegion(O_ini, O_end)
            OIII_right = SpectralRegion(O_end + 500 *u.AA, lamb[-1])
            Hb_exclusion_regions = [Hb_left, Hb_region, OIII_region, OIII_right]

            #//////// Compute continuum by fitting /////////

            Ha_continuum_fit = fit_generic_continuum(clean_spectrum, exclude_regions=Ha_exclusion_regions)(clean_spectrum.spectral_axis)
            Hb_continuum_fit = fit_generic_continuum(clean_spectrum, exclude_regions=Hb_exclusion_regions)(clean_spectrum.spectral_axis)


            #////// Normalize the spectrum by its continuum ///////////

            Ha_normalized_continuum_spec = clean_spectrum / Ha_continuum_fit
            Hb_normalized_continuum_spec = clean_spectrum / Hb_continuum_fit

            #////// Compute Equivalent Width ///////////

            Ha_EW = equivalent_width(Ha_normalized_continuum_spec, continuum=1, regions=Ha_region)
            Hb_EW = equivalent_width(Hb_normalized_continuum_spec, continuum=1, regions=Hb_region)
            OIII_EW = equivalent_width(Hb_normalized_continuum_spec, continuum=1, regions=OIII_region)

            # Save the results
            jades_output = [Ha_EW.value, Ha_EW.uncertainty, Hb_EW.value, Hb_EW.uncertainty, OIII_EW.value, OIII_EW.uncertainty]


    return mast_output, jades_output




def save_EW(index):

    #Extract file to use as input in compute_EW()
    jades_file = jades_folder / JADES_files[index]
    mast_file = mast_folder / MAST_files[index]
    z = z_array[index]

    #Get the output from compute_()
    mast_output, jades_output = compute_EW(mast_file, jades_file, z)

    #Save output appending to arrays
    ID_data.append(ID_array[index])
    z_data.append(z_array[index])

    mast_EW_Ha_data.append( mast_output[0] )
    mast_EW_dHa_data.append(mast_output[1] )
    mast_EW_Hb_data.append( mast_output[2] )
    mast_EW_dHb_data.append(mast_output[3] )
    mast_EW_O_data.append(  mast_output[4] )
    mast_EW_dO_data.append( mast_output[5] )

    jades_EW_Ha_data.append( jades_output[0] )
    jades_EW_dHa_data.append(jades_output[1] )
    jades_EW_Hb_data.append( jades_output[2] )
    jades_EW_dHb_data.append(jades_output[3] )
    jades_EW_O_data.append(  jades_output[4] )
    jades_EW_dO_data.append( jades_output[5] )



def plot_spectra(mast_file, jades_file, z, ID):


    #//////// The same as before, now in order to plot the spectrum /////////////////
    #///////// Extract MAST DATA /////////

    with fits.open(mast_file) as m_f:
        #Extract data from the file
        specdata=m_f[1].data

        #Compute rest wavelength
        mast_lambda_obs = specdata['WAVELENGTH'] * 1e-6 / 1e-10 #convert from um to AA.
        lambda_rest = mast_lambda_obs / (1.+float(z))
        lambda_rest_Angstrom = lambda_rest * u.AA #u.AA just includes the units


        #Extract flux in Jy from .fits
        flux = specdata['FLUX']
        flux_Jy = flux * u.Jy


        #Uncertainty of the flux
        flux_err = specdata['FLUX_ERROR']
        flux_err_Jy = flux_err * u.Jy
        flux_uncertainty = StdDevUncertainty(flux_err_Jy)

        #Define the spectrum over which the EW will be calculated, use class Spectrum1D from specutils
        spectrum = Spectrum(spectral_axis=lambda_rest_Angstrom,
                                 flux=flux_Jy,
                                 uncertainty=flux_uncertainty)


        #Remove NaN values
        mast_mask = np.isfinite(spectrum.flux.value)
        mast_spectrum = Spectrum(spectral_axis=lambda_rest_Angstrom[mast_mask],
                              flux=flux_Jy[mast_mask],
                              uncertainty=flux_uncertainty[mast_mask])




        with warnings.catch_warnings(): #ignore warnings
            warnings.simplefilter('ignore')

            #///////// Regions to exlude /////////

            lamb = lambda_rest_Angstrom
            #For Ha
            Ha_left = SpectralRegion(lamb[0], Ha_ini - 500 *u.AA)
            Ha_region = SpectralRegion(Ha_ini, Ha_end)
            Ha_region_excl = SpectralRegion(Ha_ini, Ha_end + 150 *u.AA )
            Ha_right = SpectralRegion(Ha_end + 350 *u.AA, lamb[-1])
            Ha_exclusion_regions = [Ha_left, Ha_region_excl, Ha_right]

            #For Hb and [OIII]
            Hb_left = SpectralRegion(lamb[0], Hb_ini - 500 *u.AA)
            Hb_region = SpectralRegion(Hb_ini, Hb_end)
            OIII_region = SpectralRegion(O_ini, O_end)
            OIII_right = SpectralRegion(O_end + 500 *u.AA, lamb[-1])
            Hb_exclusion_regions = [Hb_left, Hb_region, OIII_region, OIII_right]

            #//////// Compute continuum by fitting /////////

            mast_Ha_continuum_fit = fit_generic_continuum(mast_spectrum, exclude_regions=Ha_exclusion_regions)(mast_spectrum.spectral_axis)
            mast_Hb_continuum_fit = fit_generic_continuum(mast_spectrum, exclude_regions=Hb_exclusion_regions)(mast_spectrum.spectral_axis)


        #///////// Extract JADES DATA /////////

    with fits.open(jades_file) as j_f:
        #Extract data from the file
        specdata=j_f[1].data

        #Compute rest wavelength
        jades_lambda_obs = specdata['WAVELENGTH'] * 1e-6 / 1e-10 #convert from um to AA.
        lambda_rest = jades_lambda_obs / (1.+float(z))
        lambda_rest_Angstrom = lambda_rest * u.AA #u.AA just includes the units

        #Extract and convert flux from working units to Jy
        flux = specdata['FLUX'] #flux in erg cm-2 s-1 AA-1'
        flux_Jy = flux * jades_lambda_obs**2 / 2.99792458e-5 *u.Jy

        #Uncertainty of the flux
        flux_err = specdata['FLUX_ERR'] * jades_lambda_obs**2 / 2.99792458e-5
        flux_err_Jy = flux_err * u.Jy
        flux_uncertainty = StdDevUncertainty(flux_err_Jy)

        #Define the spectrum over which the EW will be calculated, use class Spectrum1D from specutils
        spectrum = Spectrum(spectral_axis=lambda_rest_Angstrom,
                                 flux=flux_Jy,
                                 uncertainty=flux_uncertainty)


        #Remove NaN values
        mask = np.isfinite(spectrum.flux.value)
        jades_spectrum = Spectrum(spectral_axis=lambda_rest_Angstrom[mask],
                              flux=flux_Jy[mask],
                              uncertainty=flux_uncertainty[mask])
        jades_lambda = jades_spectrum.spectral_axis

        with warnings.catch_warnings(): #ignore warnings
            warnings.simplefilter('ignore')

            #///////// Regions to exlude /////////

            lamb = lambda_rest_Angstrom
            #For Ha
            Ha_left = SpectralRegion(lamb[0], Ha_ini - 500 *u.AA)
            Ha_region = SpectralRegion(Ha_ini, Ha_end)
            Ha_region_excl = SpectralRegion(Ha_ini, Ha_end + 150 *u.AA )
            Ha_right = SpectralRegion(Ha_end + 350 *u.AA, lamb[-1])
            Ha_exclusion_regions = [Ha_left, Ha_region_excl, Ha_right]

            #For Hb and [OIII]
            Hb_left = SpectralRegion(lamb[0], Hb_ini - 500 *u.AA)
            Hb_region = SpectralRegion(Hb_ini, Hb_end)
            OIII_region = SpectralRegion(O_ini, O_end)
            OIII_right = SpectralRegion(O_end + 500 *u.AA, lamb[-1])
            Hb_exclusion_regions = [Hb_left, Hb_region, OIII_region, OIII_right]

            #//////// Compute continuum by fitting /////////

            jades_Ha_continuum_fit = fit_generic_continuum(jades_spectrum, exclude_regions=Ha_exclusion_regions)(jades_spectrum.spectral_axis)
            jades_Hb_continuum_fit = fit_generic_continuum(jades_spectrum, exclude_regions=Hb_exclusion_regions)(jades_spectrum.spectral_axis)



    """
    //////////////////////////////////////////////////////////////////
    ///////////////     Plots of the spectrum       //////////////////
    //////////////////////////////////////////////////////////////////
    """

    # ///////////////////  Plot Setting  ////////////////////////////
    output_folder=objects_folder / f'{ID}'


    #General Settings
    plt.rc('axes', titlesize=18)     # fontsize of the axes title
    plt.rc('axes', labelsize=18)     # fontsize of the x and y labels
    plt.rc('axes', linewidth=1.65 )   # width of the frame
    plt.rc('xtick', labelsize=16)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=16)    # fontsize of the tick labels
    plt.rc('legend', fontsize=12)    # legend fontsize
    plt.rc('font', size=18)          # controls default text sizes



    # Latex rendering
    plt.rc('font', **{'family': 'serif', 'serif': ['palatino']}) #Helvetica
    plt.rc('text', usetex=True)
    plt.rc('font', weight='bold')


    # --- Cosmology palette (normalized RGB) ---
    COSMO_COLORS = {
        "deep_blue": (30/255, 60/255, 120/255),
        "warm_gold": (200/255, 140/255, 0/255),
        "dark_teal": (0/255, 110/255, 110/255),
        "soft_crimson": (150/255, 40/255, 40/255),
        "cool_gray": (120/255, 120/255, 120/255),
        }

    sns.set_theme(
        style="ticks",
        context='notebook',
        rc={
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
            "axes.grid": True,
        })

    sns.set_palette([
        COSMO_COLORS["deep_blue"],
        "black",
        COSMO_COLORS["warm_gold"],
        COSMO_COLORS["cool_gray"],
        COSMO_COLORS["dark_teal"],
        COSMO_COLORS["soft_crimson"],
        ])



    #///////////////////////////// Figures ////////////////////
    """
    plt.figure()

    #Plot for Ha region, only MAST
    plt.rcParams["figure.figsize"] = (10, 10)
    plt.ylabel(r' Flux [$\mathrm{Jy}$]')
    plt.xlabel(r' $\lambda_{rest}$ [$\AA$]')
    plt.title(f'ID={ID}, z={z}')
    plt.tick_params(axis='x')
    plt.tick_params(axis='y')
    plt.minorticks_on()

    plt.step(mast_spectrum.spectral_axis, mast_spectrum.flux, label=r'MAST $(H\alpha)$')
    plt.step(mast_lambda, mast_Ha_continuum_fit, label='Continuum Fit (MAST)')
    plt.xlim(Ha_ini.value - 500 , Ha_end.value + 350 )
    plt.ylim(bottom=-1e-7)
    legend=plt.legend(loc='best',labelspacing=0.1);
    plt.setp(legend.get_texts());
    plt.savefig(f'{output_folder}/EW_{ID}_Ha_mast.pdf')

    #plt.show()
    plt.close()


    #Plot for Hb region, only MAST
    plt.figure()
    plt.rcParams["figure.figsize"] = (10, 10)
    plt.ylabel(r' Flux [$\mathrm{Jy}$]')
    plt.xlabel(r' $\lambda_{rest}$ [$\AA$]')
    plt.title(f'ID={ID}, z={z}')
    plt.tick_params(axis='x')
    plt.tick_params(axis='y')
    plt.minorticks_on()

    plt.step(mast_spectrum.spectral_axis, mast_spectrum.flux, label=r'MAST $H\beta$')
    plt.step(mast_lambda, mast_Hb_continuum_fit, label='Continuum Fit (MAST)')
    legend=plt.legend(loc='best',labelspacing=0.1)
    plt.setp(legend.get_texts())
    plt.xlim(Hb_ini.value - 500 , O_end.value + 500 )
    plt.ylim(bottom=-1e-7)
    plt.savefig(f'{output_folder}/EW_{ID}_Hb_mast.pdf')
    #plt.show()
    plt.close()

    """
    #Plot for all the spectrum, both sources
    plt.figure()
    plt.rcParams["figure.figsize"] = (10, 10)
    plt.ylabel(r' Flux [$\mathrm{erg s^{-1} cm^{-2} \AA}$]')
    plt.xlabel(r' $\lambda_{rest}$ [$\AA$]')
    plt.title(f'ID={ID}, z={z}')
    plt.tick_params(axis='x')
    plt.tick_params(axis='y')
    plt.minorticks_on()

    jades_flux_ergs = jades_spectrum.flux * 2.99792458e-5 / jades_lambda_obs[mask]**2
    mast_flux_ergs = mast_spectrum.flux * 2.99792458e-5 / mast_lambda_obs[mast_mask]**2

    plt.step(mast_spectrum.spectral_axis, mast_flux_ergs, label='MAST', color='black')
    plt.step(jades_spectrum.spectral_axis, jades_flux_ergs, label='JADES', color=COSMO_COLORS['warm_gold'])

    legend=plt.legend(loc='best',labelspacing=0.1)
    plt.setp(legend.get_texts(),fontsize='14')
    plt.savefig(f'{output_folder}/EW_{ID}_All.pdf')
    plt.show()

    #///////////////////////////// End of Figures ////////////////////



def save_spectra(index):
    #Extract file to use as input plot_spectra
    jades_file = jades_folder / JADES_files[index]
    mast_file = mast_folder / MAST_files[index]
    z = z_array[index]
    ID = ID_array[index]

    plot_spectra(mast_file, jades_file, z, ID)

#////////////////////////////////////
#//////    End of Functions    //////
#////////////////////////////////////

main()
