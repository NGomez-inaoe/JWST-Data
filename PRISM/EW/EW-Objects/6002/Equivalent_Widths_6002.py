"""
  This version implements uncertainties from the distribution of the spectrum itself around the lines
  with the std deviation of the flux values, instead of extracting directly from the .fits file.
  Another chane is the computation of the EW for each line with a single function.
  Date: March 8th
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
indx=11

Ha_ini, Ha_end = 6539 *u.AA, 6600 *u.AA
Hb_ini, Hb_end = 4838 *u.AA, 4890 *u.AA
O_ini, O_end = 4931 *u.AA, 5053 *u.AA
#//////////////////////////////////////


def main():


  #Save EW estimate
  save_EW(indx)


  #Save data in data frame
  EW_data = {
      "ID": ID_data,
      "redshift": z_data,
      "EW(Ha) MAST": mast_EW_Ha_data,
      "DEW(Ha) MAST": mast_EW_dHa_data,
      "EW(Ha) JADES": jades_EW_Ha_data,
      "DEW(Ha) JADES": jades_EW_dHa_data,
      "EW(Hb) MAST": mast_EW_Hb_data,
      "DEW(Hb) MAST": mast_EW_dHb_data,
      "EW(Hb) JADES": jades_EW_Hb_data,
      "DEW(Hb) JADES": jades_EW_dHb_data,
      "EW([OIII]) MAST": mast_EW_O_data,
      "EW([OIII]) unc MAST": mast_EW_dO_data,
      "EW([OIII]) JADES": jades_EW_O_data,
      "EW([OIII]) unc JADES": jades_EW_dO_data,
      "Ha_ini": Ha_ini.value,
      "Ha_end": Ha_end.value,
      "Hb_ini": Hb_ini.value,
      "Hb_end": Hb_end.value,
      "[OIII]_ini": O_ini.value,
      "[OIII]_end": O_end.value
    }

  ID = ID_array[indx]
  EWD = pd.DataFrame(EW_data)
  folder=objects_folder / f'{ID}'
  EWD.to_csv(f'{folder}/EW_output_{ID}_2.tsv', sep="\t", index=False)
  EWD.to_csv(f'{output_folder}/EW_output_v6_1.tsv', sep="\t", index=False, mode='a', header=False)

  print(f'Equitalent Widths for obj {ID}, saved!')
  




#//////////////////////////////////////////
#///          Functions                 ///
#//////////////////////////////////////////
def flux_stdDev(lambda_array, flux_array, regionLeft_ini, regionLeft_end, regionRight_ini, regionRight_end):

    #Compute the standard deviation of the flux around a given line
    lambda_array = np.array(lambda_array)
    flux_array = np.array(flux_array)

    # Crear máscara booleana para el rango deseado
    lambda_ini = regionLeft_ini
    lambda_end = regionLeft_end
    mask_left = (lambda_array >= lambda_ini) & (lambda_array <= lambda_end)

    lambda_ini = regionRight_ini
    lambda_end = regionRight_end
    mask_right = (lambda_array >= lambda_ini) & (lambda_array <= lambda_end)

    # Extraer los valores correspondientes
    lambda_left = lambda_array[mask_left]
    flux_left = flux_array[mask_left]

    lambda_right = lambda_array[mask_right]
    flux_right = flux_array[mask_right]

    stdDev_left = np.std(flux_left)
    stdDev_right = np.std(flux_right)
    
    

    return (stdDev_left + stdDev_right)/2



#Compute EW for the given line
def compute_line_EW(lamb, flux, flux_uncertainty, exclusion_regions, line_region):


  spectrum = Spectrum(spectral_axis=lamb,
                              flux=flux,
                              uncertainty=flux_uncertainty)

  with warnings.catch_warnings(): #ignore warnings
    warnings.simplefilter('ignore')

    #//////// Compute continuum by fitting /////////
    try:
      Continuum_fit = fit_generic_continuum(spectrum, exclude_regions=exclusion_regions)(spectrum.spectral_axis)

    except Exception:
      Continuum_fit = None

    #////// Compute Equivalent Width ///////////
    if Continuum_fit is not None:
      # Normalize the spectrum by its continuum 
      Normalized_continuum_spectrum = spectrum / Continuum_fit
    
      ew = equivalent_width(Normalized_continuum_spectrum, continuum=1, regions=line_region)
      EW = ew.value
      EW_err = ew.uncertainty.value
    else:
      EW = np.nan
      EW_err = np.nan

  return EW, EW_err


#Compute the Equivalent Widths for Halpha, Hbeta and [OIII]
def compute_EWs(lambda_rest, flux):

  
  lambda_rest_Angstrom = lambda_rest *u.AA
  
  flux_Jy = flux * u.Jy

  #/////////// Equivalent Widths ////////////

  #//// H alfa

  #Uncertainty of the flux
  flux_err_Ha = flux_stdDev(lambda_rest, flux, Ha_ini.value - 500, Ha_ini.value, Ha_end.value, Ha_end.value + 500 )
  flux_err_Jy = flux_err_Ha * np.ones(len(flux)) * u.Jy
  flux_uncertainty = StdDevUncertainty(flux_err_Jy)
  


  #  Regions to exlude
  lamb = lambda_rest_Angstrom
  Ha_left = SpectralRegion(lamb[0], Ha_ini - 500 *u.AA)
  Ha_region = SpectralRegion(Ha_ini, Ha_end)
  Ha_region_excl = SpectralRegion(Ha_ini, Ha_end + 150 *u.AA )
  Ha_right = SpectralRegion(Ha_end + 350 *u.AA, lamb[-1])
  Ha_exclusion_regions = [Ha_left, Ha_region_excl, Ha_right]

  #// Compute EW of Ha //
  Ha_EW, Ha_EW_err = compute_line_EW(lambda_rest_Angstrom, flux_Jy, flux_uncertainty, Ha_exclusion_regions, Ha_region)

  #//// H beta and OIII

  #Uncertainty of the flux
  flux_err_Hb = flux_stdDev(lambda_rest, flux, Hb_ini.value - 500, Hb_ini.value, O_end.value, O_end.value + 500)
  flux_err_Jy = flux_err_Hb * np.ones(len(flux)) * u.Jy
  flux_uncertainty = StdDevUncertainty(flux_err_Jy)

  #For Hb and [OIII]
  Hb_left = SpectralRegion(lamb[0], Hb_ini - 500 *u.AA)
  Hb_region = SpectralRegion(Hb_ini, Hb_end)
  O3_region = SpectralRegion(O_ini, O_end)
  O3_right = SpectralRegion(O_end + 500 *u.AA, lamb[-1])
  Hb_exclusion_regions = [Hb_left, Hb_region, O3_region, O3_right]

  #// Compute EW of Hb and [OIII] //
  Hb_EW, Hb_EW_err = compute_line_EW(lambda_rest_Angstrom, flux_Jy, flux_uncertainty, Hb_exclusion_regions, Hb_region)
  O3_EW, O3_EW_err = compute_line_EW(lambda_rest_Angstrom, flux_Jy, flux_uncertainty, Hb_exclusion_regions, O3_region)
  

  return Ha_EW, Ha_EW_err, Hb_EW, Hb_EW_err, O3_EW, O3_EW_err



#Use the codes above for the given MAST and JADES file
def EquivalentWidths(mast_file, jades_file, z):


  with fits.open(jades_file) as j_f:
    specdata=j_f[1].data

    #Extract and convert flux from working units to Jy
    flux_erg = specdata['FLUX'] #flux in erg cm-2 s-1 AA-1'
    mask = np.isfinite(flux_erg)
    flux_erg = flux_erg[mask]

    #Compute rest wavelength
    lambda_obs = specdata['WAVELENGTH'] * 1e-6 / 1e-10 #convert from um to AA.
    lambda_rest = lambda_obs[mask] / (1.+float(z))

    flux = flux_erg * lambda_obs[mask]**2 / 2.99792458e-5
    
    jades_EW = compute_EWs(lambda_rest, flux)
  
  
  try:
    with fits.open(mast_file) as m_f:
      #Extract data from the file
      specdata=m_f[1].data

      #Extract flux in Jy from .fits
      flux = specdata['FLUX']

      #Remove NaN values
      mask = np.isfinite(flux)
      flux = flux[mask]    

      #Compute rest wavelength
      lambda_obs = specdata['WAVELENGTH'] * 1e-6 / 1e-10 #convert from um to AA.
      lambda_rest = lambda_obs[mask] / (1.+float(z))
    
      mast_EW = compute_EWs(lambda_rest, flux)
    
  except FileNotFoundError:
    print(f'MAST file for Object {ID_array[indx]} not found')  
    mast_EW = np.nan, np.nan, np.nan, np.nan, np.nan, np.nan
    pass

  return mast_EW, jades_EW

def save_EW(index):

    #Extract file to use as input in compute_EW()
    jades_file = jades_folder / JADES_files[index]
    mast_file = mast_folder / MAST_files[index]
    z = z_array[index]

    mast_output, jades_output = EquivalentWidths(mast_file, jades_file, z)


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

#////////////////////////////////////
#//////    End of Functions    //////
#////////////////////////////////////

main()
