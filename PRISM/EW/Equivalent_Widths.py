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



#Create arrays to save the data with save_EW()
ID_data = []
z_data = []
mast_EW_Ha_data = []
jades_EW_Ha_data = []
mast_EW_Hb_data = []
jades_EW_Hb_data = []
mast_EW_O_data = []
jades_EW_O_data = []

def main():
    
    
    #Execute the functions
    
        #print the EW data for a given object in short list with index indx
    #numb = 13
    #get_info(ID_array[numb])
    
    
    
    for n in range(len(ID_array)):
        try:
            save_EW(n)
    
        except FileNotFoundError:
            print(f'File {MAST_files[n]} not found')
            continue
        
        except TypeError:
            print(f'Unkown error for {ID_array[n]}')
            continue
        
        except IndexError:
            print(f'Bad regions for {ID_array[n]}')
    
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
    
    EWD = pd.DataFrame(EW_data)
    EWD.to_csv("EW_output.tsv", sep="\t", index=False)
    
    
    
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
        
        lamb = lambda_rest_Angstrom
        
        with warnings.catch_warnings(): #ignore warnings
            warnings.simplefilter('ignore')
            
            # Normalize the spectrum by its continuum
            spec_continuum_fitted = fit_generic_continuum(mast_spectrum)(mast_spectrum.spectral_axis)
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
            
            # Normalize the spectrum by its continuum
            spec_continuum_fitted = fit_generic_continuum(jades_spectrum)(jades_spectrum.spectral_axis)
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
    
    
    


main()
    