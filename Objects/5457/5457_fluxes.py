import warnings
import numpy as np
import pandas as pd
from pathlib import Path
from astropy.io import fits
from astropy import units as u
from specutils.analysis import line_flux
from astropy.nddata import StdDevUncertainty
from specutils import Spectrum, SpectralRegion


prism_file_jades = Path('/home/nicolas/Documents/Research/PhD/JWST-Data/PRISM/Spectra1D/JADES/hlsp_jades_jwst_nirspec_goods-s-deephst-00005457_clear-prism_v1.0_x1d.fits')
#prism_file_mast = Path('/home/nicolas/Documents/Research/PhD/JWST-Data/PRISM/Spectra1D/MAST/jw01210-o001_s000000054_nirspec_clear-prism_x1d.fits')
g235m_file_jades = Path('/home/nicolas/Documents/Research/PhD/JWST-Data/Medium_Resolution/Spectra1D/JADES/hlsp_jades_jwst_nirspec_goods-s-deephst-00005457_f170lp-g235m_v1.0_x1d.fits')
#g235m_file_mast = Path('/home/nicolas/Documents/Research/PhD/JWST-Data/Mediun_Resolution/Spectra1D/MAST/jw01210-o001_s000000054_nirspec_f170lp-g235m_x1d.fits')

z = 4.8631511339788
ID = '5457'



def main():

    
    jades_prism= flux_prism()
    jades_g235m = flux_g235m()
    
    Hbeta_prism, Hbeta_prism_unc, O3_prism, O3_prism_unc = jades_prism
    O5007_g395h, O5007_g395h_unc, O4959_g395h, O4959_g395h_unc = jades_g235m

    

    only_5007_jades, only_5007_jades_unc = compute_5007_flux(O3_prism, O5007_g395h, O4959_g395h, O3_prism_unc, O5007_g395h_unc, O4959_g395h_unc)

    print(f"Hbeta flux from JADES: {Hbeta_prism:.4e} +/- {Hbeta_prism_unc:.4e} erg/s/cm^2")
    print(f"Only 5007 flux from JADES: {only_5007_jades:.4e} +/- {only_5007_jades_unc:.4e} erg/s/cm^2")

    """
    Hbeta_prism, Hbeta_prism_unc, O3_prism, O3_prism_unc = mast_prism
    O5007_g395h, O5007_g395h_unc, O4959_g395h, O4959_g395h_unc = mast_g235m

    only_5007_mast, only_5007_mast_unc = compute_5007_flux(O3_prism, O5007_g395h, O4959_g395h, O3_prism_unc, O5007_g395h_unc, O4959_g395h_unc)

    
    print(f"Hbeta flux from MAST: {Hbeta_prism:.4e} +/- {Hbeta_prism_unc:.4e} erg/s/cm^2")
    print(f"Only 5007 flux from MAST: {only_5007_mast:.4e} +/- {only_5007_mast_unc:.4e} erg/s/cm^2")
    """


def compute_5007_flux(O3_prism, O5007_g395h, O4959_g395h, O3_prism_unc, O5007_g395h_unc, O4959_g395h_unc):

    #Compute flux for only 5007
    P = O3_prism
    R = O5007_g395h / O4959_g395h

    only_5007_flux = P / (1 + R)
    only_5007_flux_unc = only_5007_flux * np.sqrt((O3_prism_unc / O3_prism)**2 + (O5007_g395h_unc / O5007_g395h)**2 + (O4959_g395h_unc / O4959_g395h)**2)

    return only_5007_flux, only_5007_flux_unc
    
#=======================================================================
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
def compute_line_flux(lamb, flux, region, z, filter):
    
    r_ini = region[0] 
    r_end = region[1] 
    vlight = 2.99792458e-5 
    
    with warnings.catch_warnings(): #ignore warnings
        warnings.simplefilter('ignore')


        if np.isfinite(r_ini) and np.isfinite(r_end):

            #Flux Uncerainty
            if filter=='PRISM':
                flux_err = flux_stdDev(lamb.value, flux.value, [4800 - 200, 4800], [5060, 5060 + 200]) 
            else:
                flux_err = flux_stdDev(lamb.value, flux.value, [4850 - 50, 4850], [5015, 5015 + 50]) 

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


    return np.abs(l_flux*1e19), l_flux_unc*1e19
#------------------------------------------------
#
#------------------------------------------------
def get_line_fluxes_prism(lambda_rest, flux):

     #Arrays for the lines
    
    Hb_region = [4818, 4905] 
    O3_region = [4934, 5065]
    
    
    flux_Jy = flux * u.Jy
    lambda_rest_AA = lambda_rest *u.AA
    
    
    Hb_flux, Hb_flux_unc = compute_line_flux(lambda_rest_AA, flux_Jy, Hb_region, z, filter='PRISM')
    O3_flux, O3_flux_unc = compute_line_flux(lambda_rest_AA, flux_Jy, O3_region, z, filter='PRISM')
    
    


    return Hb_flux, Hb_flux_unc, O3_flux, O3_flux_unc
#------------------------------------------------
#
#------------------------------------------------
def flux_prism():
    

    with fits.open(prism_file_jades) as f:
        specdata = f[1].data

        #Flux vector
        flux = specdata['FLUX']  #flux in erg cm-2 s-1 AA-1'
        mask = np.isfinite(flux)
        flux_erg = flux[mask]
        
        #Compute rest wavelength
        lambda_obs = specdata['WAVELENGTH'] * 1e-6 / 1e-10 #convert from um to AA.
        lambda_rest = lambda_obs[mask] / (1.+float(z))
    
        flux = flux_erg * lambda_obs[mask]**2 / 2.99792458e-5
        

        jades_fluxes = get_line_fluxes_prism(lambda_rest, flux) 

    
    return jades_fluxes
#------------------------------------------------------------------------
#=======================================================================
#------------------------------------------------------------------------
def get_line_fluxes_g235m(lambda_rest, flux):

     #Arrays for the lines
    
    O4959_region = [4954, 4964] 
    O5007_region = [4996, 5015]
    
    
    flux_Jy = flux * u.Jy
    lambda_rest_AA = lambda_rest *u.AA
    
    
    O4959_flux, O4959_flux_unc = compute_line_flux(lambda_rest_AA, flux_Jy, O4959_region, z, 'G235M')
    O5007_flux, O5007_flux_unc = compute_line_flux(lambda_rest_AA, flux_Jy, O5007_region, z, 'G235M')
    
    


    return O4959_flux, O4959_flux_unc, O5007_flux, O5007_flux_unc
#------------------------------------------------
#
#------------------------------------------------
def flux_g235m():
    

    with fits.open(g235m_file_jades) as f:
        specdata = f[1].data

        #Flux vector
        flux = specdata['FLUX']  #flux in erg cm-2 s-1 AA-1'
        mask = np.isfinite(flux)
        flux_erg = flux[mask]
        
        #Compute rest wavelength
        lambda_obs = specdata['WAVELENGTH'] * 1e-6 / 1e-10 #convert from um to AA.
        lambda_rest = lambda_obs[mask] / (1.+float(z))
    
        flux = flux_erg * lambda_obs[mask]**2 / 2.99792458e-5
        

        jades_fluxes = get_line_fluxes_g235m(lambda_rest, flux) 

    return jades_fluxes

#=======================================================================
main()