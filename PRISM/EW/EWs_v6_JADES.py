import math
from scipy import integrate
from scipy.integrate import simps
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate 
import pylab as P
#import mpfit
from scipy import stats
from numpy import array, ndarray, sqrt, sin, sinh, maximum
from scipy.integrate import quad
from math import sqrt as msqrt

from astropy.io import fits
from astropy.table import Table
from astropy import units as u
from matplotlib import pyplot as plt
from astropy.visualization import quantity_support
quantity_support()  # for getting units on the axes below 

from specutils import SpectralRegion
from specutils.analysis import equivalent_width 
from specutils import Spectrum1D, SpectrumCollection
import pandas as pd
		
import warnings
from specutils.fitting import fit_generic_continuum
from specutils.fitting.continuum import fit_continuum


###Leer catalogo de redshift de JADES
t=Table.read('/Users/analuisa/Documents/JWST/candidates/EWs/Clear-Prism_spectra/JADES_Prism/jades_dr3_prism_public_gs_v1.1.fits')
NIRSpec_ID=t['NIRSpec_ID']
z_Spec=t['z_Spec'] #redshift
RA=t['RA_TARG']
DEC=t['Dec_TARG']


z_Spec_2=[]
NIRSpec_ID_2=[]
RA_2=[]
DEC_2=[]
for j in range(len(z_Spec)):
	if float(z_Spec[j])>=3.37:
		z_Spec_2.append(float(z_Spec[j]))
		z_Spec_3=np.asfarray(z_Spec_2)
		NIRSpec_ID_2.append(NIRSpec_ID[j])
		RA_2.append(RA[j])
		DEC_2.append(DEC[j])

import os
folders = os.listdir('/Users/analuisa/Documents/JWST/candidates/EWs/Clear-Prism_spectra/JADES_Prism/hlsp_jades_jwst_nirspec_clear-prism_part_11_01')
# ~ print (folders)

for j in range(len(NIRSpec_ID_2)):
	for i in range(len(folders)):
		# ~ print (folders[i])
		with fits.open('/Users/analuisa/Documents/JWST/candidates/EWs/Clear-Prism_spectra/JADES_Prism/hlsp_jades_jwst_nirspec_clear-prism_part_11_01/' + folders[i], memmap=False, ignore_missing_simple=True) as f:  
			specdata = f[1].data
			hdr = f[0].header
		# ~ lamb = specdata['WAVELENGTH'] * 1e-6/1e-10 * u.AA
		# ~ flux = specdata['FLUX'] * u.Unit('erg cm-2 s-1 AA-1') 
		# ~ spec = Spectrum1D(spectral_axis=lamb, flux=flux)
		target=hdr['TARGET']
		if NIRSpec_ID_2[j]==float(target):
			if 3.37 < z_Spec_2[j] < 6.3:
				lamb = specdata['WAVELENGTH'] * 1e-6/1e-10 / (1.+z_Spec_2[j]) * u.AA
				flux = specdata['FLUX'] * u.Unit('erg cm-2 s-1 AA-1') 
				spec = Spectrum1D(spectral_axis=lamb, flux=flux)
				with warnings.catch_warnings():  # Ignore warnings
					warnings.simplefilter('ignore')
					cont_norm_spec = spec / fit_generic_continuum(spec, exclude_regions=[SpectralRegion(lamb[len(lamb)-1] - (200*u.AA), lamb[len(lamb)-1]), SpectralRegion( lamb[0], lamb[0] + (200*u.AA))])(spec.spectral_axis)
				EW_Ha=equivalent_width(cont_norm_spec, continuum=1, regions=SpectralRegion(6505 * u.AA, 6610 * u.AA)) # este es el rest_EW ya que el espectro esta en rest
				EW_Hb=equivalent_width(cont_norm_spec, continuum=1, regions=SpectralRegion(4815 * u.AA, 4917 * u.AA)) # este es el rest_EW ya que el espectro esta en rest
				EW_OIII=equivalent_width(cont_norm_spec, continuum=1, regions=SpectralRegion(4972 * u.AA, 5062 * u.AA)) # este es el rest_EW ya que el espectro esta en rest
				print (NIRSpec_ID_2[j],z_Spec_2[j],RA_2[j],DEC_2[j],EW_Ha.value,EW_Hb.value,EW_OIII.value)				
				plt.step(spec.spectral_axis, spec.flux)
				plt.ylabel(r' Flux [ergs s$^{-1}$ cm$^{-2}$ $\AA^{-1}$]', fontsize=15)
				plt.xlabel(r' $\lambda_{rest}$ [$\AA$]', fontsize=15)
				plt.title(folders[i]+'_z='+str(z_Spec_2[j]),fontsize=6)
				plt.tick_params(axis='x',labelsize=10)
				plt.tick_params(axis='y',labelsize=10)
				plt.savefig('v3/Plots_Clear-Prism_part_11_01/'+folders[i]+'.pdf')				
				plt.close()	
			else:	
				if 4.8 < z_Spec_2[j] < 9.3:
				# ~ if float(z_Spec[j])!= -1.0:
					# ~ print (NIRSpec_ID[j], target, z_Spec[j])
					lamb = specdata['WAVELENGTH'] * 1e-6/1e-10 / (1.+z_Spec_2[j]) * u.AA
					flux = specdata['FLUX'] * u.Unit('erg cm-2 s-1 AA-1') 
					spec = Spectrum1D(spectral_axis=lamb, flux=flux)
					with warnings.catch_warnings():  # Ignore warnings
						warnings.simplefilter('ignore')
						cont_norm_spec = spec / fit_generic_continuum(spec, exclude_regions=[SpectralRegion(lamb[len(lamb)-1] - (200*u.AA), lamb[len(lamb)-1]), SpectralRegion( lamb[0], lamb[0] + (200*u.AA))])(spec.spectral_axis)
					EW_Hb=equivalent_width(cont_norm_spec, continuum=1, regions=SpectralRegion(4815 * u.AA, 4917 * u.AA)) # este es el rest_EW ya que el espectro esta en rest
					EW_OIII=equivalent_width(cont_norm_spec, continuum=1, regions=SpectralRegion(4972 * u.AA, 5062 * u.AA)) # este es el rest_EW ya que el espectro esta en rest					
					plt.step(spec.spectral_axis, spec.flux)
					plt.ylabel(r' Flux [ergs s$^{-1}$ cm$^{-2}$ $\AA^{-1}$]', fontsize=15)
					plt.xlabel(r' $\lambda_{rest}$ [$\AA$]', fontsize=15)
					plt.title(folders[i]+'_z='+str(z_Spec_2[j]),fontsize=6)
					plt.tick_params(axis='x',labelsize=10)
					plt.tick_params(axis='y',labelsize=10)
					plt.savefig('v3/Plots_Clear-Prism_part_11_01/'+folders[i]+'.pdf')				
					plt.close()
					print (NIRSpec_ID_2[j],z_Spec_2[j],RA_2[j],DEC_2[j],-99.99,EW_Hb.value,EW_OIII.value)				
	


# ~ print (len(z_Spec_3))
# ~ print (DEC_2)

# ~ plt.hist(z_Spec_3, bins=30)
# ~ plt.ylabel(r' N', fontsize=15)
# ~ plt.xlabel(r' redshift', fontsize=15)
# ~ plt.tick_params(axis='x',labelsize=10)
# ~ plt.tick_params(axis='y',labelsize=10)
# ~ plt.text(11,110,'N=634',fontsize=15)
# ~ plt.title('JADES CLEAR/PRISM Spectra, GOODS-S Field, z > 3.4')
# ~ plt.show()
#634 targets con z >=3.37 de los cuales he medido 128 EWs de los cuales aorox. 90 tiene altos EWs, osea un 70%
#osea tendriamos aprox.400 targets con altos EWs
#757 targets con z >=3.0
