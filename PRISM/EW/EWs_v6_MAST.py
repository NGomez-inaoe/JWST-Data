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

#multiples objetos

f1=open('/Users/analuisa/Documents/JWST/candidates/combined_table.dat',"r")
lines=f1.readlines()
f1.close
T=[line.split() for line in lines]
# ~ print (T)
# ~ T2=np.genfromtxt('/Users/analuisa/Documents/JWST/candidates/combined_table.dat')
T2=np.genfromtxt('/Users/analuisa/Documents/JWST/candidates/combined_table.dat', dtype='unicode')
# ~ print(T2.shape)
# ~ print (T2)

ID_NIRSpec=T2[:,0]	
mast_filename=T2[:,2]
AR=T2[:,3]
DEC=T2[:,4]
z_spec=T2[:,5]	

# ~ print (ID_NIRSpec)
# ~ print (mast_filename)
# ~ print (z_spec)
from specutils import SpectralRegion

# ~ f = open("/Users/analuisa/Documents/JWST/candidates/EWs/candidates_EWs.txt", "w")

for i in range(len(z_spec)):
	# ~ if z_spec[i]!='--':
	if '3.37' < z_spec[i] < '6.3':
		# ~ print (z_spec[i], mast_filename[i])
		with fits.open('/Users/analuisa/Documents/JWST/Spectra/JWST_all_JADES/'+ mast_filename[i] + '/' + mast_filename[i]+'_x1d.fits') as f:  
			specdata = f[1].data
		lamb_obs = specdata['WAVELENGTH'] * 1e-6/1e-10 
		lamb = specdata['WAVELENGTH'] * 1e-6/1e-10 / (1.+float(z_spec[i])) * u.AA
		# ~ flux = specdata['FLUX']*2.99792458E-05 / lamb_obs**2 * u.Unit('erg cm-2 s-1 AA-1') + (2e-19 * u.Unit('erg cm-2 s-1 AA-1'))
		flux = specdata['FLUX']*2.99792458E-05 / lamb_obs**2 * u.Unit('erg cm-2 s-1 AA-1')
		# ~ print (i,mast_filename[i])
		spec = Spectrum1D(spectral_axis=lamb, flux=flux)
		plt.step(spec.spectral_axis, spec.flux)
		plt.show()
		with warnings.catch_warnings():  # Ignore warnings
			warnings.simplefilter('ignore')
			cont_norm_spec = spec / fit_generic_continuum(spec)(spec.spectral_axis)
			# ~ print (fit_generic_continuum(spec)(spec.spectral_axis))
		plt.step(cont_norm_spec.wavelength, cont_norm_spec.flux)  
			# ~ plt.step(cont_norm_spec.wavelength, fit_generic_continuum(spec)(spec.spectral_axis))  
		# ~ plt.xlim(6450 * u.AA, 6650 * u.AA)
		EW_Ha=equivalent_width(cont_norm_spec, continuum=1, regions=SpectralRegion(6550 * u.AA, 6580 * u.AA)) # este es el rest_EW ya que el espectro esta en rest
		# ~ print ('i,mast_filename[i],rest_EW_Ha',i,mast_filename[i],EW_Ha)
		print (ID_NIRSpec[i],mast_filename[i],AR[i],DEC[i],z_spec[i],EW_Ha.value)
		# ~ if EW_Ha.value < -145.:
			# ~ print (ID_NIRSpec[i],mast_filename[i],AR[i],DEC[i],z_spec[i],EW_Ha.value)
			# ~ f.writeto(EW_Ha.value)
		plt.show()
# ~ f.close()

print ('')

for i in range(len(z_spec)):
	# ~ if z_spec[i]!='--':
	if '4.8' < z_spec[i] < '9.3':
		# ~ print (z_spec[i], mast_filename[i])
		with fits.open('/Users/analuisa/Documents/JWST/Spectra/JWST_all_JADES/'+ mast_filename[i] + '/' + mast_filename[i]+'_x1d.fits') as f:  
			specdata = f[1].data
		lamb_obs = specdata['WAVELENGTH'] * 1e-6/1e-10 
		lamb = specdata['WAVELENGTH'] * 1e-6/1e-10 / (1.+float(z_spec[i])) * u.AA
		# ~ flux = specdata['FLUX']*2.99792458E-05 / lamb_obs**2 * u.Unit('erg cm-2 s-1 AA-1') + (2e-19 * u.Unit('erg cm-2 s-1 AA-1'))
		flux = specdata['FLUX']*2.99792458E-05 / lamb_obs**2 * u.Unit('erg cm-2 s-1 AA-1')
		# ~ print (i,mast_filename[i])
		spec = Spectrum1D(spectral_axis=lamb, flux=flux)
		# ~ plt.step(spec.spectral_axis, spec.flux)
		# ~ plt.show()
		with warnings.catch_warnings():  # Ignore warnings
			warnings.simplefilter('ignore')
			# ~ region_Ha=[(6480 * u.AA, 6550 * u.AA), (6575 * u.AA, 6640 * u.AA)]
			cont_norm_spec = spec / fit_generic_continuum(spec)(spec.spectral_axis)
			# ~ cont_norm_spec = spec / fit_continuum(spec, window=region_Ha)(spec.spectral_axis)
			# ~ print (fit_generic_continuum(spec)(spec.spectral_axis))
			# ~ plt.step(cont_norm_spec.wavelength, cont_norm_spec.flux)
			# ~ plt.step(cont_norm_spec.wavelength, fit_generic_continuum(spec)(spec.spectral_axis)) 
			# ~ plt.step(cont_norm_spec.wavelength, fit_continuum(spec, window=region_Ha)(spec.spectral_axis))
		# ~ plt.xlim(6450 * u.AA, 6650 * u.AA)
		# ~ print (cont_norm_spec)
		EW_OIII=equivalent_width(cont_norm_spec, continuum=1, regions=SpectralRegion(4990 * u.AA, 5024 * u.AA)) # este es el rest_EW ya que el espectro esta en rest
		print (ID_NIRSpec[i],mast_filename[i],AR[i],DEC[i],z_spec[i],EW_OIII)
		# ~ if EW_Ha.value < -145.:
			# ~ print (ID_NIRSpec[i],mast_filename[i],AR[i],DEC[i],z_spec[i],EW_Ha.value)
		# ~ plt.show()
