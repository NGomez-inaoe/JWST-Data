import math
from scipy import integrate
# ~ from scipy.integrate import simps
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




#archivo de JADES
t_JADES=Table.read('./Spectra1D/JADES/hlsp_jades_jwst_nirspec_goods-s-deephst-00003184_clear-prism_v1.0_x1d.fits', hdu="EXTRACT1D")


w_JADES=t_JADES['WAVELENGTH']
flux_JADES=t_JADES['FLUX'] #flux in Jy


#archivo de mast
# ~ fits.open('file.fits')
t=Table.read('./MAST/jw01210-o001_b000000092_nirspec_f290lp-g395h/jw01210-o001_b000000092_nirspec_f290lp-g395h_x1d.fits', hdu="EXTRACT1D")
w=t['WAVELENGTH']
flux=t['FLUX'] #flux in Jy




z=6.712776
rest_w=w/(1+z) #in micras
rest_w_Angstrom=rest_w*1e-6/1e-10 # rest wavelength in Angstrom
w_Angstrom=w*1e-6/1e-10 # observed wavelength in Angstrom

flux_ergs = 2.99792458E-05 *flux / w_Angstrom**2. #flux in erg s^-1 cm^-2 A^-1

print(flux_ergs)

plt.step(w/1e-4,flux_ergs, label='MAST Data-Archive')
plt.step(w_JADES/1e-4,flux_JADES, label='JADES Survey')

plt.xlabel(r'$\lambda_{obs}$ [$\AA$]', fontsize=22)
plt.ylabel(r'FLUX [erg s$^{-1}$ cm$^{-2}$ $\AA^{-1}$]', fontsize=22)
# ~ legend=plt.legend(loc=4,labelspacing=0.1)
#plt.title(r'ID 1594, z=4.277382', fontsize=18.)
legend=plt.legend(loc='best',labelspacing=0.1)
plt.setp(legend.get_texts(),fontsize='14')
#plt.savefig('FirstPlots/4297.pdf')

