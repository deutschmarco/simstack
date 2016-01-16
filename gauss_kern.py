import pdb
def gauss_kern(fwhm, side, pixsize): 
	''' Create a 2D Gaussian (size= side x side)''' 
	from scipy.ndimage.filters import gaussian_filter
	from shift import shift_twod
	import numpy as np
	from numpy import zeros
	from numpy import shape

	sig = fwhm / 2.355 / pixsize
	delt = zeros([side,side])
	delt[0,0]=1.0
	#delt = shift_twod(delt,int(side / 2),int(side / 2))
	ms = shape(delt)
	delt = shift_twod(delt, ms[0] / 2, ms[1] / 2)
	kern = delt
	gaussian_filter(delt, sig, output= kern)
	#gauss = np.exp(-1.0*(k_dist ** 2 / 2.0 / sig ** 2))
	kern /= np.max(kern)

	return kern