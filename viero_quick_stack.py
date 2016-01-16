import pdb
import numpy as np
from VieroLibrary import readcol
from astropy.io import fits
from simstack import stack_in_redshift_slices

def viero_quick_stack(
	map_names, 
	catalog_names, 
	noise_map_names,
	n_sources_max = None
	): 

	if n_sources_max == None: n_sources_max=50000l
	nmap = len(map_names)

	fwhm  = np.array([18.1, 25.2, 36.6])
	efwhm = np.array([17.6, 23.9, 35.2])

	#PUT DATA INTO CUBE
	nlists = len(catalog_names)
	nsources = 0 # initialize a counter
	cube = np.zeros([n_sources_max, nlists, 2]) # nsources by nlis/nts by 2 for RA/DEC
	for i in range(nlists): 
		list_name = catalog_names[i]
		ra, dec = readcol.readcol(list_name,fsep=',',twod=False)
		nsources_list=len(ra)
		if nsources_list > n_sources_max: 
			print 'too many sources in catalog: use N_SOURCES_MAX flag'
			break

		cube[0:nsources_list,i,0]=ra
		cube[0:nsources_list,i,1]=dec
		if nsources_list > nsources: 
			nsources=nsources_list

	cube=cube[0:nsources-1,:,:] # Crop it down to the length of longest list

	stacked_sed=np.zeros([nmap, nlists])
	stacked_sed_err=np.zeros([nmap,nlists])
	#STACK AT ONE WAVELENGTH AT A TIME 

	for wv in range(nmap): 
		#READ MAPS
		cmap, chd = fits.getdata(map_names[wv], 0, header = True)

		#cnoise=0
		if noise_map_names != None: 
			cnoise, nhd = fits.getdata(noise_map_names[wv], 0, header = True)

		#pdb.set_trace()
		stacked_object=stack_in_redshift_slices(
			cmap,
			chd,
			cube,
			efwhm[wv],
			cnoise=cnoise,
		#	err_ss=err_ss,
			quiet=None)
		#pdb.set_trace()

		stacked_flux = np.array(stacked_object.values.values())
		stacked_sed[wv,:] = stacked_flux
			
		#stacked_sed_err[wv,:]=err_ss

		#stacked_flux = None
	return stacked_sed
