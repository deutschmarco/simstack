import pdb
import numpy as np
import os
from VieroLibrary import readcol
from astropy.io import fits
from simstack import stack_in_redshift_slices as simstack
from sedstack import stack_in_redshift_slices as sedstack
from VieroLibrary.invert_sed import single_simple_flux_from_greybody

def viero_quick_stack(
	map_names, 
	catalog_names, 
	noise_map_names,
	efwhm,
	color_correction=None,
	n_sources_max = None,
	sedfitwavelengths = None,
	zed = 0.001
	): 

	if n_sources_max == None: n_sources_max=50000l
	nmap = len(map_names)

	#fwhm  = np.array([18.1, 25.2, 36.6])
	#efwhm = np.array([17.6, 23.9, 35.2])

	#PUT DATA INTO CUBE
	nlists = len(catalog_names)
	nsources = 0 # initialize a counter
	cube = np.zeros([n_sources_max, nlists, 2]) # nsources by nlis/nts by 2 for RA/DEC
	for i in range(nlists): 
		list_name = catalog_names[i]
		if os.path.getsize(list_name) > 0: 
			ra, dec = readcol.readcol(list_name,fsep=',',twod=False)
			nsources_list=len(ra)
			if nsources_list > n_sources_max: 
				print 'too many sources in catalog: use N_SOURCES_MAX flag'
				break
			if nsources_list > 0:
				cube[0:nsources_list,i,0]=ra
				cube[0:nsources_list,i,1]=dec
			if nsources_list > nsources: 
				nsources=nsources_list

	cube=cube[0:nsources-1,:,:] # Crop it down to the length of longest list

	stacked_sed=np.zeros([nmap, nlists])
	stacked_sed_err=np.zeros([nmap,nlists])

	if sedfitwavelengths != None:
		#FIT SEDS TO FIND FLUXES
		cmaps = [] 
		cnoise = [] 
		for wv in range(nmap): 
			#READ MAPS
			tmaps, thd = fits.getdata(map_names[wv], 0, header = True)
			if color_correction != None:
				tmaps *= color_correction[wv]
			cmaps.append(tmaps)
			if noise_map_names != None: 
				tnoise, nhd = fits.getdata(noise_map_names[wv], 0, header = True)
				if color_correction != None:
					tnoise *= color_correction[wv]
				cnoise.append(tnoise)

		stacked_object=sedstack(
			np.asarray(cmaps),
			thd,
			cube,
			sedfitwavelengths,
			efwhm,
			cnoise=np.asarray(cnoise),
		#	err_ss=err_ss,
			zed=zed,
			quiet=None)

		print 'yeaaahhh!!'
		v = stacked_object.params.valuesdict()
		beta = np.asarray(v['b'])
		for ised in range(nlists):
			Temp = np.asarray(v['T'+str(ised)])
			Lir = np.asarray(v['L'+str(ised)])
			#pdb.set_trace()
			stacked_sed[:,ised]=single_simple_flux_from_greybody(np.asarray(sedfitwavelengths), Trf = Temp, Lrf = Lir, b=beta, zin=zed)
		#pdb.set_trace()
		return [stacked_sed, v]
	else:
		#STACK AT ONE WAVELENGTH AT A TIME 
		for wv in range(nmap): 
			#READ MAPS
			cmap, chd = fits.getdata(map_names[wv], 0, header = True)

			if noise_map_names != None: 
				cnoise, nhd = fits.getdata(noise_map_names[wv], 0, header = True)

			stacked_object=simstack(
				cmap,
				chd,
				cube,
				efwhm[wv],
				cnoise=cnoise,
			#	err_ss=err_ss,
				quiet=None)

			stacked_flux = np.array(stacked_object.params.values())
			stacked_sed[wv,:] = stacked_flux
			#stacked_sed_err[wv,:]=err_ss
		#pdb.set_trace()
		return stacked_sed