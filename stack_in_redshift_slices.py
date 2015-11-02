#+
#  STACK_IN_REDSHIFT_SLICES,  cmap, hd, cube, fhwm,$
#           cnoise=cnoise,mask=mask,beam_area=beam_area,$
#           err_ss=err_ss, quiet=quiet
#
#  PURPOSE: Takes input map with header and cube containing layers and regresses 
#     with MPFITFUN to find the best stacked solution.  Noise and masking 
#     are both optional.
#
#  INPUTS:  
#     cmap : cropped map in Jy/beam.  If in MJy/sr, must supply beam_area
#     hd   : header file, must contain CD2_2 or CDELT2 to determine pixel size
#     cube : data cube, N sources per list, M lists, RA and DEC.  Lists will have 
#      different number of sources, so N = longest list, and others have zero
#      in place of values.
#     fwhm : resolution of map.  This is important to get right and a potential source
#      of systematic uncertainty (i.e., consider simulating to measure potential error)
#
#  OPTIONAL OUTPUTS:
#     mask : Optional places where you don't want to stack
#     beam_area : To convert from MJy/sr (like spitzer or BLAST) to Jy/beam (like SPIRE)
#     err_ss: Errors from mpfit 
#
#  OUTPUTS:
#     cov_ss: M stacked flux densities 
#
#  USAGE:
#     stacked_fluxes = STACK_IN_REDSHIFT_SLICES(cmap,hd,cube,fwhm,cnoise=cnoise,mask=mask,quiet=1)
#
#  CALLS:
#     NAN_FINDER
#     GAUSS_KERN 
#     CIRCLE_MASK
#     MPFITFUN
#     SIMULTANEOUS_STACK_ARRAY
#     
#
#  HISTORY:
#     CREATED FOR IDL BY marco.viero@caltech.edu     (2013-03-03)
#     TRANSLATED TO PYTHONG BY marco.viero@stanford.edu (2015-07)
#-

import pdb
import numpy as np
from astropy.wcs import WCS
from astropy.io import fits
#import mpfit
#import numpy.oldnumeric as Numeric
import lmfit

def stack_in_redshift_slices(
  cmap, 
  hd, 
  cube, 
  fwhm, 
  cnoise=None, 
  mask=None, 
  beam_area=None, 
  err_ss=None, 
  quiet=None):
  
  w = WCS(hd)
  #FIND SIZES OF MAP AND LISTS
  cms = np.shape(cmap)
  zeromask = np.zeros(cms)

  size_cube = np.shape(cube)
  nsrcmax = size_cube[0]
  nlists = size_cube[1]
  
  ind_map_zero = np.where(np.isnan(cmap) == 'True')
  nzero = len(ind_map_zero)

  if cnoise == 0: cnoise=cmap*0.0 + 1.0

  pix=hd["CD2_2"]*3600.
  if pix == 0: pix=hd["CDELT2"]*3600.

  #[STEP 0] - Calibrate maps
  if beam_area != None:
    cmap=cmap*beam_area*1e6
    cnoise=noise*beam_area*1e6

  # STEP 1  - Make Layers Cube
  cube_dimensions=[nlists,cms[0],cms[1]]
  layers = np.zeros(cube_dimensions)

  for s in range(nlists):
    ind_src = np.where(cube[:,s,0] != 0)
    ra = cube[ind_src,s,0]
    dec = cube[ind_src,s,1]
    # CONVERT FROM RA/DEC to X/Y
    tx,ty = w.wcs_world2pix(ra, dec, 1)
    # CHECK FOR SOURCES THAT FALL OUTSIDE MAP
    #ind_keep = np.where(tx[0,:] >= 0 & tx[0,:] < cms[0] & ty[0,:] >= 0 & ty[0,:] < cms[1])
    ind_keep = np.where((tx >= 0) & (tx < cms[0]) & (ty >= 0) & (ty < cms[1])
    nt0 = len(ind_keep)
    real_x=np.round(tx[ind_keep])
    real_y=np.round(ty[ind_keep])
    # CHECK FOR SOURCES THAT FALL ON ZEROS MAP
    if nzero > 0:
      tally = np.zeros(nt0)
      for d in range(nt0):
        if cmap[real_x[d],real_y[d]] != 0: 
          tally[d]=1.
      ind_nz=np.where(tally == 1)
      nt = len(ind_nz)
      real_x = real_x[ind_nz]
      real_y = real_y[ind_nz]
    else: nt = nt0
    #if keyword_set(verbose) then $
    #  print, 'of '+ strcompress(string(nt,format='(i10)'),/remove_all)+' sources in list '+strcompress(string(nt0-nt,format='(i10)'),/remove_all)+' fall outside'
    for ni in range(nt):
      layers[s, real_x[ni],real_y[ni]]+=1.0

  # STEP 2  - Convolve Layers and put in pixels
  radius = 1.1
  #kern = gauss_kern(fwhm, size, pixsize)
  sig = fwhm / 2.355 / pixsize 
  total_circles_mask = circle_mask(flattened_pixmap, radius * fwhm, pixsize)
  ind_fit = np.where(total_circles_mask >= 1 & zeromask != 0)
  nhits = len(ind_fit)
  cfits_maps = np.zeros([nlists,nhits])
  for u in range(nlists):
    #layer = np.sum(layers, axis = 0)
    layer = layers[u,:,:]
    tmap = gaussian_filter(layer, sig) 
    cfits_maps[u,:] = tmap[ind_fit]

  # STEP 3 - Regress Layers with Map (i.e., stack!)
  p0 = np.ones(nlists)
  cmap -= np.mean(cmap[ind_fit], dtype=np.float32)

  #cov_ss = mpfitfun('simultaneous_stack_array',cfits_maps, $
  #  cmap[ind_fit], cnoise[ind_fit],p0,ftol=1d-15,quiet=quiet,PERROR=err_ss)

  # STEP 4  - Bob's your uncle 
  return cov_ss
