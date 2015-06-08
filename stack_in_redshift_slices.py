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
#     CREATED BY marco.viero@caltech.edu     (2013-03-03)
#     ADDED ROUNDING OF FLOATING POING AFTER ADXY  (2014-04-07)
#     ADDED MEAN OF ind_mean SUBSET OF tmap/cmap   (2014-04-07)
#-

import numpy as np

def main(cmap, hd, cube, fwhm, cnoise=0, mask=0, beam_area=0, err_ss=0, quiet=0):
  
  #FIND SIZES OF MAP AND LISTS
  cms = np.shape(cmap)
  zeromask = np.zeros(cms)

  size_cube = np.shape(cube)
  nsrcmax = size_cube[0]
  nlists = size_cube[1]
  
  ind_map_zero = where(np.isnan(cmap) = 'True')

  if cnoise == 0: cnoise=cmap*0.0 + 1.0

  pix=hd.header['CD2_2']*3600.
  if pix == 0: pix=hd.header['CDELT2']*3600.

  #[STEP 0] - Calibrate maps
  if beam_area != 0:
    cmap=cmap*beam_area*1e6
    cnoise=noise*beam_area*1e6

  # STEP 1  - Make Layers Cube
  cube_dimensions=[nlists,ms[0],ms[1]]
  layers = np.zeros(cube_dimensions)
