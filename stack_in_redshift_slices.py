;+
;
; stacked_fluxes = STACK_IN_REDSHIFT_SLICES(cmap, cube, hd=hd, pix=pix, fwhm...)
;
; INPUTS:
;  CMAP : cropped map in Jy/beam.  If in MJy/sr, must supply beam_area
;  HD   : header file, must contain CD2_2 or CDELT2 to determine pixel size
;  CUBE : data cube, N sources per list, M lists, RA and DEC.  Lists will have 
;         different number of sources, so N = longest list, and others have zero
;         in place of values.
;  FWHM : resolution of map.  This is important to get right and a potential source
;         of systematic (i.e., consider simulating to measure potential error)
;
; OPTIONAL OUTPUTS:
;  MASK : Optional places where you don't want to stack
;  BEAM_AREA : To convert from MJy/sr (like spitzer or BLAST) to Jy/beam (like SPIRE)
;
; OUTPUTS:
;   cov_ss: M stacked flux densities 
;
; OPTIONAL OUTPUTS:
;   err_ss: Errors from mpfit (Not very informative)
;
; CREATED BY: mviero 2012-04-30
;
;-

import numpy

fun main:
  ind_map_zero = where(np.isnan(cmap) = 'True')
