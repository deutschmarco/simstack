import pdb
import numpy as np
from astropy.wcs import WCS
from shift import shift_twod
from VieroLibrary.dist_idl import dist_idl
from lmfit import Parameters, minimize, fit_report
from smoothmap import smoothmap
from gauss_kern import gauss_kern

def simultaneous_stack_array_oned(p, layers_1d, data1d, err1d = None):
  ''' Function to Minimize written specifically for lmfit '''

  v = p.valuesdict()

  len_model = len(data1d)
  nlayers = len(layers_1d)/len_model

  model = np.zeros(len_model)

  for i in range(nlayers):
    model[:] += layers_1d[i*len_model:(i+1)*len_model] * v['layer'+str(i)] 

  if err1d is None:
    return (data1d - model)
  return (data1d - model)/err1d

def simultaneous_stack_array(p, layers_2d, data, err = None):
  ''' Function to Minimize written specifically for lmfit '''

  v = p.valuesdict()

  csize = np.shape(layers_2d)

  model = np.zeros(csize[1])

  for i in range(csize[0]):
    model += layers_2d[i,:] * v['layer'+str(i)] 

  if err is None:
    return (data - model)
  return (data - model)/err

def circle_mask(pixmap,radius_in,pixres):
  ''' Makes a 2D circular image of zeros and ones'''

  radius=radius_in/pixres
  xy = np.shape(pixmap)
  xx = xy[0]
  yy = xy[1]
  beforex = np.log2(xx)
  beforey = np.log2(yy)
  if beforex != beforey:
    if beforex > beforey:
      before = beforex 
    else:
      before = beforey
  else: before = beforey
  l2 = np.ceil(before)
  pad_side = 2.0 ** l2
  outmap = np.zeros([pad_side, pad_side])
  outmap[:xx,:yy] = pixmap

  dist_array = shift_twod(dist_idl(pad_side, pad_side), pad_side/2, pad_side/2)
  circ = np.zeros([pad_side, pad_side])
  ind_one = np.where(dist_array <= radius)
  circ[ind_one] = 1.
  mask  = np.real( np.fft.ifft2( np.fft.fft2(circ) *
          np.fft.fft2(outmap)) 
          ) * pad_side * pad_side
  mask = np.round(mask)
  ind_holes = np.where(mask >= 1.0)
  mask = mask * 0.
  mask[ind_holes] = 1.
  maskout = shift_twod(mask, pad_side/2, pad_side/2)

  return maskout[:xx,:yy]

def stack_in_redshift_slices(
  cmap, 
  hd, 
  layers_radec, 
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

  size_cube = np.shape(layers_radec)
  nsrcmax = size_cube[0]
  nlists = int(size_cube[1])
  
  ind_map_zero = np.where(np.isnan(cmap))
  nzero = np.shape(ind_map_zero)[1]

  if np.sum(cnoise) == 0: cnoise=cmap*0.0 + 1.0

  pix=hd["CD2_2"]*3600.
  if pix == 0: pix=hd["CDELT2"]*3600.

  #[STEP 0] - Calibrate maps
  if beam_area != None:
    cmap=cmap*beam_area*1e6
    cnoise=noise*beam_area*1e6

  # STEP 1  - Make Layers Cube
  layers=np.zeros([nlists,cms[0],cms[1]])

  for s in range(nlists):
    ind_src = np.where(layers_radec[:,s,0] != 0)
    #print np.shape(ind_src)[1]
    if np.shape(ind_src)[1] > 0:
      ra = layers_radec[ind_src,s,0]
      dec = layers_radec[ind_src,s,1]
      # CONVERT FROM RA/DEC to X/Y
      # DANGER!!  NOTICE THAT I FLIP X AND Y HERE!! 
      #tx,ty = w.wcs_world2pix(ra, dec, 1)# WHAT IS THE DIFFERENCE BETWEEN 0 AND 1???!!!  
      ty,tx = w.wcs_world2pix(ra, dec, 0)# NOTICE I FLIPPED X AND Y AND NO LONGER TRANSPOSE! 
      # CHECK FOR SOURCES THAT FALL OUTSIDE MAP
      ind_keep = np.where((tx[0] >= 0) & (tx[0] < cms[0]) & (ty[0] >= 0) & (ty[0] < cms[1]))
      nt0 = np.shape(ind_keep)[1]
      real_x=np.floor(tx[0,ind_keep][0]).astype(int)
      real_y=np.floor(ty[0,ind_keep][0]).astype(int)
      # CHECK FOR SOURCES THAT FALL ON ZEROS MAP
      if nzero > 0:
        tally = np.zeros(nt0)
        for d in range(nt0):
          if cmap[real_x[d],real_y[d]] != 0: 
            tally[d]=1.
        ind_nz=np.where(tally == 1)
        nt = np.shape(ind_nz)[1]
        real_x = real_x[ind_nz]
        real_y = real_y[ind_nz]
      else: nt = nt0
      for ni in range(nt):
        layers[s, real_x[ni],real_y[ni]]+=1.0

  # STEP 2  - Convolve Layers and put in pixels
  radius = 1.1
  sig = fwhm / 2.355 / pix 
  flattened_pixmap = np.sum(layers,axis=0)
  total_circles_mask = circle_mask(flattened_pixmap, radius * fwhm, pix)
  ind_fit = np.where(total_circles_mask >= 1) # & zeromask != 0)
  nhits = np.shape(ind_fit)[1]
  cfits_maps = np.zeros([nlists,nhits])

  kern = gauss_kern(fwhm, np.floor(fwhm * 10), pix)
  for u in range(nlists):
    layer = layers[u,:,:] 
    #layer = np.transpose(layers[u,:,:]) ## DANGER!! Transpose NO LONGER required AFTER FLIPPING x and y!!!!!!!
    #tmap = gaussian_filter(layer, sig) 
    tmap = smoothmap(layer, kern)
    tmap -= np.mean(tmap[ind_fit])
    cfits_maps[u,:] = tmap[ind_fit]

  # STEP 3 - Regress Layers with Map (i.e., stack!)

  cmap[ind_fit] -= np.mean(cmap[ind_fit], dtype=np.float32)

  fit_params = Parameters()

  for iarg in range(nlists): 
    fit_params.add('layer'+str(iarg),value= np.random.randn())
  imap = cmap[ind_fit]
  ierr = cnoise[ind_fit]

  #cov_ss = minimize(simultaneous_stack_array, fit_params, args=(cfits_maps,), kws={'data':imap,'err':ierr})
  ##cov_ss = minimize(simultaneous_stack_array, fit_params, args=(cfits_maps,imap,ierr))

  cov_ss_1d = minimize(simultaneous_stack_array_oned, fit_params, 
    args=(np.ndarray.flatten(cfits_maps),), kws={'data1d':np.ndarray.flatten(imap),'err1d':np.ndarray.flatten(ierr)})
    #args=(np.ndarray.flatten(cfits_maps), np.ndarray.flatten(imap),np.ndarray.flatten(ierr)))
  # STEP 4  - Returns a minimizer object

  #pdb.set_trace()
  return cov_ss_1d












