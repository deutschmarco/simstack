import pdb
import numpy as np
from astropy.wcs import WCS
import lmfit
from scipy.ndimage.filters import gaussian_filter
from shift import shift_twod
from VieroLibrary.dist_idl import dist_idl
from lmfit import Parameters, minimize, fit_report

def simultaneous_stack_array(p, layers_2d, data = None, err = None):
  ''' Function to Minimize written specifically for lmfit '''

  #import numpy as np

  v = p.valuesdict()

  csize = np.shape(layers_2d)

  model = np.zeros(csize[1])

  for i in range(csize[0]):
    model += layers_2d[i,:] * v['layer'+str(i)] 

  #pdb.set_trace()

  if data is None:
    return model
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

  #print 'oneeeee'
  #pdb.set_trace()
  #[STEP 0] - Calibrate maps
  if beam_area != None:
    cmap=cmap*beam_area*1e6
    cnoise=noise*beam_area*1e6

  # STEP 1  - Make Layers Cube
  layers=np.zeros([nlists,cms[0],cms[1]])

  for s in range(nlists):
    ind_src = np.where(layers_radec[:,s,0] != 0)
    ra = layers_radec[ind_src,s,0]
    dec = layers_radec[ind_src,s,1]
    # CONVERT FROM RA/DEC to X/Y
    tx,ty = w.wcs_world2pix(ra, dec, 1) 
    tx = tx - 1.0
    ty = ty - 1.0
    # CHECK FOR SOURCES THAT FALL OUTSIDE MAP
    ind_keep = np.where((tx[0] >= 0) & (tx[0] < cms[0]) & (ty[0] >= 0) & (ty[0] < cms[1]))
    nt0 = np.shape(ind_keep)[1]
    real_x=np.rint(tx[0,ind_keep][0]).astype(int)
    real_y=np.rint(ty[0,ind_keep][0]).astype(int)
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

  #print 'twooooo'
  #pdb.set_trace()
  # STEP 2  - Convolve Layers and put in pixels
  radius = 1.1
  sig = fwhm / 2.355 / pix 
  flattened_pixmap = np.sum(layers,axis=0)
  #print 'threeee'
  #pdb.set_trace()
  total_circles_mask = circle_mask(flattened_pixmap, radius * fwhm, pix)
  #print 'fourrrrr'
  #pdb.set_trace()
  ind_fit = np.where(total_circles_mask >= 1) # & zeromask != 0)
  nhits = np.shape(ind_fit)[1]
  cfits_maps = np.zeros([nlists,nhits])

  #print 'fiveeeee'
  #pdb.set_trace()
  for u in range(nlists):
    layer = np.transpose(layers[u,:,:]) ## DANGER!!
    tmap = gaussian_filter(layer, sig) 
    tmap -= np.mean(tmap[ind_fit])
    cfits_maps[u,:] = tmap[ind_fit]

  # STEP 3 - Regress Layers with Map (i.e., stack!)
  #p0 = np.ones(nlists)

  cmap[ind_fit] -= np.mean(cmap[ind_fit], dtype=np.float32)

  #cov_ss = mpfitfun('simultaneous_stack_array',cfits_maps, $
  #  cmap[ind_fit], cnoise[ind_fit],p0,ftol=1d-15,quiet=quiet,PERROR=err_ss)

  #print 'fourrrrr'
  #pdb.set_trace()
  fit_params = Parameters()
  for iarg in range(nlists): 
    fit_params.add('layer'+str(iarg),value= 0.001)
  imap = cmap[ind_fit]
  ierr = cnoise[ind_fit]

  #cov_ss = minimize(simultaneous_stack_array, fit_params, args=(cfits_maps), kws={'data':imap,'err':ierr})
  cov_ss = minimize(simultaneous_stack_array, fit_params, args=(cfits_maps,), kws={'data':imap,'err':ierr})
  #cov_ss = minimize(simultaneous_stack_array, fit_params, args=(cfits_maps,imap))

  pdb.set_trace()
  # STEP 4  - Returns a minimizer object

  return cov_ss