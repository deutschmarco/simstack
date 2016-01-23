import pdb
import numpy as np
from astropy.wcs import WCS
from shift import shift_twod
from VieroLibrary.dist_idl import dist_idl
from lmfit import Parameters, minimize, fit_report
from smoothmap import smoothmap
from gauss_kern import gauss_kern
from invert_sed import single_simple_flux_from_greybody

def simultaneous_stack_sed_oned(p, layers_1d, data1d, wavelengths, LenLayers, zed, err1d = None):
  ''' Function to Minimize written specifically for lmfit '''

  v = p.valuesdict()

  nwv = len(wavelengths)
  LenModel = len(data1d) 
  Nlayers = len(layers_1d)/LenModel

  model = np.zeros(LenModel)

  LL = [0]
  SuperChunk = [0]
  LL.extend(LenLayers)
  cuLL = np.cumsum(LL)
  SuperChunk.extend(LenLayers * Nlayers)
  cuSuperChunk = np.cumsum(SuperChunk)
  for i in range(Nlayers):
    Temp = np.asarray(v['T'+str(i)])
    Lir = np.asarray(v['L'+str(i)])

    fluxes = single_simple_flux_from_greybody(np.asarray(wavelengths), Trf = Temp, Lrf = Lir, b=2, zin=zed)
    for iwv in range(nwv):
      model[cuLL[iwv]:cuLL[iwv+1]] +=  fluxes[0][iwv] * layers_1d[ cuSuperChunk[iwv] + i * LenLayers[iwv]: cuSuperChunk[iwv] + (i+1) * LenLayers[iwv] ] 

  if err1d is None:
    return (data1d - model)
  return (data1d - model)/err1d

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
  cmaps, 
  hd, 
  layers_radec, 
  wavelengths,
  fwhm, 
  cnoise=None, 
  mask=None, 
  beam_area=None, 
  err_ss=None, 
  zed=0.01,
  quiet=None):
  
  w = WCS(hd)
  #FIND SIZES OF MAP AND LISTS
  cms = np.shape(cmaps) # should be a cube
  #nwv = np.size(wavelengths)
  nwv = cms[0] 
  #zeromask = np.zeros(cms)

  size_cube = np.shape(layers_radec)
  nsrcmax = size_cube[0]
  nlists = int(size_cube[1])
  
  ind_map_zero = np.where(np.isnan(cmaps))
  nzero = np.shape(ind_map_zero)[1]

  if np.sum(cnoise) == 0: cnoise=cmaps*0.0 + 1.0

  pix=hd["CD2_2"]*3600.
  if pix == 0: pix=hd["CDELT2"]*3600.

  for iwv in range(nwv):
    #[STEP 0] - Calibrate maps
    if beam_area != None:
      cmaps[iwv,:,:]=cmaps[iwv,:,:]*beam_area[iwv]*1e6
      cnoise[iwv,:,:]=noise[iwv,:,:]*beam_area[iwv]*1e6

  # STEP 1  - Make Layers Cube
  layers=np.zeros([nlists,cms[1],cms[2]])

  for s in range(nlists):
    ind_src = np.where(layers_radec[:,s,0] != 0)
    if np.shape(ind_src)[1] > 0:
      ra = layers_radec[ind_src,s,0]
      dec = layers_radec[ind_src,s,1]
      ty,tx = w.wcs_world2pix(ra, dec, 0) 
      # CHECK FOR SOURCES THAT FALL OUTSIDE MAP
      ind_keep = np.where((tx[0] >= 0) & (tx[0] < cms[1]) & (ty[0] >= 0) & (ty[0] < cms[2]))
      nt0 = np.shape(ind_keep)[1]
      real_x=np.floor(tx[0,ind_keep][0]).astype(int)
      real_y=np.floor(ty[0,ind_keep][0]).astype(int)
      # CHECK FOR SOURCES THAT FALL ON ZEROS MAP
      # THIS NEEDS COMPLETE OVERHAUL, PARTICULARLY WHEN INCLUDING DIFFERENT AREA MAPS!!
      if nzero > 0:
        tally = np.zeros(nt0)
        for d in range(nt0):
          if cmaps[0,real_x[d],real_y[d]] != 0: 
            tally[d]=1.
        ind_nz=np.where(tally == 1)
        nt = np.shape(ind_nz)[1]
        real_x = real_x[ind_nz]
        real_y = real_y[ind_nz]
      else: nt = nt0
      for ni in range(nt):
        layers[s, real_x[ni],real_y[ni]]+=1.0

  # STEP 2  - Convolve Layers and put in pixels
  #all_map_layers = np.zeros(np.append(nwv,np.shape(layers)))

  cfits_flat = np.asarray([])
  cfits_flat2= np.asarray([])
  flat_maps= np.asarray([])
  flat_noise= np.asarray([])
  LenLayers= np.zeros([nwv])

  radius = 1.1
  for iwv in range(nwv):
    sig = fwhm[iwv] / 2.355 / pix 
    flattened_pixmap = np.sum(layers,axis=0)
    total_circles_mask = circle_mask(flattened_pixmap, radius * fwhm[iwv], pix)
    ind_fit = np.where(total_circles_mask >= 1) # & zeromask != 0)
    nhits = np.shape(ind_fit)[1]
    #cfits_maps = np.zeros([nlists,nhits])
    LenLayers[iwv] = nhits

    kern = gauss_kern(fwhm[iwv], np.floor(fwhm[iwv] * 10), pix)
    for u in range(nlists):
      layer = layers[u,:,:]  
      tmap = smoothmap(layer, kern)
      tmap -= np.mean(tmap[ind_fit])
      #cfits_maps[u,:] = tmap[ind_fit]
      cfits_flat = np.append(cfits_flat,np.ndarray.flatten(tmap[ind_fit]))

    #FLATTEN EVERYTHING HERE RATHER THAN BELOW
    #cfits_flat2= np.append(cfits_flat,np.ndarray.flatten(cfits_maps))

    #pdb.set_trace()
    lmap = cmaps[iwv]
    lnoise = cnoise[iwv]
    lmap[ind_fit] -= np.mean(lmap[ind_fit], dtype=np.float32)
    flat_maps = np.append(flat_maps,np.ndarray.flatten(lmap[ind_fit]))
    flat_noise = np.append(flat_noise,np.ndarray.flatten(lnoise[ind_fit]))
    #pdb.set_trace()


  # STEP 3 - Regress Layers with Map (i.e., stack!)

  fit_params = Parameters()

  for iarg in range(nlists): 
    fit_params.add('T'+str(iarg),value= 25.,vary=True,min=10.,max=150.)
    fit_params.add('L'+str(iarg),value= 1e12,min=0.,max=1e14)

  cov_ss_1d = minimize(simultaneous_stack_sed_oned, fit_params, 
    args=(cfits_flat,), kws={'data1d':flat_maps,'err1d':flat_noise,'wavelengths':wavelengths,'LenLayers':LenLayers,'zed':zed})
    
  return cov_ss_1d












