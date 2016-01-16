import pdb
import pylab as plt
from viero_quick_stack import viero_quick_stack

map_path='/Users/marco/Documents/Publications/my_old_publications/k_band_stacking_paper/lorenzo_code/SPIRE/stacking_code/test_data/'
catalog_path=map_path
nlists=3
list_names=[]
for a in range(3):
	list_names.append(catalog_path+'test_cat_'+str(a)+'.csv')

map_suffix=['PSW','PMW','PLW']
wavelength=[250,350,500]
fwhm =[18.1, 25.2, 36.6]
efwhm=[17.6, 23.9, 35.2] # want to the measured effective FWHM later

map_files=[]
noise_files=[]
for m in map_suffix: map_files.append(map_path+'uds_cutout_flux_' + m + '.fits')
for m in map_suffix: noise_files.append(map_path+'uds_cutout_noise_' + m + '.fits')

stacked_fluxes =  None
n_sources_max = None
stacked_fluxes = viero_quick_stack(
	map_files, 
	list_names,
	noise_files)

wv1=0 # here wv1=0 means start at 250
wv2=2 # and wv2=2 means end at 500
nwv=wv2-wv1+1
#xt=textoidl('\lambda [\mum]')
#yt='stacked flux [Jy]'
#xr=[100,1000]
#yr=[1e-3,3e-2]
#lot,[0],[0],xr=xr,yr=yr,/xs,/ys,/xl,/yl,xtitle=xt,ytitle=yt,charsize=1.8
#or j=0,nlists-1 do begin
#  oplot, [wavelength[wv1:wv2]],[stacked_sed[*,j]],psym=-1.*j
#  errplot, [wavelength[wv1:wv2]],[stacked_sed[*,j]-stacked_sed_err[*,j]],$
#     [stacked_sed[*,j]+stacked_sed_err[*,j]]
for j in range(nlists):
	plt.plot(wavelength,stacked_fluxes[:,j])
plt.show()
