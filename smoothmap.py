import pdb
import numpy as np
from shift import shift_twod

def smoothmap(mapin, psfin):

	s = np.shape(mapin)
	mnx = s[0]
	mny = s[1]

	s = np.shape(psfin)
	pnx = s[0]
	pny = s[1]

	psf_x0 = pnx/2
	psf_y0 = pny/2
	psf = psfin
	px0 = psf_x0
	py0 = psf_y0
	# pad psf
	psfpad = np.zeros([mnx, mny])
	psfpad[0:pnx,0:pny] = psf

	# shift psf so that centre is at (0,0)
	psfpad = shift_twod(psfpad, -px0, -py0)

	smmap = np.real( np.fft.ifft2( np.fft.fft2(mapin) *
    	np.fft.fft2(psfpad)) 
		) # * mnx * mny

	return smmap