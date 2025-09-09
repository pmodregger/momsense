import numpy as np
import pylab as plt
import pyFAI
import fabio
from scipy.signal import find_peaks



#%% load data
frame = fabio.open('ID11_martensite_averaged_frame.tif').data
mask = fabio.open('eiger_mask_E-08-0144_20240205.edf').data
poni = pyFAI.load('CeO2_20240920_1158_office.poni')

#%% find peak positions
twotheta_grad, curve = poni.integrate1d_ng(frame,1000,
                                 unit='2th_deg',
                                 correctSolidAngle=True,
                                 normalization_factor = 1,
                                 method=('full','csc','cython'),
                                 radial_range=(4,15),
                                 polarization_factor=None,mask=mask)

peakinds = find_peaks(curve,height=0.01*curve.max(),distance=80)[0]
peakpos = twotheta_grad[peakinds]
peakint = curve[peakinds]

plt.semilogy(twotheta_grad,curve/curve.max())
plt.semilogy(peakpos,peakint/curve.max(),'*')
plt.grid()
plt.show()



       
#%% call momsense
Ntheta = 200
twotheta_range = 2. # in deg

import momsense

Nphotons, M0, uM0, M1, uM1, M2, ueps = momsense.momsense(frame,poni,peakpos,
                                        azimuthal_range=(-22.5,22.5),twotheta_range=0.5,
                                        mask=mask,Ntheta=200,verbose=1)



