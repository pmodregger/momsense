import numpy as np
import pylab as plt
import pyFAI


#%% poni file 
poni_file = 'test.poni'
poni_data = pyFAI.load(poni_file)

poni_dist = 0.24

#%% parameters
angular_pixel_size = np.arctan(poni_data.pixel1/poni_data.dist)

peakpos = np.array((5.,7.5,10.,12.5,15.,17.5,20.)) # 2theta in deg
peakwidth = 5*angular_pixel_size # in rad

azimuthal_range = np.array((-180,180))

#%% generate test frame
x = (np.arange(2048)-1024)*angular_pixel_size
XX, YY = np.meshgrid(x,x)
R = np.sqrt(XX**2+YY**2)

frame = np.zeros_like(R) 
sumpeak_org = np.zeros(len(peakpos))
for ispeak in range(len(peakpos)):

    diffring = 100*np.exp(-0.5*( (R-np.tan(np.deg2rad(peakpos[ispeak]))   )/peakwidth)**2  ) * np.exp(-(peakpos[ispeak]/peakpos[0]))
    sumpeak_org[ispeak] = np.sum(diffring)
    
    frame += diffring
    
       
#%% call momsense
Ntheta = 200
twotheta_range = 2. # in deg

import momsense

Nphotons, M0, uM0, M1, uM1, M2, ueps = momsense.momsense(frame,poni_data,peakpos,azimuthal_range=(-180,180),twotheta_range=1.,Ntheta=200,verbose=1)



