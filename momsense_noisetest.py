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
    
       
#%% repeat over noise realisations
Ntheta = 200
twotheta_range = 2. # in deg
ispeak = 0

radial_range = (peakpos[ispeak]-twotheta_range,peakpos[ispeak]+twotheta_range)
Dtheta = np.diff(radial_range)[0]/Ntheta
normfactor = angular_pixel_size**2 / np.sin(np.deg2rad(peakpos[ispeak]))/np.deg2rad((azimuthal_range[1]-azimuthal_range[0]))/np.deg2rad(Dtheta)



L = 100
Nphotons = np.zeros((L))
M0 = np.zeros((L))
M1 = np.zeros((L))
M2 = np.zeros((L))
uM0 = np.zeros((L))
uM1 = np.zeros((L))
ueps = np.zeros((L))

for isl in range(L):
    
    # add photon shot noise to frame
    noise_frame = np.random.poisson(lam=frame).astype('float')
    # azimuthal integration    
    twotheta_grad_peak, curve_peak = poni_data.integrate1d_ng(noise_frame,Ntheta,
                                                 unit='2th_deg',azimuth_range=azimuthal_range,
                                                 correctSolidAngle=True,
                                                 normalization_factor = normfactor,
                                                 radial_range=radial_range,
                                                 method=('full','csc','cython'),
                                                 polarization_factor=None)
    # plt.plot(twotheta_grad_peak,curve_peak)
    # plt.show()
    
    Nphotons[isl] = np.sum(curve_peak)
    Dx = np.diff(twotheta_grad_peak[:2])[0]
    M0[isl] = np.sum(curve_peak)*Dx
    M1[isl] = np.sum(curve_peak*twotheta_grad_peak)*Dx/M0[isl]
    M2[isl] = np.sum(curve_peak*(twotheta_grad_peak-M1[isl])**2)*Dx/M0[isl]


#% results
uM0_teo = np.sqrt(Dx*np.mean(M0))
uM0_num = np.std(M0)

uM1_teo = np.deg2rad(np.sqrt(Dx*np.mean(M2)/np.mean(M0))) # in 2theta in deg
uM1_num = np.deg2rad(np.std(M1))

ueps_teo = uM1_teo/2/np.tan(np.deg2rad(peakpos[ispeak]))
ueps_num = uM1_num/2/np.tan(np.deg2rad(peakpos[ispeak]))

print('uM0 teo: ',uM0_teo)
print('uM0 num: ',uM0_num)

print('uM1 teo in theta and rad: ',uM1_teo)
print('uM1 num in theta and rad: ',uM1_num)

print('ueps teo: ',ueps_teo)
print('ueps num: ',ueps_num)
