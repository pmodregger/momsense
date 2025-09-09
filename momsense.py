import numpy as np
import pyFAI
import pylab as plt

def momsense(frame,poni,peakpos,azimuthal_range=(-180,180),twotheta_range=1.,Ntheta=200,
             verbose=1,mask=None):
    
    """
    
    Takes 1 diffraction pattern and estimates the M0, M1 sensitivity for all peaks given in peakpos
    assuming photon shot noise
    
    USAGE:
        Nphotons, M0, uM0, M1, uM1, M2, ueps = momsense(frame,poni,peakpos,azimuthal_range=(-180,180),
                                                        twotheta_range=1.,Ntheta=200,verbose=1)
        
    INPUT:
        frame  - 2D numpy array of the diffraction frame
        poni   - poni structure for pyFAI
        peakpos - list of 2theta peak positions to be analysed
        
    
    """
    
    angular_pixel_size = np.arctan(poni.pixel1/poni.dist)
    Npeaks = len(peakpos)
    Nphotons = np.zeros(Npeaks)
    M0 = np.zeros(Npeaks)
    M1 = np.zeros(Npeaks)
    M2 = np.zeros(Npeaks)
    uM0 = np.zeros(Npeaks)
    uM1 = np.zeros(Npeaks)
    ueps = np.zeros(Npeaks)
    
    for ispeak in range(len(peakpos)):
    
        # determine the correct normalization factor
        radial_range = (peakpos[ispeak]-twotheta_range,peakpos[ispeak]+twotheta_range)
        Dtheta = np.diff(radial_range)[0]/Ntheta
        normfactor = angular_pixel_size**2 / np.sin(np.deg2rad(peakpos[ispeak]))/np.deg2rad((azimuthal_range[1]-azimuthal_range[0]))/np.deg2rad(Dtheta)

        # azimuthal integration
        twotheta_grad, curve = poni.integrate1d_ng(frame,Ntheta,
                                 unit='2th_deg',azimuth_range=azimuthal_range,
                                 correctSolidAngle=True,
                                 normalization_factor = normfactor,
                                 radial_range=radial_range,method=('full','csc','cython'),
                                 polarization_factor=None,mask=mask)

        if verbose >= 2:
            plt.plot(twotheta_grad,curve)
            plt.show()

        # determine moments
        Nphotons[ispeak] = np.sum(curve)
        Dx = np.diff(twotheta_grad[:2])[0]
        M0[ispeak] = np.sum(curve)*Dx
        M1[ispeak] = np.sum(curve*twotheta_grad)*Dx/M0[ispeak]
        M2[ispeak] = np.sum(curve*(twotheta_grad-M1[ispeak])**2)*Dx/M0[ispeak]

        # determine sensitivities    
        uM0[ispeak] = np.sqrt(Dx*np.mean(M0[ispeak]))
        uM1[ispeak] = np.deg2rad(np.sqrt(Dx*M2[ispeak]/M0[ispeak])) # in 2theta in rad
        ueps[ispeak] = uM1[ispeak]/2/np.tan(np.deg2rad(peakpos[ispeak]/2))
        
        # show results if desired
        if verbose:
            print('\n')
            print('Peak: ',ispeak)
            print('peakpos / deg: ',peakpos[ispeak])
            print('Nphotons: ',Nphotons[ispeak])
            print('M0: ',M0[ispeak])
            print('u(M0): ',uM0[ispeak])
            print('M1 / 2theta & rad: ',np.deg2rad(M1[ispeak]))
            print('u(M1) / 2theta & rad: ',uM1[ispeak])
            print('M2 / 2theta & deg: ',M2[ispeak])
            print('u(eps): ',ueps[ispeak])
            
            
    return Nphotons, M0, uM0, M1, uM1, M2, ueps
            
    