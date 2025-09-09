# momsense
Rapid estimation of angular moment sensitivity in x-ray diffraction from a single diffraction pattern

# Reference
Users are kindly asked to cite the following reference:

"Ultimate Sensitivity in X-ray Diffraction: Angular Moments vs. Shot Noise"
Peter Modregger, Felix Wittwer, Ahmar Khaliq, Niklas Pyrlik, James A. D. Ball, Jan Garrevoet, Gerald Falkenberg, Alexander Liehr, Michael Stuckelberger,
Journal of Applied Crystallography 58 (2025)
doi.org/10.1107/S1600576725006715

# Example usage
import momsense
Nphotons, M0, uM0, M1, uM1, M2, ueps = momsense.momsense(frame,poni,peakpos,azimuthal_range=(-180,180),twotheta_range=1.,Ntheta=200,verbose=1)

# Application examples included in
momsense_num_example.py - numerical example
momsense_exp_example.py - example from ID11 (ESRF)

# Dependencies
main library uses pyFAI
additional usage of find_peaks
