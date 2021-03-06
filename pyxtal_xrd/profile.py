import numpy as np
import numba as nb

@nb.njit(nb.f8[:](nb.f8[:], nb.f8, nb.f8, nb.f8, nb.f8, nb.i8), cache = True)
def mod_pseudo_voigt(x, fwhm, A, eta_h, eta_l, N):
    
    """
    A modified split-type pseudo-Voigt function for profiling peaks
    - Izumi, F., & Ikeda, T. (2000). 
    """
    
    tmp = np.zeros((N))
    for xi, dx in enumerate(x):
        if dx < 0:
            A = A
            eta_l = eta_l
            eta_h = eta_h
        else:
            A = 1/A
            eta_l = eta_h
            eta_h = eta_l

        tmp[xi] = ((1+A)*(eta_h + np.sqrt(np.pi*np.log(2))*(1-eta_h))) /\
            (eta_l + np.sqrt(np.pi*np.log(2)) * (1-eta_l) + A*(eta_h +\
            np.sqrt(np.pi*np.log(2))*(1-eta_h))) * (eta_l*2/(np.pi*fwhm) *\
            (1+((1+A)/A)**2 * (dx/fwhm)**2)**(-1) + (1-eta_l)*np.sqrt(np.log(2)/np.pi) *\
            2/fwhm *np.exp(-np.log(2) * ((1+A)/A)**2 * (dx/fwhm)**2))
    return tmp

@nb.njit(nb.f8[:](nb.f8, nb.f8[:], nb.f8), cache = True)
def gaussian(theta2, alpha, fwhm):

    """
    Gaussian function for profiling peaks
    """

    tmp = ((alpha - theta2)/fwhm)**2
    return np.exp(-4*np.log(2)*tmp)

@nb.njit(nb.f8[:](nb.f8, nb.f8[:], nb.f8), cache = True)
def lorentzian(theta2, alpha, fwhm):

    """
    Lorentzian function for profiling peaks
    """

    tmp = 1 + 4*((alpha - theta2)/fwhm)**2
    return 1/tmp

def pseudo_voigt(theta2, alpha, fwhm, eta):

    """
    Original Pseudo-Voigt function for profiling peaks
    - Thompson, D. E. Cox & J. B. Hastings (1986). 
    """
    
    L = lorentzian(theta2, alpha, fwhm)
    G = gaussian(theta2, alpha, fwhm)
    return eta * L + (1 - eta) * G



class Profile:

    """
    This class applies a profiling function to simulated or experimentally obtained
    XRD spectra.

    Parameters
    ----------
    method: str
        Type of function used to profile
    res: float
        resolution of the profiling array in degree
    user_kwargs: dict
        The parameters for the profiling method.
    """

    def __init__(self, method='mod_pseudo-voigt', res = 0.01, user_kwargs=None):
        
        self.method = method
        self.user_kwargs = user_kwargs
        self.res = res       
        kwargs = {}

        if method == 'mod_pseudo-voigt':
           _kwargs = {
                        'U': 5.776410E-03,
                        'V': -1.673830E-03,
                        'W': 5.668770E-03,
                        'A': 1.03944,
                        'eta_h': 0.504656,
                        'eta_l': 0.611844,
                     }
        elif method == 'gaussian' or method == 'lorentzian' or method == 'pseudo-voigt':
           _kwargs = {
                        'FWHM': 0.02
                     }
        
        else:
           msg = method + " isn't supported."
           raise NotImplementedError(msg)

        kwargs.update(_kwargs)

        if user_kwargs is not None:
           kwargs.update(user_kwargs)

        self.kwargs = kwargs

    def get_profile(self, two_thetas, intensities, min2theta, max2theta):

        """
        Performs profiling with selected function, resolution, and parameters 
       
        Parameters
        ----------
        two_thetas: 1d float array 
           simulated/measured 2 theta values
        intensities: 
           simulated/measures peaks
        """
    
        N = int((max2theta-min2theta)/self.res)
        px = np.linspace(min2theta, max2theta, N) 
        py = np.zeros((N))

        for two_theta, intensity in zip(two_thetas, intensities):
            if self.method == 'gaussian':
               fwhm = self.kwargs['FWHM']
               tmp = gaussian(two_theta, px, fwhm)
            
            elif self.method == 'lorentzian':
               fwhm = self.kwargs['FWHM'] 
               tmp = lorentzian(two_theta, px, fwhm) 
            
            elif self.method == 'pseudo-voigt':
               try:
                  fwhm_g = self.kwargs['FWHM-G'] 
                  fwhm_l = self.kwargs['FWHM-L']
               except:
                  fwhm_g = self.kwargs['FWHM']
                  fwhm_l = self.kwargs['FWHM'] 
                
               fwhm = (fwhm_g**5 + 2.69269*fwhm_g**4*fwhm_l + 2.42843*fwhm_g**3*fwhm_l**2 +
                       4.47163*fwhm_g**2*fwhm_l**3 + 0.07842*fwhm_g*fwhm_l**4 + fwhm_l**5)**(1/5)
               
               eta = 1.36603*fwhm_l/fwhm - 0.47719*(fwhm_l/fwhm)**2 + 0.11116*(fwhm_l/fwhm)**3

               tmp = pseudo_voigt(two_theta, px, fwhm, eta)

            elif self.method == 'mod_pseudo-voigt':
               U = self.kwargs['U']
               V = self.kwargs['V']
               W = self.kwargs['W']
               A = self.kwargs['A']
               eta_h = self.kwargs['eta_h']
               eta_l = self.kwargs['eta_l']
               
               fwhm = np.sqrt(U*np.tan(np.pi*two_theta/2/180)**2 + V*np.tan(np.pi*two_theta/2/180) + W)
               x = px - two_theta
               tmp = mod_pseudo_voigt(x, fwhm, A, eta_h, eta_l, N)
            
            py += intensity * tmp

        py /= np.max(py)

        self.spectra = np.vstack((px,py))
        return self.spectra
