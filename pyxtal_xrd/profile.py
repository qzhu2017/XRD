import numpy as np
import numba as nb


@nb.njit(nb.f8[:](nb.f8[:], nb.f8, nb.f8, nb.f8, nb.f8, nb.i8), cache = True)
def pseudo_voigt(x, fwhm, A, eta_h, eta_l, N):
    
    """
    A modified split-type pseudo-Voigt function for profiling peaks
    (Izumi, F., & Ikeda, T. (2000). 
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

def gaussian(theta2, alpha, fwhm):

    """
    Gaussian function for profiling peaks
    """

    tmp = ((alpha - theta2)/fwhm)**2
    return np.exp(-4*np.log(2)*tmp)

def lorentzian(theta2, alpha, fwhm):

    """
    Lorentzian function for profiling peaks
    """

    tmp = 1 + 4*((alpha - theta2)/fwhm)**2
    return 1/tmp



class Profile:

    """
    This class applies a profiling function to simulated or experimentally obtained
    XRD spectra.

    Parameters
    ----------
    method: str
        Type of function used to profile
    res: float
        resolution of the profiling array
    user_kwargs: dict
        The parameters for the profiling method.
    """

    def __init__(self, method='pseudo_voigt', res = 0.01, user_kwargs=None):
        
        self.method = method
        self.N = int(1/res)
        self.user_kwargs = user_kwargs
        
        kwargs = {}

        if method == 'pseudo_voigt':
           _kwargs = {
                        'U': 5.776410E-03,
                        'V': -1.673830E-03,
                        'W': 5.668770E-03,
                        'A': 1.03944,
                        'eta_h': 0.504656,
                        'eta_l': 0.611844,
                     }
        elif method == 'gaussian' or method == 'lorentzian':
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

    def get_profile(self, two_thetas, intensities):

       """
       Performs profiling with selected function, resolution, and parameters 
       
       Parameters
       ----------
       two_thetas: 1d float array 
           simulated/measured 2 theta values
       intensities: 
           simulated/measures peaks
       """

       px = np.linspace(np.min(two_thetas) - 5, np.max(two_thetas) + 5, self.N) 
       py = np.zeros((self.N))

       for two_theta, intensity in zip(two_thetas, intensities):
           if self.method == 'gaussian':
              fwhm = self.kwargs['FWHM']
              tmp = gaussian(two_theta, px, fwhm)
           
           elif self.method == 'lorentzian':
              fwhm = self.kwargs['FWHM'] 
              tmp = lorentzian(two_theta, px, fwhm) 

           elif self.method == 'pseudo_voigt':
              U = self.kwargs['U']
              V = self.kwargs['V']
              W = self.kwargs['W']
              A = self.kwargs['A']
              eta_h = self.kwargs['eta_h']
              eta_l = self.kwargs['eta_l']
              
              fwhm = np.sqrt(U*np.tan(np.pi*two_theta/2/180)**2 + V*np.tan(np.pi*two_theta/2/180) + W)
              x = px - two_theta
              tmp = pseudo_voigt(x, fwhm, A, eta_h, eta_l, self.N)
           
           py += intensity * tmp

       py /= np.max(py)

       self.spectra = np.vstack((px,py))
