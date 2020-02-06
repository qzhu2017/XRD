#!/usr/bin/env  python
# encoding: utf-8

"""
This file to handle various profile functions
"""

import numpy as np

def gaussian_profile(I0, theta2, alpha, fwhm):
    tmp = ((alpha - theta2)/fwhm)**2
    return I0 * np.exp(-4*np.log(2)*tmp)

def lorentzian_profile(I0, theta2, alpha,fwhm):
    tmp = 1 + 4*((alpha - theta2)/fwhm)**2
    return I0 * 1/tmp

def split_type(x, I0, N, fwhm, A, eta_h, eta_l):
    tmp = np.zeros((N))
    for k, dx in enumerate(x):
        if dx >= 0:
            A = 1/A

        tmp[k] = ((1+A)*(eta_h + np.sqrt(np.pi*np.log(2))*(1-eta_h))) /\
                 (eta_l + np.sqrt(np.pi*np.log(2)) * (1-eta_l) + A*(eta_h + \
                 np.sqrt(np.pi*np.log(2))*(1-eta_h))) * (eta_l*2/(np.pi*fwhm) * \
                 (1+((1+A)/A)**2 * (dx/fwhm)**2)**(-1) + (1-eta_l)*np.sqrt(np.log(2)/np.pi)* \
                 2/fwhm *np.exp(-np.log(2) * ((1+A)/A)**2 * (dx/fwhm)**2))

    return I0 * tmp


class profile():

    """

    Parameters
    ----------
    method: str
        Type of minimization scheme, e.g.: 'LBFGS'.
    user_kwargs: dict
        The arguments for the optimization method. These arguments are 
        passed with a dict.
    """

    def __init__(self, method, user_kwargs):
        self.method = method
        self.kwargs = {'method': method}    
        if method == 'split-type':
            # from ****
            _kwargs = {
                        'U': 5.776410E-03 # FWHM parameter, U
                        'V': -1.673830E-03 # FWHM parameter, V
                        'W': 5.668770E-03 # FWHM parameter, W
                        'A': 1.03944 # Asymmetry parameter, a1
                        'eta_h': 0.504656 # Mixing parameter, eta_H0
                        'eta_l': 0.611844  # Mixing parameter, eta_L0
                      }
        elif method == 'gaussian':
            _kwargs = {'V': }

        elif method == 'lorentzian':
            _kwargs = {'X': }

        else:
            msg = f"The {method} is not implemented yet."
            raise NotImplementedError(msg)

        self.kwargs.update(_kwargs)

        if user_kwargs is not None:
            self.kwargs.update(user_kwargs)




