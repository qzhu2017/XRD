import numpy as np
import scipy.integrate as integrate
from scipy import interpolate

class Similarity(object):
    """
    Class to compute the similarity between two diffraction patterns
    """

    def __init__(self, f, g, N = None, x_range = None, r_range = None, weight = {'function': 'triangle', 'params': 0.6}):
        
        """
        Args:

        f: spectra1 (2D array)
        g: spectra2 (2D array)
        x_range: the range of x values used to compute similarity ([x_min, x_max])
        N: number of sampling points for the processed spectra
        weight: weight function used to compute the similarity (dictionary)
        """

        self.fx, self.fy = f[0], f[1]
        self.gx, self.gy = g[0], g[1]
        self.N = N
        self.x_range = x_range
        if r_range == None:
            self.r_range = [-1,1]
        else:
            self.r_range = r_range
        self.weight = weight
        
        self.r = np.linspace(self.r_range[0], self.r_range[1], self.N)
        self.preprocess()
        function = self.weight['function']
        if function == 'triangle':
            self.triangleFunction()
        else:
            msg = function + 'is not supported'
            raise NotImplementedError(msg)

    def calculate(self):
        
        """
        Compute the similarity between the pair of spectra f, g
        """

        fx_r = self.fx + self.r
        fy_r = self.fy
        gx_r = self.gx + self.r
        gy_r = self.gy

        fg_dx = np.linspace(self.fx[0], gx_r[-1], self.N)
        ff_dx = np.linspace(self.fx[0], fx_r[-1], self.N)
        gg_dx = np.linspace(self.gx[0], gx_r[-1], self.N)

        xCorrfg = integrate.trapz(self.fy*gy_r, fg_dx)
        aCorrff = integrate.trapz(self.fy*fy_r, ff_dx)
        aCorrgg = integrate.trapz(self.gy*gy_r, gg_dx)

        xCorrfg_w = integrate.trapz(self.w*xCorrfg, self.r)
        aCorrff_w = integrate.trapz(self.w*aCorrff, self.r)
        aCorrgg_w = integrate.trapz(self.w*aCorrgg, self.r)
    
        return xCorrfg_w / np.sqrt(aCorrff_w * aCorrgg_w)

    def preprocess(self):

        """
        Preprocess the input spectra f and g
        """

        if self.x_range == None:
            x_min = min(np.min(self.fx), np.min(self.gx))
            x_max = max(np.max(self.fx), np.max(self.gx))
            self.x_range = [x_min,x_max]

        f_inter = interpolate.interp1d(self.fx, self.fy, 'cubic', fill_value = 'extrapolate')
        g_inter = interpolate.interp1d(self.gx, self.gy, 'cubic', fill_value = 'extrapolate')
        fgx_new = np.linspace(self.x_range[0], self.x_range[1], self.N)
        fy_new = f_inter(fgx_new)
        gy_new = g_inter(fgx_new)

        self.fx, self.fy = fgx_new, fy_new
        self.gx, self.gy = fgx_new, gy_new

    def triangleFunction(self):
        
        """
        Function to weight correlations
        """
        
        w = np.zeros((self.N))
        l = self.weight['params']
        for i in range(self.r.shape[0]):
            r = np.abs(self.r[i])
            if r < 1:
                tf = lambda r,l : 1 - r/l
                w[i] = tf(r,l)
            else:
                w[i] = 0
        self.w = w
