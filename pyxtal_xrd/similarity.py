import numpy as np
import scipy.integrate as integrate
import matplotlib.pyplot as plt
from scipy import interpolate

class Similarity(object):
    """
    Class to compute the similarity between two diffraction patterns
    """

    def __init__(self, f, g, N = None, x_range = [-1,1], l = 1, weight = 'cosine'):
        
        """
        Args:

        f: spectra1 (2D array)
        g: spectra2 (2D array)
        N: number of sampling points for the processed spectra
        x_range: the range of x values used to compute similarity ([x_min, x_max])
        l: cutoff value for shift (real)
        weight: weight function 'triangle' or 'cosine' (str)
        """
        self.fx, self.fy = f[0], f[1]
        self.gx, self.gy = g[0], g[1]
        self.N = N
        self.x_range = x_range
        self.l = abs(l)
        self.weight = weight
        self.r = np.linspace(-self.l, self.l, self.N)
        self.preprocess()
        if self.weight == 'triangle':
            self.triangleFunction()
        elif self.weight == 'cosine':
            self.cosineFunction()
        else:
            msg = function + 'is not supported'
            raise NotImplementedError(msg)

        self.showPlot()
    
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

        return np.abs(xCorrfg_w / np.sqrt(aCorrff_w * aCorrgg_w))

    def preprocess(self):

        """
        Preprocess the input spectra f and g
        """
        if self.x_range == None:
            x_min = max(np.min(self.fx), np.min(self.gx))
            x_max = min(np.max(self.fx), np.max(self.gx))
            self.x_range = [x_min,x_max]
        f_inter = interpolate.UnivariateSpline(self.fx, self.fy)#, 'cubic', fill_value = 'extrapolate')
        g_inter = interpolate.UnivariateSpline(self.gx, self.gy)#, 'cubic', fill_value = 'extrapolate')

        fgx_new = np.linspace(self.x_range[0], self.x_range[1], self.N)
        fy_new = f_inter(fgx_new)
        gy_new = g_inter(fgx_new)

        self.fx, self.fy = fgx_new, fy_new
        self.gx, self.gy = fgx_new, gy_new
        
    def triangleFunction(self):
        
        """
        Triangle function to weight correlations
        """
        
        w = np.zeros((self.N))
        l = self.l
        for i in range(self.r.shape[0]):
            r = np.abs(self.r[i])
            if r <= l:
                tf = lambda r,l : 1 - r/l
                w[i] = tf(r,l)
            else:
                w[i] = 0
        self.w = w

    def cosineFunction(self):

        """
        cosine function to weight correlations
        """
        
        w = np.zeros((self.N))
        l = self.l
        for i in range(self.r.shape[0]):
            r = np.abs(self.r[i])
            if r <= l:
                tf = lambda r,l : 0.5 * (np.cos(np.pi * r/l) + 1)
                w[i] = tf(r,l)
            else:
                w[i] = 0
        self.w = w

    def showPlot(self):
        fig1 = plt.figure(1,figsize=(15,6))
        frame1=fig1.add_axes((.1,.3,.8,.6))
    
        plt.plot(self.fx,self.fy,label='pxrd')
        plt.plot(self.gx,-self.gy,label='vesta')
        plt.legend()
        #Residual plot
        residuals = self.gy-self.fy
        frame2=fig1.add_axes((.1,.1,.8,.2))        
        plt.plot(self.gx,residuals,'.r', markersize = 0.5)
        plt.show()

