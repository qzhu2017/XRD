import numpy as np
import numba as nb
import scipy.integrate as integrate
import matplotlib.pyplot as plt
from scipy import interpolate

class Similarity(object):
    """
    Class to compute the similarity between two diffraction patterns
    """

    def __init__(self, f, g, N = None, x_range = None, l = None, weight = 'cosine'):
        
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
        if N is None:
            self.N = max([len(self.fx), len(self.gx)])
        else:
            self.N = N
        self.x_range = x_range
        self.preprocess()
        if l is None:
            self.l = (self.x_range[1] - self.x_range[0])/10
        else:
            self.l = abs(l)
        self.weight = weight
        self.r = np.linspace(-self.l, self.l, self.N)
        if self.weight == 'triangle':
            self.triangleFunction()
        elif self.weight == 'cosine':
            self.cosineFunction()
        else:
            msg = function + 'is not supported'
            raise NotImplementedError(msg)


    def calculate(self):
        self.S = self._calculate(self.r,self.w,self.d,self.Npts,self.fy,self.gy)
        return self.S
    
    @staticmethod
    @nb.njit(nb.f8(nb.f8[:], nb.f8[:], nb.f8, nb.i8, nb.f8[:], nb.f8[:]),nopython=True)
    def _calculate(r,w,d,Npts,fy,gy):
        
        """
        Compute the similarity between the pair of spectra f, g
        """
        
        xCorrfg_w = 0
        aCorrff_w = 0
        aCorrgg_w = 0
        count = 0

        for r0, w0 in zip(r, w):
            Corrfg, Corrff, Corrgg = 0, 0, 0
            for i in range(Npts):
                shift = int(round(r0/d))
                if 0 <= i + shift <= Npts-1:
                    Corrfg += fy[i]*gy[i+shift]
                    Corrff += fy[i]*fy[i+shift]
                    Corrgg += gy[i]*gy[i+shift]

            xCorrfg_w += w0*Corrfg 
            aCorrff_w += w0*Corrff
            aCorrgg_w += w0*Corrgg

        return np.abs(xCorrfg_w / np.sqrt(aCorrff_w * aCorrgg_w))


    def preprocess(self):

        """
        Preprocess the input spectra f and g
        """
        if self.x_range == None:
            x_min = max(np.min(self.fx), np.min(self.gx))
            x_max = min(np.max(self.fx), np.max(self.gx))
            self.x_range = [x_min,x_max]
        f_inter = interpolate.interp1d(self.fx, self.fy, 'cubic', fill_value = 'extrapolate')
        g_inter = interpolate.interp1d(self.gx, self.gy, 'cubic', fill_value = 'extrapolate')

        fgx_new = np.linspace(self.x_range[0], self.x_range[1], self.N)
        fy_new = f_inter(fgx_new)
        gy_new = g_inter(fgx_new)

        self.fx, self.fy, self.gy = fgx_new, fy_new, gy_new
        self.Npts = len(self.fx)
        self.d = (self.fx[-1] - self.fx[0])/self.Npts

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
        plt.title("{:6f}".format(self.S))
        plt.show()

