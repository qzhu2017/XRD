import numpy as np
import scipy.integrate as integrate


class Similarity(object):
    """
    Class to compute the similarity between two diffraction patterns
    """

    def __init__(self,  fpeaks,f2thetas,  gpeaks,g2thetas, scaling):
        
        """
        Needs explanation
        """

        self.fpeaks = fpeaks
        self.f2thetas = f2thetas
        self.gpeaks = gpeaks
        self.g2thetas = g2thetas
        self.scaling = scaling
        
        try:
            self.fpeaks.shape[0] == self.gpeaks.shape[0]
        except:
            print("Patterns are not the same shape")
        
        self.N = fpeaks.shape[0]
        self.r = np.linspace(-1, 1, self.N) 

        self.triangleFunction()
        
    def calculate(self):
        fpeaks_r = self.fpeaks
        f2thetas_r = self.f2thetas + self.r
        gpeaks_r = self.gpeaks
        g2thetas_r = self.g2thetas + self.r

        x_fg = np.concatenate((self.f2thetas, g2thetas_r))
        x_fg = np.linspace(x_fg[0],x_fg[-1],self.N)
        x_ff = np.concatenate((self.f2thetas, f2thetas_r))
        x_ff = np.linspace(x_ff[0],x_ff[-1],self.N)
        x_gg = np.concatenate((self.g2thetas, g2thetas_r))
        x_gg = np.linspace(x_gg[0],x_gg[-1],self.N)

        xCorrfg = integrate.trapz(self.fpeaks*gpeaks_r,x_fg)
        aCorrff = integrate.trapz(self.fpeaks*fpeaks_r,x_ff)
        aCorrgg = integrate.trapz(self.gpeaks*gpeaks_r,x_gg)

        xCorrfg_w = integrate.trapz(self.w*xCorrfg, self.r)
        aCorrff_w = integrate.trapz(self.w*aCorrff, self.r)
        aCorrgg_w = integrate.trapz(self.w*aCorrgg, self.r)
    
        return xCorrfg_w / np.sqrt(aCorrff_w * aCorrgg_w)

    def triangleFunction(self):
        w = np.zeros((self.N))
        l = self.scaling
        for i in range(self.r.shape[0]):
            r = np.abs(self.r[i])
            if r < 1:
                tf = lambda r,l : 1 - r/l
                w[i] = tf(r,l)
            else:
                w[i] = 0
        self.w = w
