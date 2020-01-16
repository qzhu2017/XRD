import numpy as np
import scipy.integrate as integrate


class Similarity(object):
    """
    Class to compute the similarity between two diffraction patterns
    """

    def __init__(self,  fpeaks,  gpeaks, scaling):
        
        """
        Needs explanation
        """

        self.fpeaks = fpeaks
        self.gpeaks = gpeaks
        self.scaling = scaling
        
        try:
            self.fpeaks.shape[0] == self.gpeaks.shape[0]
        except:
            print("Patterns are not the same shape")
        
        self.N = fpeaks.shape[0]
        self.r = np.linspace(-1, 1, self.N) 
        
        self.triangleFunction()
        
    def calculate(self):
        
        gpeaks_r = self.gpeaks
        fpeaks_r = self.fpeaks
        xCorrfg = np.correlate(self.fpeaks, gpeaks_r)
        aCorrff = np.correlate(self.fpeaks, fpeaks_r)
        aCorrgg = np.correlate(self.gpeaks, gpeaks_r)

        xCorrfg_w = integrate.trapz(self.w * xCorrfg, self.r)
        aCorrff_w = integrate.trapz(self.w * aCorrff, self.r)
        aCorrgg_w = integrate.trapz(self.w * aCorrgg, self.r)

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
