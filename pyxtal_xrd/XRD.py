import numpy as np
import matplotlib.pyplot as plt
import json
import os
import collections
from ase.io import read
from database.element import Element

def create_index():
    hkl_index = []
    for i in [-1,0,1]:
        for j in [-1,0,1]:
            for k in [-1,0,1]:
                hkl = np.array([i,j,k])
                if sum(hkl*hkl)>0:
                    hkl_index.append(hkl)
    return np.array(hkl_index)

class XRD(object):
    """a XRD class
    """

    def __init__(self, crystal, wavelength=1.54184, max2theta=180, 
                 preferred_orientation = False, march_parameter = None):
       
        """ Class to compute the powder XRD.
        
        Parameters
        ----------
        crystal: ase atoms object
        wavelength: float 
        max2theta: float
        preferred_orientation: boolean
        march_parameter: float
        """
        self.profiling = None
        self.wavelength = wavelength
        self.max2theta = np.radians(max2theta)
        self.name = crystal.get_chemical_formula()
        self.preferred_orientation = preferred_orientation
        self.march_parameter = march_parameter
        self.all_dhkl(crystal)
        self.intensity(crystal)
        self.pxrdf()     
    
        
    def by_hkl(self, hkl):
        """ d for any give abitray [h,k,l] index """
        # this is a simple print statement, does not need to be optimized

        id1 = np.where(np.all(self.hkl_list == np.array(hkl), axis=1 ))
        if id1 is None:
           print('This hkl is not in the given 2theta range')
        else:
           print('  2theta     d_hkl     hkl       Intensity')
           for i in id1[0]:
                print('%8.3f  %8.3f   [%2d %2d %2d] %8.2f' % \
                (np.degrees(self.theta2[i]), self.d_hkl[i], \
                 self.hkl_list[i,0], self.hkl_list[i,1], self.hkl_list[i,2], \
                 self.xrd_intensity[i] ))
    
    def all_dhkl(self, crystal):
        """ 3x3 representation -> 1x6 (a, b, c, alpha, beta, gamma)"""
        rec_matrix = crystal.get_reciprocal_cell()
        d_min = self.wavelength/np.sin(self.max2theta/2)/2

        # This block is to find the shortest d_hkl, 
        # for all basic directions (1,0,0), (0,1,0), (1,1,0), (1,-1,0) 
        hkl_index = create_index()
        multiple = np.round(1/np.linalg.norm(np.dot(hkl_index, rec_matrix), axis=1)/d_min)
        hkl_max = np.max(np.abs(np.einsum('i,ij->ij', multiple, hkl_index)), axis=1)
        h1, k1, l1 = int(hkl_max[0]), int(hkl_max[1]), int(hkl_max[2])

        h = np.arange(-h1,h1+1)
        k = np.arange(-k1,k1+1)
        l = np.arange(-l1,l1+1)

        hkl = np.array((np.meshgrid(h,k,l))).transpose()
        hkl_list = np.reshape(hkl, [len(h)*len(k)*len(l),3])
        hkl_list = hkl_list[np.where(hkl_list.any(axis=1))[0]]
        d_hkl = 1/np.linalg.norm(np.dot(hkl_list, rec_matrix), axis=1)

        shortlist = d_hkl > (d_min)
        d_hkl = d_hkl[shortlist]
        hkl_list = hkl_list[shortlist]
        sintheta = self.wavelength/2/d_hkl

        self.theta = np.arcsin(sintheta)
        self.hkl_list = hkl_list
        self.d_hkl = d_hkl

    def intensity(self, crystal):

        """
        This function calculates all that is necessary to find the intensities.
        This scheme is based off of pymatgen
        Needs improvement from different correction factors.
        """
        # open a json file with atomic scattering parameters, should eventuall go to Element class

        with open(os.path.join(os.path.dirname(__file__),
                       "atomic_scattering_params.json")) as f:
                        ATOMIC_SCATTERING_PARAMS = json.load(f)

        d0 = (1/2/self.d_hkl)**2

        # obtiain scattering parameters, atomic numbers, and occus (need to look into occus)
        coeffs = []
        zs = []

        for elem in crystal.get_chemical_symbols():
            if elem == 'D':
                elem = 'H'
            c = ATOMIC_SCATTERING_PARAMS[elem]
            z = Element(elem).z
            coeffs.append(c)
            zs.append(z) 

        coeffs = np.array(coeffs)
        self.peaks = {}
        two_thetas = []

        # self.march_parameter = 1
        TWO_THETA_TOL = 1e-5 # tolerance to find repeating angles
        SCALED_INTENSITY_TOL = 1e-5 # threshold for intensities
        
        ind = 0
        intense = []
        angle = []
        count = 0

        for hkl, s2, theta, d_hkl in zip(self.hkl_list, d0, self.theta, self.d_hkl):
            count+=1
            
            # calculate the scattering factor sf
            g_dot_r = np.dot(crystal.get_scaled_positions(), np.transpose([hkl])).T[0]
            sf = zs - 41.78214 * s2 * np.sum(coeffs[:, :, 0] * np.exp(-coeffs[:, :, 1] * s2), axis=1)
            
            # calculate the structure factor f
            f = np.sum(sf * np.exp(2j * np.pi * g_dot_r))
            
            # calculate the lorentz polarization factor lf
            lf = (1 + np.cos(2 * theta) ** 2) / (np.sin(theta) ** 2 * np.cos(theta))

            # calculate the preferred orientation factor
            if self.preferred_orientation != False:
                G = self.march_parameter
                po = ((G * np.cos(theta))**2 + 1/G * np.sin(theta)**2)**(-3/2) 
            else:
                po = 1
    
            # calculate the intensity I
            I = (f * f.conjugate()).real
            
            # calculate 2*theta
            two_theta = np.degrees(2 * theta)
            
            # find where the scattered angles are equal
            ind = np.where(np.abs(np.subtract(two_thetas, two_theta)) < TWO_THETA_TOL)

            # append intensity, hkl plane, and thetas to lists
            if len(ind[0]) > 0:
                 self.peaks[two_thetas[ind[0][0]]][0] += I * lf * po
                 self.peaks[two_thetas[ind[0][0]]][1].append(tuple(hkl))
            else:
                self.peaks[two_theta] = [I * lf * po, [tuple(hkl)],d_hkl]
                two_thetas.append(two_theta)

        # obtain important intensities (defined by SCALED_INTENSITY_TOL)
        # and corresponding 2*theta, hkl plane + multiplicity, and d_hkl
        # print(peaks.keys())
        max_intensity = max([v[0] for v in self.peaks.values()])
        x = []
        y = []
        hkls = []
        d_hkls = []
        count = 0 
        for k in sorted(self.peaks.keys()):
            count +=1
            v = self.peaks[k]
            fam = self.get_unique_families(v[1])
            if v[0] / max_intensity * 100 > SCALED_INTENSITY_TOL:
                x.append(k)
                y.append(v[0])
                
                hkls.append([{"hkl": hkl, "multiplicity": mult}
                             for hkl, mult in fam.items()])
                d_hkls.append(v[2])

        self.theta2 = x
        self.xrd_intensity = y
        self.hkl_list = hkls
        self.d_hkl = d_hkls

 
    def get_profile(self, theta2, xrd_intensity, N, **kwargs):
    
        """
        args:
            theta2 (1D array): simulated theta values 
            xrd_intensity (1D array): simulated peaks 
            N (int): Resolution for profiling arrays 
            Keyword Arguments (dict): The user can choose between theta dependent OR independent fwhm.
                Key 1 (string): Function used to profile spectra 
                                (gaussian, lorentzian, or split-type).
                Theta independent:
                    Key 2 (float): constant value for fwhm
                Theta dependent:
                    Key 1 = gaussian:
                        Key 2 (float/int): adjustable parameter (V) for theta dep. fwhm
                    Key 1 = lorentzian:
                        key 2 (float/int): adjustable parameter (X) for theta dep. fwhm
                    Key 1 = split-type:
                        Key 2-4 (float/int): adjustable parameters (U, V, W) for theta dep. fwhm
                        Key 5 (float/int): assymetry parameter (A)
                        key 6 (float/int): mixing parameter (eta_l)
                        key 7 (float/int): mixing parameter (eta_h)
        stores:
            self.spectra: x and y values of profiling function (2D array)
        """

        # profile parameters
        tail = 5
        profile = kwargs['function']
        gpeaks = np.zeros((N))
        g2thetas = np.linspace(np.min(theta2) - tail, np.max(theta2) + tail, N)
        for i,j in zip(range(len(theta2)),range(len(xrd_intensity))):
            if profile == 'gaussian':
                try:
                    V = kwargs['V']
                    fwhm = V*np.tan(np.pi*theta2[i]/2/180)
                except:
                    fwhm = kwargs['FWHM']
                tmp = self.gaussian_profile(xrd_intensity[i],theta2[i],g2thetas,fwhm)

            elif profile == 'lorentzian':
                try:
                    X = kwargs['X']
                    fwhm = X/np.cos(np.pi*theta2[i]/2/180)
                except:
                    fwhm = kwargs['FWHM']
                tmp = self.lorentzian_profile(xrd_intensity[i],theta2[i],g2thetas,fwhm)
            
            elif profile == 'split-type': 
                try:
                    U = kwargs['U']
                    V = kwargs['V']
                    W = kwargs['W']
                    fwhm = np.sqrt(U*np.tan(np.pi*theta2[i]/2/180)**2 + V*np.tan(np.pi*theta2[i]/2/180) + W)
                except:
                    fwhm = kwargs['FWHM']

                x = g2thetas - theta2[i]    
                tmp = self.split_type(x, xrd_intensity[i], N, fwhm, **kwargs) 

            else:
                msg = profile + ' is not supported'
                raise NotImplementedError(msg)
            gpeaks += tmp 

        gpeaks /= np.max(gpeaks)
        self.spectra = np.vstack((g2thetas, gpeaks))

    def gaussian_profile(self, I0, theta2, alpha, fwhm):
        tmp = ((alpha - theta2)/fwhm)**2
        return I0 * np.exp(-4*np.log(2)*tmp)
    
    def lorentzian_profile(self, I0, theta2, alpha,fwhm):
        tmp = 1 + 4*((alpha - theta2)/fwhm)**2
        return I0 * 1/tmp

    def split_type(self, x, I0, N, fwhm, **kwargs):
        tmp = np.zeros((N))
        for k, dx in enumerate(x):
            if dx < 0:
                A = kwargs['A']
                eta_l = kwargs['eta_l']
                eta_h = kwargs['eta_h']

            else:
                A = 1/kwargs['A'] 
                eta_l = kwargs['eta_h']
                eta_h = kwargs['eta_l']

            tmp[k] = ((1+A)*(eta_h + np.sqrt(np.pi*np.log(2))*(1-eta_h))) /\
                     (eta_l + np.sqrt(np.pi*np.log(2)) * (1-eta_l) + A*(eta_h + \
                     np.sqrt(np.pi*np.log(2))*(1-eta_h))) * (eta_l*2/(np.pi*fwhm) * \
                     (1+((1+A)/A)**2 * (dx/fwhm)**2)**(-1) + (1-eta_l)*np.sqrt(np.log(2)/np.pi)* \
                     2/fwhm *np.exp(-np.log(2) * ((1+A)/A)**2 * (dx/fwhm)**2))

        return I0 * tmp

    def pxrdf(self):
        """
        Group the equivalent hkl planes together by 2\theta angle
        N*6 arrays, Angle, d_hkl, h, k, l, intensity
        """
        
        rank = range(len(self.theta2)) #np.argsort(self.theta2)
        PL = []
        last = 0
        for i in rank:
            if self.xrd_intensity[i] > 0.01:
                angle = self.theta2[i]
                if abs(angle-last) < 1e-4:
                    PL[-1][-1] += self.xrd_intensity[i]
                else:
                    PL.append([angle, self.d_hkl[i], \
                             self.hkl_list[i][0]["hkl"][0], self.hkl_list[i][0]["hkl"][1], \
                             self.hkl_list[i][0]["hkl"][2], self.xrd_intensity[i]])
                last = angle

        PL = (np.array(PL))
        PL[:,-1] = PL[:,-1]/max(PL[:,-1])
        self.pxrd = PL
        # print(PL[0],PL[-1])
    
    def plot_pxrd(self, filename=None, minimum_I = 0.01, show_hkl=True):
        """ plot PXRD """

        plt.figure(figsize=(20,10))

        if self.profiling != None:
            plt.plot(self.gtwo_thetas,self.gpeaks,'g-',label = str(self.profiling) + ' profiling')

        dx = np.degrees(self.max2theta)
        for i in self.pxrd:
            plt.bar(i[0],i[-1], color='b', width=dx/180)
            if i[-1] > minimum_I:
               if show_hkl:
                  label = self.draw_hkl(i[2:5])
                  plt.text(i[0]-dx/40, i[-1], label[0]+label[1]+label[2])
        
        ax=plt.gca()
        plt.grid()
        plt.xlim(0,dx)
        plt.xlabel('2θ')
        plt.ylabel('Intensity')
        plt.title('PXRD of '+self.name+ ', $\lambda$='+str(self.wavelength)+'$\AA$')
        
        if filename is None:
           plt.show()
        else:
           plt.savefig(filename)
           plt.close()
        
    def get_unique_families(self,hkls):
        """
        Returns unique families of Miller indices. Families must be permutations
        of each other.
        Args:
            hkls ([h, k, l]): List of Miller indices.
        Returns:
            {hkl: multiplicity}: A dict with unique hkl and multiplicity.
        """

       # TODO: Definitely can be sped up.
        def is_perm(hkl1, hkl2):
            h1 = np.abs(hkl1)
            h2 = np.abs(hkl2)
            return all([i == j for i, j in zip(sorted(h1), sorted(h2))])

        unique = collections.defaultdict(list)
        for hkl1 in hkls:
            found = False
            for hkl2 in unique.keys():
                if is_perm(hkl1, hkl2):
                    found = True
                    unique[hkl2].append(hkl1)
                    break
            if not found:
                unique[hkl1].append(hkl1)

        pretty_unique = {}
        for k, v in unique.items():
            pretty_unique[sorted(v)[-1]] = len(v)

        return pretty_unique

    @staticmethod
    def draw_hkl(hkl):
        """turn negative numbers in hkl to overbar"""
        hkl_str= []
        for i in hkl:
            if i<0:
               label = str(int(-i))
               label = r"$\bar{" + label + '}$'
               hkl_str.append(str(label))
            else:
               hkl_str.append(str(int(i)))

        return hkl_str

from optparse import OptionParser
import pandas as pd
from tabulate import tabulate

if __name__ == "__main__":
    #-------------------------------- Options -------------------------
    parser = OptionParser()
    parser.add_option("-m", "--hkl", dest="hkl", metavar='hkl index',
                      help="show hkl_index info, e.g., [1,0,0]")
    parser.add_option("-a", "--angle", dest="max2theta", default=180, type='float',
                      help="2theta angle range, default=180", metavar="angle")
    parser.add_option("-t", "--transform", dest="trans", metavar="files",
                      help="export file in different format")
    parser.add_option("-p", "--plot", dest="plot", default='yes',
                      help="plot pxrd, default: yes", metavar="plot")
    parser.add_option("-w", "--wavelength", dest="wavelength", default=1.54184, type='float',
                      help="wavelength: 1.54184", metavar="wave")
    parser.add_option("-c", "--crystal", dest="structure",default='',
                      help="crystal from file, cif or poscar, REQUIRED", metavar="crystal")
    parser.add_option("-f", "--full", dest="full",default='no',
                      help="show full hkl reflections", metavar="full")
    parser.add_option("-i", "--intensity", dest="minimum_I",default=0.01, type='float',
                      help="the minimum intensity to show, default 0.01", metavar="intensity")


    (options, args) = parser.parse_args()    
    if options.structure.find('cif') > 0:
       fileformat = 'cif'
    else:
       fileformat = 'vasp'

    test = read(options.structure, format=fileformat)
    if options.plot == 'yes' or options.hkl is not None:
       xrd = XRD(test, wavelength=options.wavelength, \
                       max2theta=options.max2theta)   
       if options.full in  ['no', 'No', 'NO']:
          col_name = {'2theta': xrd.pxrd[:,0], \
                      'd_hkl':  xrd.pxrd[:,1], \
                      'h': xrd.pxrd[:,2], \
                      'k': xrd.pxrd[:,3], \
                      'l': xrd.pxrd[:,4], \
                      'Intensity':xrd.pxrd[:,5]}
       else:
          rank1 = xrd.xrd_intensity > options.minimum_I
          col_name = {'2theta':    np.degrees(xrd.theta2[rank1]), \
                      'd_hkl':     xrd.d_hkl[rank1],\
                      'h':         xrd.hkl_list[rank1,0], \
                      'k':         xrd.hkl_list[rank1,1], \
                      'l':         xrd.hkl_list[rank1,2], \
                      'Intensity': xrd.xrd_intensity[rank1] }

       df = pd.DataFrame(col_name)
       print(tabulate(df, headers='keys')) #, tablefmt='psql'))

       if options.plot == 'yes':
          xrd.plot_pxrd(filename=options.structure+'.png', minimum_I = options.minimum_I)
