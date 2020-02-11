import numpy as np
import matplotlib.pyplot as plt
import json
import os
import collections
from ase.io import read
from pyxtal_xrd.database.element import Element
from pyxtal_xrd.profile import Profile

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

    def __init__(self, crystal, wavelength=1.54184, 
                 thetas = [0, 180], 
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
        self.min2theta = np.radians(thetas[0])
        self.max2theta = np.radians(thetas[1])
        self.name = crystal.get_chemical_formula()
        self.preferred_orientation = preferred_orientation
        self.march_parameter = march_parameter
        self.all_dhkl(crystal)
        self.intensity(crystal)
        self.pxrdf()     
        

    def get_profile(self, method='pseudo_voigt', res=0.01):
        self.spectra = Profile(method, res).get_profile(self.theta2, self.xrd_intensity, 
                               np.degrees(self.min2theta), np.degrees(self.max2theta))

    def by_hkl(self, hkl):
        
        # this is a simple print statement, does not need to be optimized

        """ d for any give abitray [h,k,l] index """
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
        fp = os.path.join(os.path.dirname(__file__), "database/atomic_scattering_params.json")
        with open(fp, 'r') as f:
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

    def plotly_pxrd(self, minimum_I = 0.01, html=None):
        import plotly.graph_objects as go
        from plotly.subplots import make_subplots

        """
        interactive plot for pxrd powered by plotly
        Args:
        xrd: xrd object
        html: html filename (str)
        """

        
        x, y, labels = [], [], []
        for i in range(len(self.pxrd)):
            theta2, d, h, k, l, intensity = self.pxrd[i]
            h, k, l = int(h), int(k), int(l)
            if intensity > minimum_I:
                label = '<br>2&#952;: {:6.2f}<br>d: {:6.4f}</br>hkl: ({:d}{:d}{:d})'.format(theta2, d, h, k, l)
                x.append(theta2)
                y.append(-0.1)
                labels.append(label)

        trace1 = go.Bar(x=x, y=y, text=labels,
                        hovertemplate = "%{text}",
                        width=0.2, name='Index')
        trace2 = go.Scatter(x=self.spectra[0], y=self.spectra[1], name='Profile')

        fig = go.Figure(data=[trace2, trace1])
        fig.update_layout(xaxis_title = '2&#952; ({:.4f} &#8491;)'.format(self.wavelength),
                          yaxis_title = 'Intensity',
                          title = 'PXRD of '+self.name)

        if html is None:
            return fig.to_html()
        else:
            fig.write_html(html)
