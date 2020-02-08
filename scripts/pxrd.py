from pyxtal_xrd.XRD import XRD
from ase.io import read
from optparse import OptionParser

if __name__ == "__main__":
    #-------------------------------- Options -------------------------
    parser = OptionParser()
    parser.add_option("-a", "--angle", dest="max2theta", default=180, type='float',
                      help="2theta angle range, default=180", metavar="angle")
    parser.add_option("-w", "--wavelength", dest="wavelength", default=1.54184, type='float',
                      help="wavelength: 1.54184", metavar="wave")
    parser.add_option("-c", "--crystal", dest="structure",default='',
                      help="crystal from file, cif or poscar, REQUIRED", metavar="crystal")

    (options, args) = parser.parse_args()    
    if options.structure.find('cif') > 0:
       fileformat = 'cif'
    else:
       fileformat = 'vasp'


    U = 5.776410E-03 # FWHM parameter, U
    V = -1.673830E-03 # FWHM parameter, V
    W = 5.668770E-03 # FWHM parameter, W
    A = 1.03944 # Asymmetry parameter, a1
    eta_h = 0.504656 # Mixing parameter, eta_H0
    eta_l = 0.611844  # Mixing parameter, eta_L0
    profile = {'function':'split-type', 'theta_dependence': True, 'U': U, 'V':V, 'W':W, 'A':A, 'eta_h':eta_h, 'eta_l':eta_l}

    test = read(options.structure, format=fileformat)
    xrd = XRD(test, wavelength=options.wavelength, max2theta=options.max2theta)   
    xrd.get_profile(xrd.theta2, xrd.xrd_intensity, 1000, **profile)
    xrd.plotly_pxrd(html='1.html')


