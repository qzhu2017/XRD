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


    test = read(options.structure, format=fileformat)
    xrd = XRD(test, wavelength=options.wavelength, max2theta=options.max2theta)   
    xrd.get_profile(res=0.01)
    xrd.plotly_pxrd(html='1.html')


