from pyxtal_xrd.XRD import XRD
from ase.io import read
from optparse import OptionParser
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from pyxtal_xrd.similarity import Similarity

if __name__ == "__main__":
    #-------------------------------- Options -------------------------
    parser = OptionParser()
    parser.add_option("-a", "--angle", dest="thetas", default='[0, 120]', type='str',
                      help="2theta angle range, default=180", metavar="angle")
    parser.add_option("-w", "--wavelength", dest="wavelength", default=1.54184, type='float',
                      help="wavelength: 1.54184", metavar="wave")
    parser.add_option("-c", "--crystal", dest="structure",default='',
                      help="crystal from file, cif or poscar, REQUIRED", metavar="crystal")
    parser.add_option("-l", "--shift", dest="shift",default=None, type='float',
                      help="shift for similarity", metavar="shift")


    (options, args) = parser.parse_args()    

    files = options.structure.replace('[','').replace(']','').split(',')
    strucs = []
    xrds = []
    for file in files:
        if file.find('cif') > 0:
            fileformat = 'cif'
        else:
            fileformat = 'vasp'
        strucs.append(read(file, format=fileformat))
    thetas = options.thetas.replace('[','').replace(']','')
    t = [float(i) for i in thetas.split(',')]

    for struc in strucs:
        xrd = XRD(struc, wavelength=options.wavelength, thetas=t)   
        xrd.get_profile(res=0.05)
        xrds.append(xrd)

    if len(xrds) == 1:
        xrd.plotly_pxrd(html='1.html')
    elif len(xrds) == 2:
        S = Similarity(xrds[0].spectra, xrds[1].spectra, l=options.shift, weight='triangle')
        S.calculate()
        title = 'PXRD Similarity {:6.4f} with shift {:6.4f}'.format(S.S, S.l)
        traces = []
        for i, xrd in enumerate(xrds):
            traces.append(go.Scatter(x=xrd.spectra[0], y=xrd.spectra[1], name=str(files[i])))
        fig = go.Figure(data=traces)
        fig.update_layout(xaxis_title = '2&#952; ({:.4f} &#8491;)'.format(options.wavelength),
                          yaxis_title = 'Intensity',
                          title_text = title, 
                          title_x=0.5)
        print(title)
        fig.write_html('1.html')
