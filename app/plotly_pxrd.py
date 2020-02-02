import plotly.graph_objects as go
from plotly.subplots import make_subplots

def plot_pxrd(xrd, minimum_I = 0.01, html=None):
    """
    interactive plot for pxrd powered by plotly
    Args:
    xrd: xrd object
    html: html filename (str)
    """

    x, y, labels = [], [], []
    for i in range(len(xrd.pxrd)):
        theta2, d, h, k, l, intensity = xrd.pxrd[i]
        h, k, l = int(h), int(k), int(l)
        if intensity > minimum_I:
            label = '<br>2&#952;: {:6.2f}<br>d: {:6.4f}</br>hkl: ({:d}{:d}{:d})'.format(theta2, d, h, k, l)
            x.append(theta2)
            y.append(-0.1)
            labels.append(label)
    trace1 = go.Bar(x=x, y=y, text=labels, 
                    hovertemplate = "%{text}",
                    width=0.2, name='Index')
    trace2 = go.Scatter(x=xrd.spectra[0], y=xrd.spectra[1], name='Profile')
    fig = go.Figure(data=[trace2, trace1])
    fig.update_layout(xaxis_title = '2&#952; ({:.4f} &#8491;)'.format(xrd.wavelength),
                      yaxis_title = 'Intensity')

    if html is None:
        return fig.to_html()
    else:
        fig.write_html(html)

if __name__ == "__main__":

    from XRD import crystal, XRD

    wavelength = 1.54056
    max2theta = 90
    N = 10000
    file1 = 'POSCAR-NaCl'
    struct1 = crystal('POSCAR', filename=file1)
    xrd1 = XRD(struct1, wavelength, max2theta) 
    xrd1.get_profile(xrd1.theta2, xrd1.xrd_intensity,N)
    plot_pxrd(xrd1, html='1.html')


