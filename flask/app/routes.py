import os
import plotly.graph_objects as go
from flask import render_template, flash, session, Markup
from app import app
from app.forms import MainForm
from werkzeug.utils import secure_filename
from pyxtal_xrd.XRD import XRD
from pyxtal_xrd.similarity import Similarity
from ase.io import read

@app.route('/', methods=['GET', 'POST'])
@app.route('/index', methods=['GET', 'POST'])
def index():
    form = MainForm()
    if form.validate_on_submit():
        if form.upload.data:
            process_upload(form)
        if not session.get("SAVEPATH"): # new session
            flash(Markup('<b>ERROR</b>: Please upload an\
                input file.'), 'danger')
            return render_template('index.html',
                title='Single',
                form=form)
        else:
            process_form(form)
            return render_template('index.html', 
                title='Single',
                form=form,
                plot=plot())
    # initial page visit
    return render_template('index.html',
        title='Single',
        form=form)

@app.route('/comparison', methods=['GET', 'POST'])
def comparison():
    form = MainForm()
    if form.validate_on_submit():
        if form.upload.data:
            process_upload(form)
        if form.upload2.data:
            process_upload(form, True)
        if not session.get("SAVEPATH")\
            or not session.get("SAVEPATH2"): # new session
            flash(Markup('<b>ERROR</b>: Please upload\
                <b>two</b> input files.'), 'danger')
            return render_template('comparison.html',
                title='Comparison',
                form=form)
        else:
            process_form(form, True)
            return render_template('comparison.html', 
                title='Comparison',
                form=form,
                plot=compare())
    # initial page visit
    return render_template('comparison.html',
        title='Comparison',
        form=form)

def process_upload(form, comp=False):
    """
    Save upload, check validity, and update session.
    """
    # Save uploaded file to disk
    if comp:
        f = form.upload2.data
    else:
        f = form.upload.data
    savepath = os.path.join(app.instance_path, 
        'uploads', secure_filename(f.filename))
    f.save(savepath)

    # Check readability
    try:
        read(savepath) # attempt ase.io.read

        # Update session keys
        if comp:
            session["FILENAME2"] = f.filename
            session["SAVEPATH2"] = savepath
        else:
            session["FILENAME"] = f.filename
            session["SAVEPATH"] = savepath
            flash(Markup('<b>{}</b> successfully\
                processed.').format(session.get("FILENAME")), 
                'success')
    except:
        flash(Markup('<b>ERROR</b>: Unable to read\
            <b>{}</b>. Please try again or a different\
            file.').format(f.filename), 'danger')

def process_form(form, comp=True):
    """
    Advanced form validation and session update.
    """
    # Retrieve form data
    max2theta = form.max2theta.data
    min2theta = form.min2theta.data

    # Advanced form-level validation
    if min2theta > max2theta:
        min2theta = 0 # use default
        flash(Markup('WARNING: 2&theta;<sub>min</sub>\
            <i>greater</i> than\
            2&theta;<sub>max</sub>&mdash;defaulting\
            2&theta;<sub>min</sub> to 0&deg;.'), 'warning')

    # Update session keys
    session["WAVELENGTH"] = form.wavelength.data
    session["MIN2THETA"] = min2theta
    session["MAX2THETA"] = max2theta
    session["RES"] = form.res.data
    session["METHOD"] = form.method.data
    session["FWHM"] = form.fwhm.data
    session["U"] = form.u.data
    session["V"] = form.v.data
    session["W"] = form.w.data
    session["A"] = form.a.data
    session["ETA_H"] = form.eta_h.data
    session["ETA_L"] = form.eta_l.data
    if comp:
        session["SHIFT"] = form.shift.data
        
def plot():
    """
    Process and return PXRD plotly.
    """
    method = session.get("METHOD")
    if method == 'gaussian' or method == 'lorentzian':
        kwargs = {
                    'FWHM': session.get("FWHM")
                }
    elif method == 'pseudo_voigt':
        kwargs = {
                    'U': session.get("U"), 
                    'V': session.get("V"),
                    'W': session.get("W"),
                    'A': session.get("A"),
                    'eta_h': session.get("ETA_H"),
                    'eta_l': session.get("ETA_L"),
                }

    struct = read(session.get("SAVEPATH"))
    xrd = XRD(struct,
        wavelength=session.get("WAVELENGTH"),
        thetas=[session.get("MIN2THETA"),
            session.get("MAX2THETA")]) 
    xrd.get_profile(method=method,
        res=session.get("RES"),
        user_kwargs=kwargs)
    flash(Markup('Showing <b>{}</b> with <i>{}</i>\
        profiling.').format(
            session.get("FILENAME"),
            method), 'info')
    return xrd.plotly_pxrd()

def compare():
    """
    Process and return comparison PXRD plotly.
    """
    method = session.get("METHOD")
    if method == 'gaussian' or method == 'lorentzian':
        kwargs = {
                    'FWHM': session.get("FWHM")
                }
    elif method == 'pseudo_voigt':
        kwargs = {
                    'U': session.get("U"), 
                    'V': session.get("V"),
                    'W': session.get("W"),
                    'A': session.get("A"),
                    'eta_h': session.get("ETA_H"),
                    'eta_l': session.get("ETA_L"),
                }

    files = [session.get("FILENAME"),
            session.get("FILENAME2")]
    structs = [read(session.get("SAVEPATH")),
                read(session.get("SAVEPATH2"))]
    xrds = []

    for struct in structs:
        xrd = XRD(struct,
                wavelength=session.get("WAVELENGTH"),
                thetas=[session.get("MIN2THETA"),
                session.get("MAX2THETA")])
        xrd.get_profile(method=method,
            res=session.get("RES"),
            user_kwargs=kwargs)
        xrds.append(xrd)

    S = Similarity(xrds[0].spectra,
        xrds[1].spectra,
        l=session.get("SHIFT"))

    S.calculate()
    title = 'PXRD Similarity {:6.3f} with shift\
        {:6.3f}'.format(S.S, S.l)
    traces = []

    for i, xrd in enumerate(xrds):
        traces.append(go.Scatter(x=xrd.spectra[0],
            y=xrd.spectra[1],
            name=str(files[i])))
    
    fig = go.Figure(data=traces)
    fig.update_layout(xaxis_title = '2&#952; ({:.4f}\
        &#8491;)'.format(session.get("WAVELENGTH")),
        yaxis_title = 'Intensity',
        title_text = title, 
        title_x=0.5)
    flash(Markup('Comparing <b>{}</b> and <b>{}</b> with\
        <i>{}</i> profiling.').format(
            session.get("FILENAME"),
            session.get("FILENAME2"),
            method), 'info')
    return fig.to_html()
