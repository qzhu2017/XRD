import os
from flask import render_template, flash, session, Markup
from app import app
from app.forms import CalcForm
from werkzeug.utils import secure_filename
from pyxtal_xrd.XRD import XRD
from ase.io import read

@app.route('/', methods=['GET', 'POST'])
@app.route('/index', methods=['GET', 'POST'])
def index():
    form = CalcForm()
    if form.validate_on_submit(): # passed basic validation
        if form.upload.data:
            process_upload(form)
        if not session.get("SAVEPATH"): # new session
            flash('ERROR: Please upload a .CIF\
                    or -POSCAR file!', 'danger')
            return render_template('index.html',
                title='Calculator',
                form=form)
        else:
            process_form(form)
            return render_template('index.html', 
                title='Calculator',
                form=form,
                plot=plot(form))
    # initial page visit
    return render_template('index.html',
        title='Calculator',
        form=form)

def process_upload(form):
    """
    Save upload, check validity, and update session.
    """
    # Save uploaded file to disk
    f = form.upload.data
    savepath = os.path.join(app.instance_path, 
        'uploads', secure_filename(f.filename))
    f.save(savepath)

    # Check readability
    try:
        read(savepath) # attempt ase.io.read

        # Update session keys
        session["FILENAME"] = f.filename
        session["SAVEPATH"] = savepath
        flash('Successfully processed XRD for {}'.format(
            session.get("FILENAME")), 'success')
    except:
        flash('ERROR: Unable to read file.\
            Please try again or another file.',
            'danger')

def process_form(form):
    """
    Advanced form validation and session update.
    """
    # Retrieve form data
    wavelength = form.wavelength.data
    max2theta = form.max2theta.data
    min2theta = form.min2theta.data
    res = form.res.data
    profile = form.profile.data

    # Advanced form-level validation
    if min2theta > max2theta:
        min2theta = 0 # use default
        flash(Markup('WARNING: 2&theta;<sub>min</sub>\
            greater than\
            2&theta;<sub>max</sub>&mdash;defaulting\
            2&theta;<sub>min</sub> to 0&deg;.'), 'warning')

    # Update session keys
    session["WAVELENGTH"] = wavelength
    session["MIN2THETA"] = min2theta
    session["MAX2THETA"] = max2theta
    session["RES"] = res
        
def plot(form):
    """
    Process and return PXRD plotly.
    """
    struct = read(session.get("SAVEPATH"))
    xrd = XRD(struct,
        wavelength=session.get("WAVELENGTH"),
        thetas=[session.get("MIN2THETA"),
            session.get("MAX2THETA")]) 

    xrd.get_profile(res=session.get("RES"))
    return xrd.plotly_pxrd()
