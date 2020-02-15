import os
from flask import render_template, flash, session
from app import app
from app.forms import CalcForm
from werkzeug.utils import secure_filename
from pyxtal_xrd.XRD import XRD
from ase.io import read

@app.route('/', methods=['GET', 'POST'])
@app.route('/index', methods=['GET', 'POST'])
def index():
    form = CalcForm()
    if form.validate_on_submit():
        # Brand new session
        if not session.get("SAVEPATH"):
            if form.upload.data:
                save_upload(form)
                return render_template('index.html', 
                    title='Calculator',
                    form=form,
                    plot=plot(form))
            else:
                flash('Please upload a .CIF or -POSCAR file!',
                    'danger')
                return render_template('index.html',
                    title='Calculator',
                    form=form)
        # Existing upload
        else:
            if form.upload.data: # newer upload
                save_upload(form) 
            return render_template('index.html', 
                title='Calculator',
                form=form,
                plot=plot(form))
    return render_template('index.html',
        title='Calculator',
        form=form)

def save_upload(form):
    """
    Convenience function to save uploads and update session
    """
    # Save uploaded file to disk
    f = form.upload.data
    savepath = os.path.join(app.instance_path, 
        'uploads', secure_filename(f.filename))
    f.save(savepath)
    # Update session data
    session["FILENAME"] = f.filename
    session["SAVEPATH"] = savepath
        
def plot(form):
    """
    Convenience function to process and return PXRD plotly
    """
    # Retrieve form data
    wavelength = form.wavelength.data
    max2theta = form.theta.data
    min2theta = 0 # form.theta.data1
    N = 10000

    struct = read(session.get("SAVEPATH"))
    xrd = XRD(struct,
        wavelength=wavelength,
        thetas=[min2theta, max2theta]) 
    xrd.get_profile(res=0.01)

    flash('XRD for {} plotted below.'.format(
        session.get("FILENAME")), 'success')

    return xrd.plotly_pxrd()

# retrieve file extension
# ext = os.path.splitext(savepath)[1]
