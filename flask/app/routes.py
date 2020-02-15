import os
from flask import render_template, flash
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
        # Save uploaded files to disk
        f = form.upload.data
        save_path = os.path.join(app.instance_path, 'uploads', secure_filename(f.filename))
        f.save(save_path)
        # ext = os.path.splitext(save_path)[1] # retrieve file extension

        # Retrieve form data
        wavelength = form.wavelength.data
        max2theta = form.theta.data
        min2theta = 0 # form.theta.data1
        N = 10000
        # Alert w/ submission data
        flash('XRD calculated for {}, λ={} Å, 2θ={}°'.format(f.filename, wavelength, max2theta))

        struct = read(save_path)
        xrd = XRD(struct, wavelength=wavelength, thetas=[0, max2theta]) 
        xrd.get_profile(res=0.01)
        plot = xrd.plotly_pxrd()

        return render_template('index.html', title='Calculator', form=form, plot=plot)
    return render_template('index.html', title='Calculator', form=form)
