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
        # Save uploaded CIF to disk
        f = form.cif.data
        save_path = os.path.join(app.instance_path, 'CIF', secure_filename(f.filename))
        f.save(save_path)
        
        # Retrieve form data
        wavelength = form.wavelength.data
        max2theta = form.theta.data
        N = 10000
        # Alert w/ submission data
        flash('XRD calculated for {}, λ={} Å, 2θ={}°'.format(f.filename, wavelength, max2theta))

        # Calculate and plot
        '''QZ: here is an example if you use split-type'''
        U = 5.776410E-03 # FWHM parameter, U
        V = -1.673830E-03 # FWHM parameter, V
        W = 5.668770E-03 # FWHM parameter, W
        A = 1.03944 # Asymmetry parameter, a1
        eta_h = 0.504656 # Mixing parameter, eta_H0
        eta_l = 0.611844  # Mixing parameter, eta_L0
        profile = {'function':'split-type', 'theta_dependence': True, 'U': U, 'V':V, 'W':W, 'A':A, 'eta_h':eta_h, 'eta_l':eta_l}

        struct = read(save_path, format='cif')
        xrd = XRD(struct, wavelength=wavelength, max2theta=max2theta) 
        xrd.get_profile(xrd.theta2, xrd.xrd_intensity, N, **profile)
        plot = xrd.plotly_pxrd()

        return render_template('index.html', title='Calculator', form=form, plot=plot)
    return render_template('index.html', title='Calculator', form=form)
