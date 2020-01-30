import os
from flask import render_template, flash
from app import app
from app.forms import CalcForm
from werkzeug.utils import secure_filename

@app.route('/', methods=['GET', 'POST'])
@app.route('/index', methods=['GET', 'POST'])
def index():
    form = CalcForm()
    if form.validate_on_submit():
        f = form.cif.data
        filename = secure_filename(f.filename)
        f.save(os.path.join(
            app.instance_path, 'CIF', filename
        ))
        flash('XRD Calculation for {}, Wavelength={}Å, 2θ={}°'.format(f.filename, form.wavelength.data, form.theta.data))
        return render_template('index.html', title='Calculator', form=form)
    return render_template('index.html', title='Calculator', form=form)
