from flask import render_template, flash
from app import app
from app.forms import CalcForm

@app.route('/', methods=['GET', 'POST'])
@app.route('/index', methods=['GET', 'POST'])
def index():
    form = CalcForm()
    if form.validate_on_submit():
        flash('XRD Calculation for CIF, Wavelength={}Å, 2θ={}°'.format(form.wavelength.data, form.theta.data))
        return render_template('index.html', title='Calculator', form=form)
    return render_template('index.html', title='Calculator', form=form)
