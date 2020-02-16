from flask_wtf import FlaskForm
from flask_wtf.file import FileField, FileRequired, FileAllowed
from wtforms import FloatField, SelectField, SubmitField
from wtforms.validators import DataRequired, NumberRange, ValidationError

# from ase.io import read
# import os
# from app import app
# from werkzeug.utils import secure_filename

class CalcForm(FlaskForm):
    upload = FileField(
        label='CIF/POSCAR',
        description='Upload a .CIF or -POSCAR file')
    wavelength = FloatField(
        label='&lambda; (&#8491;)',
        validators=[
            DataRequired(),
            NumberRange(
                min=0.1,
                max=5,
                message='Must be between %(min)s Å and %(max)s Å')],
        description='X-ray wavelength in angstroms',
        default=1.54056)
    min2theta = FloatField(
        label='2&theta;<sub>min</sub> (&deg;)',
        validators=[
            NumberRange(
                min=0,
                max=180,
                message='Must be between %(min)s° and %(max)s°')],
        description='Minimum diffraction angle in degrees',
        default=0)
    max2theta = FloatField(
        label='2&theta;<sub>max</sub> (&deg;)',
        validators=[
            DataRequired(),
            NumberRange(
                min=5,
                max=180,
                message='Must be between %(min)s° and %(max)s°')],
        description='Maximum diffraction angle in degrees',
        default=90)
    res = FloatField(
        label='Resolution',
        validators=[
            DataRequired(),
            NumberRange(
                min=1e-3,
                max=1,
                message='Must be between %(min)s and %(max)s')
            ],
        description='Profiling resolution',
        default=0.01)
    profile = SelectField(
        label='Profile',
        choices=[('gaussian', 'Gaussian'), ('lorentzian', 'Lorentzian'), ('split-type', 'Split-type')],
        description='Profiling function applied to simulated XRD pattern')
    submit = SubmitField('Visualize')
    
# Create separate classes for each profiling branch w/ member parameters (see "Field Enclosures")
