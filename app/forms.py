from flask_wtf import FlaskForm
from flask_wtf.file import FileField, FileRequired, FileAllowed
from wtforms import FloatField, SubmitField
from wtforms.validators import DataRequired, NumberRange

class CalcForm(FlaskForm):
    cif = FileField(
        label='CIF',
        validators=[FileRequired(), FileAllowed(
            ['cif'],
            message='.CIF (Crystallographic Information Files) only!')],
        description='Upload a Crystallographic Information File (.CIF)')
    wavelength = FloatField(
        label='&lambda; (&#8491;)',
        validators=[DataRequired(), NumberRange(
            min=0.1,
            max=50,
            message='Must be between %(min)s Å and %(max)s Å'
        )],
        description='Wavelength in angstroms',
        default=1.54056)
    theta = FloatField(
        label='2&theta;<sub>max</sub> (&deg;)',
        validators=[DataRequired(), NumberRange(
            min=28,
            max=180,
            message='Must be between %(min)s° and %(max)s°'
        )],
        description='Maximum diffraction angle in degrees',default=90)
    submit = SubmitField('Calculate')
    
# Create separate classes for each profiling branch w/ member parameters (see "Field Enclosures")
