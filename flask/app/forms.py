from flask_wtf import FlaskForm
from flask_wtf.file import FileField, FileRequired, FileAllowed
from wtforms import FloatField, SelectField, SubmitField
from wtforms.validators import DataRequired, NumberRange, ValidationError
# from wtforms import StringField, FieldList, FormField

class CalcForm(FlaskForm):
    upload = FileField(
        label='Input File (CIF, POSCAR, etc.)',
        description='See <a href="https://wiki.fysik.dtu.dk/ase/ase/io/io.html?highlight=formats#file-input-and-output" target="_blank">table</a> for a list of readable formats')
    wavelength = FloatField(
        label='<i>&lambda;</i> (&#8491;)',
        validators=[
            DataRequired(),
            NumberRange(
                min=0.1,
                max=5,
                message='Must be between %(min)s Å and %(max)s Å')],
        description='X-ray wavelength in angstroms',
        default=1.54056)
    min2theta = FloatField(
        label='2<i>&theta;</i><sub>min</sub> (&deg;)',
        validators=[
            NumberRange(
                min=0,
                max=180,
                message='Must be between %(min)s° and %(max)s°')],
        description='Minimum diffraction angle in degrees',
        default=0)
    max2theta = FloatField(
        label='2<i>&theta;</i><sub>max</sub> (&deg;)',
        validators=[
            DataRequired(),
            NumberRange(
                min=5,
                max=180,
                message='Must be between %(min)s° and %(max)s°')],
        description='Maximum diffraction angle in degrees',
        default=90)
    res = FloatField(
        label='Resolution (&deg;)',
        validators=[
            DataRequired(),
            NumberRange(
                min=1e-3,
                max=1,
                message='Must be between %(min)s° and %(max)s°')
            ],
        description='Resolution in degrees',
        default=0.01)
    method = SelectField(
        label='Profiling Function',
        choices=[('gaussian', 'Gaussian'),
            ('lorentzian', 'Lorentzian'),
            ('pseudo_voigt', 'Pseudo-Voigt')],
        description='Applied to simulated XRD pattern')
    fwhm = FloatField(
        label='FWHM',
        validators=[
            DataRequired(),
            NumberRange(
                min=1e-3,
                max=1,
                message='Must be between %(min)s and %(max)s')
            ],
        description='Full width at half maximum',
        default=0.02)
    u = FloatField(
        label='<i>U</i>',
        validators=[
            DataRequired(),
            NumberRange(
                min=1e-3,
                max=1,
                message='Must be between %(min)s and %(max)s')
            ],
        description='Pseudo-Voigt parameter',
        default=5.776410E-03)
    v = FloatField(
        label='<i>V</i>',
        validators=[
            DataRequired(),
            # NumberRange(
            #     min=1e-3,
            #     max=1,
            #     message='Must be between %(min)s and %(max)s')
            ],
        description='Pseudo-Voigt parameter',
        default=-1.673830E-03)
    w = FloatField(
        label='<i>W</i>',
        validators=[
            DataRequired(),
            NumberRange(
                min=1e-3,
                max=1,
                message='Must be between %(min)s and %(max)s')
            ],
        description='Pseudo-Voigt parameter',
        default=5.668770E-03)
    a = FloatField(
        label='<i>A</i>',
        validators=[
            DataRequired(),
            NumberRange(
                min=0.1,
                max=5,
                message='Must be between %(min)s and %(max)s')
            ],
        description='Pseudo-Voigt parameter',
        default=1.03944)
    eta_h = FloatField(
        label='<i>&eta;</i><sub>h</sub>',
        validators=[
            DataRequired(),
            NumberRange(
                min=0.1,
                max=1,
                message='Must be between %(min)s and %(max)s')
            ],
        description='Pseudo-Voigt parameter',
        default=0.504656)
    eta_l = FloatField(
        label='<i>&eta;</i><sub>l</sub>',
        validators=[
            DataRequired(),
            NumberRange(
                min=0.1,
                max=1,
                message='Must be between %(min)s and %(max)s')
            ],
        description='Pseudo-Voigt parameter',
        default=0.611844)
    submit = SubmitField('Visualize')

class CompForm(FlaskForm):
    upload = FileField(
        label='2<sup>nd</sup> Input File',
        description='Second input to compare')
    shift = FloatField(
        label='Shift',
        validators=[NumberRange(
                min=0,
                max=1,
                message='Must be between %(min)s and %(max)s')
            ],
        description='Shift for similarity',
        default=2.0)

# Create separate classes for each profiling branch
# w/ member parameters (see "Field Enclosures")
# class IMForm(FlaskForm):
#     protocol = SelectField(choices=[('aim', 'AIM'), ('msn', 'MSN')])
#     username = StringField()

# class ContactForm(FlaskForm):
#     first_name  = StringField()
#     last_name   = StringField()
#     im_accounts = FieldList(FormField(IMForm))
    