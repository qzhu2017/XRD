from flask_wtf import FlaskForm
from flask_wtf.file import FileField, FileRequired, FileAllowed
from wtforms import FloatField, SubmitField
from wtforms.validators import DataRequired, NumberRange, ValidationError

# from ase.io import read
# import os
# from app import app
# from werkzeug.utils import secure_filename

class CalcForm(FlaskForm):
    # struct = [] # try new attribute

    upload = FileField(
        label='CIF/POSCAR',
        validators=[
            FileRequired(),
            # FileAllowed(
            #     ['cif', ''],
            #     message='.CIF (Crystallographic Information Files) and -POSCAR only!')
            ],
        description='Upload a Crystallographic Information File (.CIF) or -POSCAR')
    wavelength = FloatField(
        label='&lambda; (&#8491;)',
        validators=[
            DataRequired(),
            NumberRange(
                min=0.1,
                max=5,
                message='Must be between %(min)s Å and %(max)s Å')],
        description='Wavelength in angstroms',
        default=1.54056)
    theta = FloatField(
        label='2&theta;<sub>max</sub> (&deg;)',
        validators=[
            DataRequired(),
            NumberRange(
                min=5,
                max=180,
                message='Must be between %(min)s° and %(max)s°')],
        description='Maximum diffraction angle in degrees',
        default=90)
    submit = SubmitField('Calculate')

    # Tried introspective validator
    # def validate_upload(self, upload):
    #     try:
    #         f = upload.data
    #         savepath = os.path.join(app.instance_path, 'uploads', secure_filename(f.filename))
    #         f.save(savepath)
    #         struct.append(read(savepath))
    #     except:
    #         raise ValidationError('.CIF and -POSCAR files only!')
    
# Create separate classes for each profiling branch w/ member parameters (see "Field Enclosures")
