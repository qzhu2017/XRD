from flask_wtf import FlaskForm
from flask_wtf.file import FileField, FileRequired, FileAllowed
from wtforms import TextAreaField, StringField, BooleanField, SubmitField
from wtforms.validators import DataRequired

class CalcForm(FlaskForm):
    cif = FileField('CIF', validators=[FileRequired(), FileAllowed(['cif'], '.CIF (Crystallographic Information Files) only!')])
    wavelength = StringField('Wavelength (&#8491;)', validators=[DataRequired()])
    theta = StringField('Max. 2&theta; value', validators=[DataRequired()])
    submit = SubmitField('Calculate')
    
# Create separate classes for each profiling branch w/ member parameters
