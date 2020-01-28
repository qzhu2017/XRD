from flask_wtf import FlaskForm
from wtforms import TextAreaField, StringField, BooleanField, SubmitField
from wtforms.validators import DataRequired

class CalcForm(FlaskForm):
    cif = TextAreaField('CIF', validators=[DataRequired()])
    wavelength = StringField('Wavelength (&#8491;)', validators=[DataRequired()])
    theta = StringField('Max. 2&theta; value', validators=[DataRequired()])
    submit = SubmitField('Calculate')
    