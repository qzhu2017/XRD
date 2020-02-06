import os

class Config(object):
    SECRET_KEY = os.environ.get('SECRET_KEY') or '<X}9~B%~c*Nm77WX'
    