import os

class Config(object):
    SECRET_KEY = os.environ.get('SECRET_KEY') or 'tR4rcW_A2JrGmZWWj-UtLA'

    MAIL_SERVER = 'smtp.gmail.com'
    MAIL_PORT = 587
    MAIL_USE_TLS = 1
    MAIL_USERNAME = 'vxrd.info'
    MAIL_PASSWORD = os.environ.get('MAIL_PASSWORD')
    ADMINS = ['vxrd.info@gmail.com']
