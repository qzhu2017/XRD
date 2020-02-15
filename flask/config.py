import os

class Config(object):
    SECRET_KEY = os.environ.get('SECRET_KEY') or 'tR4rcW_A2JrGmZWWj-UtLA'
    