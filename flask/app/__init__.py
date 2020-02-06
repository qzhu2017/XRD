from flask import Flask
from config import Config
from flask_bootstrap import Bootstrap
from flask_uploads import patch_request_class

app = Flask(__name__)
app.config.from_object(Config)
bootstrap = Bootstrap(app)
patch_request_class(app, size=16777216)

from app import routes, errors
