# XRD
X-ray diffraction calculations.

## VXRD: Virtual X-Ray Diffraction
In order to run and view VXRD locally, run the following shell commands:
```bash
$ cd XRD
$ pip install -r requirements.txt
$ cd ./flask
$ flask run
```
It's **important** to rerun `pip install -r requirements.txt` if this repository has been updated in case there are any new dependencies.

If everything is setup correctly, you should see the following output:
```bash
 * Serving Flask app "vxrd.py"
 * ...
 * Running on http://127.0.0.1:5000/ (Press CTRL+C to quit)
```
Then, open your web browser and enter the following URL:
`http://localhost:5000/`

When finished, press `CTRL+C` in your terminal to shutdown the web-app.

### VXRD: JSmol
In order to see the 3D structure visualized with JSmol, you'll need to unzip `jsmol.zip` into the following directory:
```bash
$ unzip jsmol.zip ./flask/app/static/
```
