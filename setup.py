from distutils.core import setup
import setuptools  # noqa
from os import path

this_directory = path.abspath(path.dirname(__file__))
with open(path.join(this_directory, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

exec(open('pyxtal_xrd/version.py').read())

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name="pyxtal_xrd",
    version=__version__,
    author="Dean Sayred and Qiang Zhu",
    author_email="qiang.zhu@unlv.edu",
    description="Python code for crystal XRD simulation",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/qzhu2017/XRD",
    packages=['pyxtal_xrd', 
              'pyxtal_xrd.database', 
              ],
    package_data={'pyxtal.database': ['*.csv', '*.json'],
                 },

    #scripts=['scripts/pyxtal_atom', 
    #         'scripts/pyxtal_test', 
    #         'scripts/pyxtal_symmetry',
    #         'scripts/pyxtal_molecule',
    #         ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    install_requires=[
        'numpy>=1.13.3', 
        'scipy>=1.1.0', 
        'matplotlib>=2.0.0',
        'ase>=3.18.0',
        ],
    python_requires='>=3.6.1',
    license='MIT',
)
