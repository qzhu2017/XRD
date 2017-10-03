# XRD
X-ray diffraction calculations.

This branch is initally planned for my Computational Physics class in 2017 fall.
https://github.com/qzhu2017/2017-cmp

Currently, there are three classes,
- Element
- crystal
- XRD

One could load the crystal from 
- dictionary
- POSCAR
- CIF file (to add in near future)

To perform XRD calculation, one needs to provide the following info
- crystal structure
- wavelength (default is Cu-Ka: 1.54184 \AA)
- maximum 2\theta value (defult: 180 degree)

The atomic scattering factor is calculated from 9-parameter equation by Don Cromer and J. Mann.

More detailed usage could be found in the jupyter notebook.
