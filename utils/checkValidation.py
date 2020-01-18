import numpy as np 
import matplotlib.pyplot as plt
from XRD import crystal, Element, XRD
import sys
import os
from similarity import Similarity

def getListOfFiles(dirName):
    # create a list of file and sub directories 
    # names in the given directory 
    listOfFile = os.listdir(dirName)
    allFiles = []
    # Iterate over all the entries
    for entry in listOfFile:
        # Create full path
        fullPath = os.path.join(dirName, entry)
        # If entry is a directory then get the list of files in this directory 
        if os.path.isdir(fullPath):
            allFiles = allFiles + getListOfFiles(fullPath)
        else:
            allFiles.append(fullPath)
                
    return allFiles  

path = '/Users/sayred1/Research/dXRD/XRD/datatest'
allFiles = getListOfFiles(path)
allFiles.sort()
poscarFiles = []
cifFiles = []
diffFiles = []

for file in allFiles:
    if 'POSCAR' in file:
        poscarFiles.append(file)
    elif 'cif' in file:
        cifFiles.append(file)
    elif 'diff' in file:
        diffFiles.append(file)

"""
- Take the POSCAR file, load it into pxrd
- plot pxrd diffraction dataset against DIFF file
"""

wavelength = 1.54056
max2theta = 90
fwhm = 0.9
N = 10000
profile = 'gaussian'
sim = []
count = 1
for poscardata, diffdata in zip(poscarFiles,diffFiles):
    
    """
    Run .cif files through XRD, get profile
    """
    struct = crystal('POSCAR',filename=poscardata)
    xrd1 = XRD(struct, wavelength, max2theta)   
    xrd1.get_profile(xrd1.theta2, xrd1.xrd_intensity,N, profile, fwhm)
    fpeaks = xrd1.gpeaks
    
    """
    Load the diffraction data 
    """
    diff = np.loadtxt(diffdata,str,delimiter='\n')
    size = diff.shape[0]
    xval = []
    yval = []
    i = 0
    while i < size:
        if '2-THETA' in diff[i] and 'INTENSITY' in diff[i]:
            for j in range(i+1, size):
                try:
                    xval.append(float(diff[j].split()[0]))
                    yval.append(float(diff[j].split()[1]))
                except:
                    break
        i+=1
    
    """
    Get profile for diffraction data
    """
    xval = np.array(xval)
    yval = np.array(yval)
    yval/= np.max(yval)
    xrd2 = XRD(struct, wavelength, max2theta)
    xrd2.get_profile(xval, yval,N,profile,fwhm)
    gpeaks = xrd2.gpeaks
    
    S = Similarity(fpeaks, gpeaks, 0.6).calculate()
    sim.append(S)
    print(count,S)
    count +=1
with open('sim-validation.txt', 'w') as f:
    f.write(sim)


