# python module to check the similarity between two loaded diffraction profiles
# a test script to calculate the similarity between to structures

import numpy as np
from pyxtal_xrd.XRD import Element, crystal, XRD
from pyxtal_xrd.profile import Profile
from pyxtal_xrd.similarity import Similarity

class Compare:

    def __init__(self, path_to_file1, file_type1, path_to_file2, file_type2):
        
        self.struct1 = crystal(file_type1, filename = path_to_file1)
        self.struct2 = crystal(file_type2, filename = path_to_file2)
    
        self.getPeaks()
        self.Profile()
        self.calculateSimilarity()

    def getPeaks(self):

        self.xrd1 = XRD(self.struct1)
        self.xrd2 = XRD(self.struct2)

    def Profile(self):

        profile1 = Profile()
        profile1.get_profile(self.xrd1.theta2, self.xrd1.xrd_intensity)
        self.f = profile1.spectra

        profile2 = Profile()
        profile2.get_profile(self.xrd2.theta2, self.xrd2.xrd_intensity)
        self.g = profile2.spectra
        
    def calculateSimilarity(self):
        
        similarity = Similarity(self.f, self.g)
        S = similarity.calculate()
        print(S)
