
import numpy as np
import matplotlib.pyplot as plt
from XRD import crystal, Element, XRD
from pyquery import PyQuery as pq
import mechanize
import os
import string
import sys


spaceGroup = np.loadtxt("sgs.csv",str) # load list of all possible space groups
unaccSpaceGroup = [] # list for space groups not accounted for 
seeds = []
counts = []

os.chdir('data') # dump files here
count = 0
counts.append(count)
seed = 0


# while count < np.shape(spaceGroup)[0]:
for group in spaceGroup:  
    count += 1
    
    # if count != counts[-1]:
    #     seed = 0

    # initialize/open main page for search
    br = mechanize.Browser() 
    br.open("http://rruff.geo.arizona.edu/AMS/")
    br.set_handle_robots(False)
    
    # select and write space group to cellparam control in br form
    br.select_form(nr=0)
    br.set_all_readonly(False)
    br["CellParam"] = "sg="+group

    # submit or search, brings to next page
    response1 = br.submit()

    # select form
    br.select_form(nr=0)
    br.set_all_readonly(False)
    
    # if result only returns a single file, skip
    if int(br.get_value("hid1")) == 1:
        unaccSpaceGroup.append(group)
        print(group,count)
        continue
    
    # random seed to choose file from list of structures
    seed = np.random.randint(0,int(br.get_value("hid1")))
    """    
    while seed in seeds:
        seed = np.random.randint(0,int(br.get_value("hid1")))
    
    seeds.append(seed)
    
    """
    
    file = 1 # to load cif data
    br.find_control("check[]").items[seed].selected=True # this form chooses the crystal
    br.find_control("down").items[file].selected=True # this form chooses the file
    name = "sg"+str(count)
    # submit form and write data to file
    response2 = br.submit(name = "downloadSelected")
    data1 = response2.read()
    cif = open(name+".cif","wb")
    cif.write(data1)
    cif.close()
    br.back()
    
    file = 2 # to load diff data
    br.select_form(nr=0)
    br.set_all_readonly(False)
    br.find_control("check[]").items[seed].selected=True # this form chooses the crystal
    br.find_control("down").items[file].selected=True # this form chooses the file
    response3 = br.submit()
    data2 = response3.read()
    diff = open(name+".diff","wb")
    diff.write(data2)
    diff.close()
    
    """
    Here we load the files to check the composition of the crystal
    """
    
    # structure = crystal('cif', filename = name+".cif")
