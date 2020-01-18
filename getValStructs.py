i == 1:
        with open(path + 'single.txt', 'a') as f:
            f.write(name + '\n')
        continue 

        
    found = False
    tmpseed = 0
    file = 0
    numAtoms  = []
    
    while tmpseed < top:
        br.find_control("check[]").items[tmpseed].selected=True # this form chooses the crystal
        br.find_control("down").items[file].selected=True # this form chooses the file
        response0 = br.submit()
        data0 = response0.read()
            
        with open(path+name+".amg","wb") as amgfile:
            amgfile.write(data0)
        with open(path+name+".amg","r") as amgfile:
            try:
                content = amgfile.readlines()
            except:
                numAtoms.append(0)
                br.back()
                os.remove(path+name+".amg")
                break
        content = [x.strip() for x in content]
        
        count = 0
        for j in content: 
            tmp = j.split()
            try:
                o = tmp.index("occ")
                numAtoms.append(0)
                break
            except:
                count += 1
            try: 
                a = tmp.index("atom")
                numAtoms.append(len(content) - count)
                break
            except:
                pass
                    
        
        tmpseed+=1
        br.back()
        br.select_form(nr=0)
        br.set_all_readonly(False)
        os.remove(path+name+".amg")
    
    if sum(numAtoms) == 0:
        with open(path+'partialOcc.txt', 'a') as f:
            f.write(name+ '\n')
        continue
    else:
        minNum = numAtoms[0]
        i = 0
        while minNum == 0:
            i+=1
            minNum = numAtoms[i]
        for i in range(1,len(numAtoms)):
            if numAtoms[i] < minNum and numAtoms[i] != 0:
                minNum = numAtoms[i]
    
    seed  = numAtoms.index(minNum)
    file = 1 
    br.select_form(nr=0)
    br.set_all_readonly(False)
    br.find_control("check[]").items[seed].selected=True # this form chooses the crystal
    br.find_control("down").items[file].selected=True # this form chooses the file
    response1 = br.submit()
    data1 = response1.read()
    with open(path+name+".cif","wb") as ciffile:    
        ciffile.write(data1)    
    
    br.back()
    
    file = 2 
    br.select_form(nr=0)
    br.set_all_readonly(False)
    br.find_control("check[]").items[seed].selected=True # this form chooses the crystal
    br.find_control("down").items[file].selected=True # this form chooses the file
    response2 = br.submit()
    data2 = response2.read()
    with open(path+name+".diff","wb") as diffile:
        diffile.write(data2)
    
    structure = mg.Structure.from_file(path+name+".cif")
    structure.to(filename=path+name+"-POSCAR")
