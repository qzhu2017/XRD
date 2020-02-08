# a test script to calculate the similarity between two structures 

from scripts.compare import Compare

path = './dataset/'

Compare(path+'A-1-POSCAR', 'POSCAR', path+'A-1-POSCAR', 'POSCAR')
