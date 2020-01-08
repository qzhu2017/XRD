from math import acos, pi, ceil
import numpy as np
import matplotlib.pyplot as plt
import re
from math import sin,cos,sqrt, degrees
import time
import sys
import json
import os
import collections

def angle(a,b):
    """ calculate the angle between vector a and b """
    return acos(np.dot(a,b)/np.linalg.norm(a)/np.linalg.norm(b))


class Element:
    def __init__(self, input_value):
        self.input = input_value

        # list with atomic number z, short name, full name, valence, 
                # valence electrons, covalent radius, good bonds, Maximum CN:
        self.elements_list = [
            (1, 'H', 'Hydrogen',    1.0, 1, 0.31),
            (2, 'He', 'Helium',     0.5, 2, 0.28),
            (3, 'Li', 'Lithium',    1.0, 1, 1.28),
            (4, 'Be', 'Beryllium',  2.0, 2, 0.96),
            (5, 'B', 'Boron',       3.0, 3, 0.84),
            (6, 'C', 'Carbon',      4.0, 4, 0.70),
            (7, 'N', 'Nitrogen',    3.0, 5, 0.71),
            (8, 'O', 'Oxygen',      2.0, 6, 0.66),
            (9, 'F', 'Fluorine',    1.0, 7, 0.57),
            (10, 'Ne', 'Neon',      0.5, 8, 0.58),
            (11, 'Na', 'Sodium',    1.0, 1, 1.66),
            (12, 'Mg', 'Magnesium', 2.0, 2, 1.41),
            (13, 'Al', 'Aluminium', 3.0, 3, 1.21),
            (14, 'Si', 'Silicon',   4.0, 4, 1.11),
            (15, 'P', 'Phosphorus', 3.0, 5, 1.07),
            (16, 'S', 'Sulfur',     2.0, 6, 1.05),
            (17, 'Cl', 'Chlorine',  1.0, 7, 1.02),
            (18, 'Ar', 'Argon',     0.5, 8, 1.06),
            (19, 'K', 'Potassium',  1.0, 1, 2.03),
            (20, 'Ca', 'Calcium',   2.0, 2, 1.76),
            (21, 'Sc', 'Scandium',  3.0, 3, 1.70),
            (22, 'Ti', 'Titanium',  4.0, 4, 1.60),
            (23, 'V', 'Vanadium',   4.0, 5, 1.53),
            (24, 'Cr', 'Chromium',  3.0, 6, 1.39),
            (25, 'Mn', 'Manganese', 4.0, 5, 1.39),
            (26, 'Fe', 'Iron',      3.0, 3, 1.32),
            (27, 'Co', 'Cobalt',    3.0, 3, 1.26),
            (28, 'Ni', 'Nickel',    2.0, 3, 1.24),
            (29, 'Cu', 'Copper',    2.0, 2, 1.32),
            (30, 'Zn', 'Zinc',      2.0, 2, 1.22),
            (31, 'Ga', 'Gallium',   3.0, 3, 1.22),
            (32, 'Ge', 'Germanium', 4.0, 4, 1.20),
            (33, 'As', 'Arsenic',   3.0, 5, 1.19),
            (34, 'Se', 'Selenium',  2.0, 6, 1.20),
            (35, 'Br', 'Bromine',   1.0, 7, 1.20),
            (36, 'Kr', 'Krypton',   0.5, 8, 1.16),
            (37, 'Rb', 'Rubidium',  1.0, 1, 2.20),
            (38, 'Sr', 'Strontium', 2.0, 2, 1.95),
            (39, 'Y', 'Yttrium',    3.0, 3, 1.90),
            (40, 'Zr', 'Zirconium', 4.0, 4, 1.75),
            (41, 'Nb', 'Niobium',   5.0, 5, 1.64),
            (42, 'Mo', 'Molybdenum',4.0, 6, 1.54),
            (43, 'Tc', 'Technetium',4.0, 5, 1.47),
            (44, 'Ru', 'Ruthenium', 4.0, 3, 1.46),
            (45, 'Rh', 'Rhodium',   4.0, 3, 1.42),
            (46, 'Pd', 'Palladium', 4.0, 3, 1.39),
            (47, 'Ag', 'Silver',    1.0, 2, 1.45),
            (48, 'Cd', 'Cadmium',   2.0, 2, 1.44),
            (49, 'In', 'Indium',    3.0, 3, 1.42),
            (50, 'Sn', 'Tin',       4.0, 4, 1.39),
            (51, 'Sb', 'Antimony',  3.0, 5, 1.39),
            (52, 'Te', 'Tellurium', 2.0, 6, 1.38),
            (53, 'I', 'Iodine',     1.0, 7, 1.39),
            (54, 'Xe', 'Xenon',     0.5, 8, 1.40),
            (55, 'Cs', 'Caesium',   1.0, 1, 2.44),
            (56, 'Ba', 'Barium',    2.0, 2, 2.15),
            (57, 'La', 'Lanthanum', 3.0, 3, 2.07),
            (58, 'Ce', 'Cerium',    4.0, 3, 2.04),
            (59,'Pr','Praseodymium',3.0, 3, 2.03),
            (60, 'Nd', 'Neodymium', 3.0, 3, 2.01),
            (61, 'Pm', 'Promethium',3.0, 3, 1.99),
            (62, 'Sm', 'Samarium',  3.0, 3, 1.98),
            (63, 'Eu', 'Europium',  3.0, 3, 1.98),
            (64, 'Gd', 'Gadolinium',3.0, 3, 1.96),
            (65, 'Tb', 'Terbium',   3.0, 3, 1.94),
            (66, 'Dy', 'Dysprosium',3.0, 3, 1.92),
            (67, 'Ho', 'Holmium',   3.0, 3, 1.92),
            (68, 'Er', 'Erbium',    3.0, 3, 1.89),
            (69, 'Tm', 'Thulium',   3.0, 3, 1.90),
            (70, 'Yb', 'Ytterbium', 3.0, 3, 1.87),
            (71, 'Lu', 'Lutetium',  3.0, 3, 1.87),
            (72, 'Hf', 'Hafnium',   4.0, 3, 1.75),
            (73, 'Ta', 'Tantalum',  5.0, 3, 1.70),
            (74, 'W', 'Tungsten',   4.0, 3, 1.62),
            (75, 'Re', 'Rhenium',   4.0, 3, 1.51),
            (76, 'Os', 'Osmium',    4.0, 3, 1.44),
            (77, 'Ir', 'Iridium',   4.0, 3, 1.41),
            (78, 'Pt', 'Platinum',  4.0, 3, 1.36),
            (79, 'Au', 'Gold',      1.0, 3, 1.36),
            (80, 'Hg', 'Mercury',   2.0, 3, 1.32),
            (81, 'Tl', 'Thallium',  3.0, 3, 1.45),
            (82, 'Pb', 'Lead',      4.0, 4, 1.46),
            (83, 'Bi', 'Bismuth',   3.0, 5, 1.48),
            (84, 'Po', 'Polonium',  2.0, 6, 1.40),
            (85, 'At', 'Astatine',  1.0, 7, 1.50),
            (86, 'Rn', 'Radon',     0.5, 8, 1.50),
            (87, 'Fr', 'Francium',  1.0, 1, 2.60),
            (88, 'Ra', 'Radium',    2.0, 2, 2.21),
            (89, 'Ac', 'Actinium',  3.0, 3, 2.15),
            (90, 'Th', 'Thorium',   4.0, 3, 2.06),
            (91,'Pa','Protactinium',4.0, 3, 2.00),
            (92, 'U', 'Uranium',    4.0, 3, 1.96),
            (93, 'Np', 'Neptunium', 4.0, 3, 1.90),
            (94, 'Pu', 'Plutonium', 4.0, 3, 1.87),
            (95, 'Am', 'Americium', 4.0, 3, 1.80),
            (96, 'Cm', 'Curium',    4.0, 3, 1.69),
            (97, 'Bk', 'Berkelium', 4.0, 3, None),
            (98,'Cf','Californium', 4.0, 3, None),
            (99,'Es','Einsteinium', 4.0, 3, None),
            (100, 'Fm', 'Fermium',  4.0, 3, None),
            (101,'Md','Mendelevium',4.0, 3, None),
            (102, 'No', 'Nobelium', 4.0, 3, None),
            (103, 'Lr','Lawrencium',4.0, 3, None),
            (104,'Rf','Rutherfordium',4.0,3,None),
            (105, 'Db', 'Dubnium',  2.0, 3, None),
        ]
       

        self.sf = None
        self.z = None
        self.short_name = None
        self.long_name = None
        self.valence = None
        self.valence_electrons = None
        self.covalent_radius = None

        pos = None

        try:
            int(self.input)
            self.z = self.input

            for i, el in enumerate(self.elements_list):
                if el[0] == self.z:
                    pos = i
                    self.short_name = el[1]
                    self.long_name = el[2]
                    break
        except ValueError:
            self.short_name = self.input
            for i, el in enumerate(self.elements_list):
                if el[1] == self.short_name:
                    pos = i
                    self.z = el[0]
                    self.long_name = el[2]
                    break

            if not self.z:
                self.short_name = None
                self.long_name = self.input
                for i, el in enumerate(self.elements_list):
                    if el[2] == self.long_name:
                        pos = i
                        self.z = el[0]
                        self.short_name = el[1]
                        break
                if not self.z:
                    self.long_name = None

        if pos is not None:
            self.valence = self.elements_list[pos][3]
            self.valence_electrons = self.elements_list[pos][4]
            self.covalent_radius = self.elements_list[pos][5]

    def get_all(self, pos):
        els = []
        for el in self.elements_list:
            els.append(el[pos])
        return els

    def get_sf(self, pos):
        with open(os.path.join(os.path.dirname(__file__),
               "atomic_scattering_params.json")) as f:
                ATOMIC_SCATTERING_PARAMS = json.load(f)
                els = ATOMIC_SCATTERING_PARAMS[pos]
        return els

    def all_z(self):
        return self.get_all(0)
    def all_short_names(self):
        return self.get_all(1)
    def all_long_names(self):
        return self.get_all(2)
    def all_valences(self):
        return self.get_all(3)
    def all_valence_electrons(self):
        return self.get_all(4)
    def all_covalent_radii(self):
        return self.get_all(5)
    def get_sf(self):
        return self.get_sf()

class crystal(object):
    """a class of crystal structure. 
    Attributes:
        cell_para: a,b,c, alpha, beta, gamma
        cell_matrix: 3*3 matrix
        rec_matrix: reciprocal of cell matrix
        atom_type:  elemental type (e.g. Na Cl)
        composition: chemical composition (e.g., [1,1])
        coordinate: atomic positions (e.g., [[0,0,0],[0.5,0.5,0.5]])
    """

    def __init__(self, fileformat='POSCAR', filename=None, \
                 lattice=None, atom_type=None, composition=None, coordinate=None):
        """Return a structure object with the proper structures info"""
        if fileformat == 'POSCAR':
           self.from_POSCAR(filename)
        elif fileformat == 'cif':
           self.from_cif(filename)
        else:
           self.from_dict(lattice, atom_type, composition, coordinate)
    
    def from_cif(self, filename):
        cif_struc = cif(filename)
        lattice = self.para2matrix(cif_struc.cell_para)
        composition = cif_struc.composition
        coordinate = cif_struc.coordinate
        atom_type = cif_struc.atom_type
        self.from_dict(lattice, atom_type, composition, coordinate)

    def from_POSCAR(self, filename):

        f = open(filename)

        tag = f.readline()
        lattice_constant = float(f.readline().split()[0])

        # Now the lattice vectors
        a = []
        for ii in range(3):
            s = f.readline().split()
            floatvect = float(s[0]), float(s[1]), float(s[2])
            a.append(floatvect)
        lattice = np.array(a) * lattice_constant

        # Number of atoms. 
        atom_type = f.readline().split()
        comp = f.readline().split()
        composition = []
        if len(atom_type)==len(comp):
           for num in comp:
               composition.append(int(num))
        else:
           print('Value Error POSCAR symbol and composition is inconsistent')
        ac_type = f.readline().split()
        # Check if atom coordinates are cartesian or direct
        cartesian = ac_type[0].lower() == "c" or ac_type[0].lower() == "k"
        tot_natoms = sum(composition)
        coordinate = np.empty((tot_natoms, 3))
        for atom in range(tot_natoms):
            ac = f.readline().split()
            coordinate[atom] = (float(ac[0]), float(ac[1]), float(ac[2]))
        # Done with all reading
        f.close()
        if cartesian:
            coordinate *= lattice_constant
        cell_para = []
        self.coordinate = np.array(composition)
        self.from_dict(lattice, atom_type, composition, coordinate)

    def from_dict(self, lattice, atom_type, composition, coordinate):
        self.cell_matrix = np.array(lattice) 
        self.atom_type = atom_type
        self.composition = np.array(composition)
        self.coordinate = np.array(coordinate)
        self.cell_para = self.matrix2para(self.cell_matrix)
        self.rec_matrix = self.rec_lat(self.cell_matrix)
        self.name = ''
        for ele, num in zip(self.atom_type, self.composition):
            self.name += ele
            if num > 1:
               self.name += str(num)
    #def show(self, L=2):
    #    """show crystal structure"""
    #    
    #    for i in range(-L, L+1):
    #        for j in range(-L, L+1):
    #            for k in range(-L, L+1):
    #                for m in self.coordinate:
    #                    
    #                    sphere(pos=vector(m[0], m[1], m[2]), radius=R)
    
    @staticmethod
    def rec_lat(matrix):
        """ calculate the reciprocal lattice """
        rec_lat = np.zeros([3,3])
        V = np.linalg.det(matrix)
        rec_lat[0] = np.cross(matrix[1], matrix[2])/V
        rec_lat[1] = np.cross(matrix[2], matrix[0])/V
        rec_lat[2] = np.cross(matrix[0], matrix[1])/V
        return  rec_lat #* 2 * pi

    @staticmethod
    def matrix2para(matrix):
        """ 3x3 representation -> 1x6 (a, b, c, alpha, beta, gamma)"""
        cell_para = np.zeros(6)
        cell_para[0] = np.linalg.norm(matrix[0])
        cell_para[1] = np.linalg.norm(matrix[1])
        cell_para[2] = np.linalg.norm(matrix[2])
    
        cell_para[5] = angle(matrix[0], matrix[1])
        cell_para[4] = angle(matrix[0], matrix[2])
        cell_para[3] = angle(matrix[1], matrix[2])

        return cell_para

    @staticmethod
    def para2matrix(cell_para):
        """ 1x6 (a, b, c, alpha, beta, gamma) -> 3x3 representation -> """
        matrix = np.zeros([3,3])
        matrix[0][0] = cell_para[0]
        matrix[1][0] = cell_para[1]*cos(cell_para[5])
        matrix[1][1] = cell_para[1]*sin(cell_para[5])
        matrix[2][0] = cell_para[2]*cos(cell_para[4])
        matrix[2][1] = cell_para[2]*cos(cell_para[3])*sin(cell_para[4])
        matrix[2][2] = sqrt(cell_para[2]**2 - matrix[2][0]**2 - matrix[2][1]**2)
        
        return matrix
    
class cif(object):
    """a class of cif reader
    Attributes:
        wavelength: default: 1.54181a, namely Cu-Ka
        max2theta: the range of 2theta angle
        intensity: intensities for all hkl planes
        pxrd: powder diffraction data
    """

    def __init__(self, filename):
        """Return a XRD object with the proper info"""
        self.from_file(filename)
        self.parse_cell()
        self.parse_atom()
        self.apply_symops()

    def from_file(self, filename):
        cif = np.genfromtxt(filename, dtype=str, delimiter='\n')
        
        # 3 modes in each flag:  
        # 0: not started; 
        # 1: reading; 
        # 2: done
        flags = {'cell':0, 'symops':0, 'atom':0}

        atom = {}
        cell = {}
        symops = {'string':[], 'matrix':[]}

        for lines in cif:

            if 'loop_' in lines:  
                #if a _loop lines starts, the current reading flag switch to 0
                for item in flags.keys():
                    if flags[item] == 1:
                        flags[item] = 2

            elif '_cell_length_' in lines or '_cell_angle_' in lines:
                #_cell_length_a          4.77985

                flags['cell'] = 1
                cell_str = lines.split()
                item = cell_str[0].replace(' ','')
                value = float(cell_str[1].split("(")[0])
                cell[item] = value

            elif '_symmetry_equiv_pos_as_xyz' in lines:
                #_symmetry_equiv_pos_as_xyz
                flags['symops'] = 1
      
            elif '_space_group_symop_operation_xyz' in lines:
                #_space_group_symop_operation_xyz
                flags['symops'] = 1
                
            elif flags['symops'] == 1:
                #1, 'x, y, z'
                #    x, -y, z
                raw_line = lines.strip().strip("'").split(' ', 1)
                if raw_line[0].isdigit():     
                    sym_str = raw_line[1].strip("'")
                else:
                    sym_str = lines.strip().strip("'").replace(' ', '')
                sym_str = sym_str.replace("'","")
                symops['string'].append(sym_str)
                symops['matrix'].append(self.xyz2sym_ops(sym_str))

            elif '_atom_site' in lines: 
                flags['atom'] = 1
                atom_str = lines.replace(' ','')
                item = atom_str
                atom[item] = []

            elif flags['atom'] == 1:
                raw_line = lines.split()
                for i, item in enumerate(atom.keys()):
                    raw_text = raw_line[i]
                    
                    if item.find('fract')>0:
                       value = float(raw_text.split("(")[0])
                    elif item.find('symbol')>0:
                       m_symbol = re.compile("([A-Z]+[a-z]*)")
                       value = str(m_symbol.findall(raw_text)).strip("[]").strip("''")
                       #print(raw_text, value)
                    else:
                       value = raw_text
                       
                    atom[item].append(value)

            elif flags['cell'] + flags['symops'] + flags['atom'] == 6:
                break

        self.cell = cell
        self.atom = atom
        self.symops = symops
   
    def parse_cell(self):
        cell_para = np.zeros(6)
        cell = self.cell
        for item in cell.keys():
            if item.find('_length_a') > 0:
                cell_para[0] = cell[item]
            elif item.find('_length_b') > 0:
                cell_para[1] = cell[item]
            elif item.find('_length_c') > 0:
                cell_para[2] = cell[item]
            elif item.find('_angle_alpha') > 0:
                cell_para[3] = np.radians(cell[item])
            elif item.find('_angle_beta') > 0:
                cell_para[4] = np.radians(cell[item])
            elif item.find('_angle_gamma') > 0:
                cell_para[5] = np.radians(cell[item])
        self.cell_para = cell_para

    def parse_atom(self):
        atom = self.atom
        N_atom = len(atom['_atom_site_fract_x'])
        cif_xyz = np.zeros([N_atom, 3])

        for item in atom.keys():
            if item.find('_fract_x') > 0:
                cif_xyz[:,0] = np.array(atom[item])
            elif item.find('_fract_y') > 0:
                cif_xyz[:,1] = np.array(atom[item])
            elif item.find('_fract_z') > 0:
                cif_xyz[:,2] = np.array(atom[item])

        self.cif_xyz = cif_xyz

    #generates all coordinates from rotation matrices and translation vectors
    def apply_symops(self):
        fract_xyz = self.cif_xyz
        symops_matrix = self.symops['matrix']
        atom_type = self.atom['_atom_site_type_symbol']
        sym_coordinates = {}
        
        for item in atom_type:
            sym_coordinates[item] = []


        for ii,item in enumerate(atom_type):
            for mat_vec in symops_matrix:
                sym_temp = np.dot(mat_vec[0], fract_xyz[ii].transpose()) + mat_vec[1]
                sym_coordinates[item].append(sym_temp)
        self.coordinate, self.composition, self.atom_type = \
                      self.remove_duplicate(sym_coordinates)

    #remove equivalent points and keep the unique ones
    #get the numbers of atoms per species
    @staticmethod
    def remove_duplicate(sym_coordinates):
        coordinate = []
        composition = []
        atom_type = []
        for item in sym_coordinates.keys():
            atom_type.append(item)
            raw_equiv = np.array(sym_coordinates[item])
            raw_equiv = raw_equiv - np.floor(raw_equiv)
            raw_equiv = np.around(raw_equiv, 4)
            raw_equiv = np.unique(raw_equiv, axis=0)
            composition.append(len(raw_equiv))
            if coordinate == []:
                coordinate = raw_equiv
            else:
                coordinate = np.concatenate((coordinate,raw_equiv),axis=0)

        return coordinate, composition, atom_type


    #function generates rotation matrices and translation vectors from equivalent points
    @staticmethod
    def xyz2sym_ops(string):
        #rotational matrix dictionary
        rot_dic = {}
        rot_dic['x'] = np.array([1.0,0,0])
        rot_dic['y'] = np.array([0,1.0,0])
        rot_dic['z'] = np.array([0,0,1.0])
        parts = string.strip().replace(' ','').lower().split(',')
        rot_mat = []
        rot_temp = np.array([0.,0.,0.])
        trans_vec = np.array([0.,0.,0.])
        #use re module to read xyz strings
        m_rot = re.compile(r"([+-]?)([\d\.]*)/?([\d\.]*)([x-z])")
        m_trans = re.compile(r"([+-]?)([\d\.]+)/?([\d\.]*)(?![x-z])")
        for jj,item in enumerate(parts):
            #rotation matrix
            for ii,m in enumerate(m_rot.finditer(item)):
                coef = -1 if m.group(1) == '-' else 1
                if m.group(2) != '':
                    if m.group(3) != '':
                        coef *= float(m.group(2))/float(m.group(3))
                    else:
                        coef *= float(m.group(2))
                if ii == 0:                  
                    rot_temp = rot_dic[m.group(4)]*coef
                else:
                    rot_temp += rot_dic[m.group(4)]*coef
            rot_mat.append(rot_temp)
            #translation vector
            for m in m_trans.finditer(item):
                coef = -1 if m.group(1) == '-' else 1
                if m.group(3) != '':
                    coef = float(m.group(2))/float(m.group(3))
                else:
                    coef = float(m.group(2))
                trans_vec[jj] = 1.0*coef
        return (rot_mat, trans_vec)
         

class XRD(object):
    """a class of crystal structure. 
    Attributes:
        cell_para: a,b,c, alpha, beta, gamma
        cell_matrix: 3*3 matrix
        rec_matrix: reciprocal of cell matrix
        atom_type:  elemental type (e.g. Na Cl)
        composition: chemical composition (e.g., [1,1])
        coordinate: atomic positions (e.g., [[0,0,0],[0.5,0.5,0.5]])
    """

    def __init__(self, crystal, wavelength=1.54184, max2theta=180, profiling=None, 
                 fwhm = None, preferred_orientation = False, march_parameter = None):
        """Return a XRD object with the proper info"""
        self.wavelength = wavelength
        self.max2theta = np.radians(max2theta)
        self.name = crystal.name
        self.profiling = profiling
        self.fwhm = fwhm
        self.preferred_orientation = preferred_orientation
        self.march_parameter = march_parameter
        self.all_dhkl(crystal)
        self.intensity(crystal)
        self.pxrdf()     
    
        
    def by_hkl(self, hkl):
        
        # this is a simple print statement, does not need to be optimized

        """ d for any give abitray [h,k,l] index """
        id1 = np.where(np.all(self.hkl_list == np.array(hkl), axis=1 ))
        if id1 is None:
           print('This hkl is not in the given 2theta range')
        else:
           print('  2theta     d_hkl     hkl       Intensity')
           for i in id1[0]:
                print('%8.3f  %8.3f   [%2d %2d %2d] %8.2f' % \
                (np.degrees(self.theta2[i]), self.d_hkl[i], \
                 self.hkl_list[i,0], self.hkl_list[i,1], self.hkl_list[i,2], \
                 self.xrd_intensity[i] ))
    
    def all_dhkl(self, crystal):
        """ 3x3 representation -> 1x6 (a, b, c, alpha, beta, gamma)"""
        d_min = self.wavelength/sin(self.max2theta/2)/2
        # This block is to find the shortest d_hkl, 
        # for all basic directions (1,0,0), (0,1,0), (1,1,0), (1,-1,0) and so on, 26 in total 
        hkl_max = np.array([1,1,1])
        hkl_index = np.array([[[-1,-1,-1]],[[-1,-1,0]],[[-1,-1,1]],[[-1,0,-1]],[[-1,0,0]],[[-1,0,1]],[[-1,1,-1]],[[-1,1,0]],[[-1,1,1]],
                     [[0,-1,-1]],[[0,-1,0]],[[0,-1,1]],[[0,0,-1]],[[0,0,1]],[[0,1,-1]],[[0,1,0]],[[0,1,1]],
                     [[1,-1,-1]],[[1,-1,0]],[[1,-1,1]],[[1,0,-1]],[[1,0,0]],[[1,0,1]],[[1,1,-1]],[[1,1,0]],[[1,1,1]]])


        for index in hkl_index:
            d = float(np.linalg.norm( np.dot(index, crystal.rec_matrix),axis = 1))
            multiple = 1/d/d_min
            index *= round(multiple)
            for i in range(len(hkl_max)):
                if hkl_max[i] < index[0,i]:
                    hkl_max[i] = index[0,i]
        
        h1, k1, l1 = hkl_max

        h = np.arange(-h1,h1+1)
        k = np.arange(-k1,k1+1)
        l = np.arange(-l1,l1+1)

        hkl = np.array((np.meshgrid(h,k,l))).transpose()
        hkl_list = np.reshape(hkl, [len(h)*len(k)*len(l),3])
        hkl_list = hkl_list[np.where(hkl_list.any(axis=1))[0]]
        d_hkl = 1/np.linalg.norm( np.dot(hkl_list, crystal.rec_matrix), axis=1)

        shortlist = d_hkl > (d_min)
        d_hkl = d_hkl[shortlist]
        hkl_list = hkl_list[shortlist]
        sintheta = self.wavelength/2/d_hkl

        self.theta = np.arcsin(sintheta)
        self.hkl_list = hkl_list
        self.d_hkl = d_hkl

    def intensity(self, crystal):

        """
        This function calculates all that is necessary to find the intensities.
        This scheme is based off of pymatgen
        Needs improvement from different correction factors.
        """

        # open a json file with atomic scattering parameters, should eventuall go to Element class

        with open(os.path.join(os.path.dirname(__file__),
                       "atomic_scattering_params.json")) as f:
                        ATOMIC_SCATTERING_PARAMS = json.load(f)

        d0 = (1/2/self.d_hkl)**2

        # obtiain scattering parameters, atomic numbers, and occus (need to look into occus)
        coeffs = []
        zs = []
        occus = []
        
        for elem,N_elem in zip(crystal.atom_type,crystal.composition):
            for N in range(N_elem):
                c = ATOMIC_SCATTERING_PARAMS[elem]
                z = Element(elem).z
                coeffs.append(c)
                zs.append(z)
                occus.append(1) # HOW TO GENERALIZE OCCUPANCIES TERM

        coeffs = np.array(coeffs)
        self.peaks = {}
        two_thetas = []

        # self.march_parameter = 1

        TWO_THETA_TOL = 1e-5 # tolerance to find repeating angles
        SCALED_INTENSITY_TOL = 1e-3 # threshold for intensities
        
        ind = 0
        intense = []
        angle = []

        rawI = []
        rawtwo_theta = []
        for hkl, s2, theta, d_hkl in zip(self.hkl_list, d0, self.theta, self.d_hkl):
            
            # calculate the scattering factor sf
            g_dot_r = np.dot(crystal.coordinate, np.transpose([hkl])).T[0]
            sf = zs - 41.78214 * s2 * np.sum(coeffs[:, :, 0] * np.exp(-coeffs[:, :, 1] * s2), axis=1)
            
            # calculate the structure factor f
            f = np.sum(sf * occus * np.exp(2j * pi * g_dot_r))
            
            # calculate the lorentz polarization factor lf
            lf = (1 + cos(2 * theta) ** 2) / (sin(theta) ** 2 * cos(theta))

            # calculate the preferred orientation factor
            if self.preferred_orientation != False:
                G = self.march_parameter
                po = ((G * np.cos(theta))**2 + 1/G * np.sin(theta)**2)**(-3/2) 
            else:
                po = 1
    
            # calculate the intensity I
            I = (f * f.conjugate()).real
            
            # calculate 2*theta
            two_theta = degrees(2 * theta)
            
            # find where the scattered angles are equal
            ind = np.where(np.abs(np.subtract(two_thetas, two_theta)) < TWO_THETA_TOL)

            # append intensity, hkl plane, and thetas to lists
            if len(ind[0]) > 0:
                self.peaks[two_thetas[ind[0][0]]][0] += I * lf * po
                self.peaks[two_thetas[ind[0][0]]][1].append(tuple(hkl))
            else:
                self.peaks[two_theta] = [I * lf * po, [tuple(hkl)],d_hkl]
                two_thetas.append(two_theta)

        # obtain important intensities (defined by SCALED_INTENSITY_TOL)
        # and corresponding 2*theta, hkl plane + multiplicity, and d_hkl
        # print(peaks.keys())
        max_intensity = max([v[0] for v in self.peaks.values()])
        x = []
        y = []
        hkls = []
        d_hkls = []
        for k in sorted(self.peaks.keys()):
            v = self.peaks[k]
            fam = self.get_unique_families(v[1])
            if v[0] / max_intensity * 100 > SCALED_INTENSITY_TOL:
                x.append(k)
                y.append(v[0])
                
                hkls.append([{"hkl": hkl, "multiplicity": mult}
                             for hkl, mult in fam.items()])
                d_hkls.append(v[2])
               
        self.theta2 = x
        self.xrd_intensity = y
        self.hkl_list = hkls
        self.d_hkl = d_hkls

        if self.profiling != None:
            self.get_profile(max_intensity)
 
    def get_profile(self, max_intensity):
    
        """
        Here gaussian and lorentzian profiling functions are smeared over the obtained
        intensities (self.xrd_intensity) from previous calculations.
          
        Individual profiles of each intensity are superimposed, giving a net profiling function
        for the structure.
  
        The profile will be plotted with the original plot and stored as a class variable as 
        self.gpeaks.
        """

        # profile parameters
        N = 1500
        tail = 5 

        assert self.fwhm != None, "User must include a value for FWHM when profiling!"
        assert isinstance(self.fwhm, float) or isinstance(self.fwhm, int), "User must include a value for FWHM that is a number!" 
        
        # initiate profiling arrays
        self.gpeaks = np.zeros(N)
        self.gtwo_thetas = np.linspace(min(self.theta2)-tail,max(self.theta2)+tail,N)
        # gpeaks = []
        # gtwo_thetas = []

        # loop over each 2theta and intensity obtained earlier
        for theta,peak in zip(self.theta2,self.xrd_intensity):
        # for i,j in zip(range(len(self.theta2)), range(len(self.xrd_intensity))):
            # if i == 0
            # elem = [np.where(v[0] == self.xrd_intensity) for v in self.peaks.values()]
            # print(elem)
            # tmp = np.linspace(theta-tail,theta+tail,N)
            if self.profiling == 'gaussian':
                profile = self.gaussian_profile(peak,theta)
            elif self.profiling == 'lorentzian':
                profile = self.lorentzian_profile(peak,theta)
            elif self.profiling == 'psuedo_voigt':
                eta = 0.5 
                profile = eta * self.lorentzian_profile(peak,theta) + (1-eta) * self.gaussian_profile(peak,theta)
            else:
                raise NotImplementedError
            # gpeaks.append(profile)
            # gtwo_thetas.append(tmp)
            # plt.plot(tmp,profile)
            # plt.plot(theta,peak,'ko')
            # add to total profile
            self.gpeaks += profile
            # print(self.gpeaks)
        self.gpeaks/=np.max(self.gpeaks) #max_intensity
        # plt.show()
        # self.gpeaks = np.concatenate(gpeaks,axis = 0)
        # print(gtwo_thetas)
        # self.gtwo_thetas = np.concatenate(gtwo_thetas,axis = 0)
        # print(self.gtwo_thetas)
        # self.gpeaks /= max_intensity
        # self.xrd_intensity = [i/max(self.xrd_intensity) for i in self.xrd_intensity]
        # plt.plot(self.theta2,self.xrd_intensity,'o')
        # plt.plot(self.gtwo_thetas,self.gpeaks)
        # plt.show()
    def gaussian_profile(self, maxI, max_theta):
        tmp = ((self.gtwo_thetas - max_theta)/self.fwhm)**2
        return maxI * np.exp(-4*np.log(2)*tmp)
    
    def lorentzian_profile(self, maxI, max_theta):
        tmp = 1 + 4*((self.gtwo_thetas - max_theta)/self.fwhm)**2
        return maxI * 1/tmp

        
    def pxrdf(self):
        """
                Group the equivalent hkl planes together by 2\theta angle
        N*6 arrays, Angle, d_hkl, h, k, l, intensity
        """
        
        rank = range(len(self.theta2)) #np.argsort(self.theta2)
        PL = []
        last = 0
        for i in rank:
            if self.xrd_intensity[i] > 0.01:
                angle = self.theta2[i]
                if abs(angle-last) < 1e-4:
                    PL[-1][-1] += self.xrd_intensity[i]
                else:
                    PL.append([angle, self.d_hkl[i], \
                             self.hkl_list[i][0]["hkl"][0], self.hkl_list[i][0]["hkl"][1], \
                             self.hkl_list[i][0]["hkl"][2], self.xrd_intensity[i]])
                last = angle

        PL = (np.array(PL))
        PL[:,-1] = PL[:,-1]/max(PL[:,-1])
        self.pxrd = PL
        # print(PL[0],PL[-1])
    
    def plot_pxrd(self, filename=None, minimum_I = 0.01, show_hkl=True):
        """ plot PXRD """



        plt.figure(figsize=(20,10))

        if self.profiling != None:
            plt.plot(self.gtwo_thetas,self.gpeaks,'g-',label = str(self.profiling) + ' profiling')

        dx = np.degrees(self.max2theta)
        for i in self.pxrd:
            plt.bar(i[0],i[-1], color='b', width=dx/180)
            if i[-1] > minimum_I:
               if show_hkl:
                  label = self.draw_hkl(i[2:5])
                  plt.text(i[0]-dx/40, i[-1], label[0]+label[1]+label[2])
        
        ax=plt.gca()
        plt.grid()
        plt.xlim(0,dx)
        plt.xlabel('2θ')
        plt.ylabel('Intensity')
        plt.title('PXRD of '+self.name+ ', $\lambda$='+str(self.wavelength)+'$\AA$')
        
        if filename is None:
           plt.show()
        """ 
        else:
           plt.savefig(filename)
           plt.close()
        """
    def get_unique_families(self,hkls):
        """
        Returns unique families of Miller indices. Families must be permutations
        of each other.
        Args:
            hkls ([h, k, l]): List of Miller indices.
        Returns:
            {hkl: multiplicity}: A dict with unique hkl and multiplicity.
        """

       # TODO: Definitely can be sped up.
        def is_perm(hkl1, hkl2):
            h1 = np.abs(hkl1)
            h2 = np.abs(hkl2)
            return all([i == j for i, j in zip(sorted(h1), sorted(h2))])

        unique = collections.defaultdict(list)
        for hkl1 in hkls:
            found = False
            for hkl2 in unique.keys():
                if is_perm(hkl1, hkl2):
                    found = True
                    unique[hkl2].append(hkl1)
                    break
            if not found:
                unique[hkl1].append(hkl1)

        pretty_unique = {}
        for k, v in unique.items():
            pretty_unique[sorted(v)[-1]] = len(v)

        return pretty_unique
    #def plot_Laue(self, filename=None, projection=[0,0,1]):
    #    """ plot  Laue graphs"""
    #    maxI = max(self.xrd_intensity)
    #    for hkl,i in zip(self.hkl_list, self.xrd_intensity):
    #        if i/maxI > 0.01:
    #           if np.dot(hkl, np.array(projection))==0:
    #              xyz = np.dot(hkl,self.rec_matrix)
    #              angle1 = angle(xyz, projection)
    #              r = np.linalg.norm(xyz)
    #              label = self.draw_hkl(hkl)
    #              x,y = r*np.cos(angle1), r*np.sin(angle1)
    #              plt.scatter(x, y, c='b', s=i/maxI*50)
    #              plt.text(x, y, label[0]+label[1]+label[2])
   
    #    ax=plt.gca()
    #    ax.set_aspect('equal')
    #    ax.set_xticks([])
    #    ax.set_yticks([])

    #    plt.title('The simulated XRD of '+self.name)
    #    if filename is None:
    #       plt.show()
    #    else:
    #       plt.savefig(filename)

    @staticmethod
    def draw_hkl(hkl):
        """turn negative numbers in hkl to overbar"""
        hkl_str= []
        for i in hkl:
            if i<0:
               label = str(int(-i))
               label = r"$\bar{" + label + '}$'
               hkl_str.append(str(label))
            else:
               hkl_str.append(str(int(i)))

        return hkl_str

from optparse import OptionParser
import pandas as pd
from tabulate import tabulate

if __name__ == "__main__":
    #-------------------------------- Options -------------------------
    parser = OptionParser()
    parser.add_option("-m", "--hkl", dest="hkl", metavar='hkl index',
                      help="show hkl_index info, e.g., [1,0,0]")
    parser.add_option("-a", "--angle", dest="max2theta", default=180, type='float',
                      help="2theta angle range, default=180", metavar="angle")
    parser.add_option("-t", "--transform", dest="trans", metavar="files",
                      help="export file in different format")
    parser.add_option("-p", "--plot", dest="plot", default='yes',
                      help="plot pxrd, default: yes", metavar="plot")
    parser.add_option("-w", "--wavelength", dest="wavelength", default=1.54184, type='float',
                      help="wavelength: 1.54184", metavar="wave")
    parser.add_option("-c", "--crystal", dest="structure",default='',
                      help="crystal from file, cif or poscar, REQUIRED", metavar="crystal")
    parser.add_option("-f", "--full", dest="full",default='no',
                      help="show full hkl reflections", metavar="full")
    parser.add_option("-i", "--intensity", dest="minimum_I",default=0.01, type='float',
                      help="the minimum intensity to show, default 0.01", metavar="intensity")



    (options, args) = parser.parse_args()    
    if options.structure.find('cif') > 0:
       fileformat = 'cif'
    else:
       fileformat = 'POSCAR'

    test = crystal(fileformat, filename=options.structure)
    if options.plot == 'yes' or options.hkl is not None:
       xrd = XRD(test, wavelength=options.wavelength, \
                       max2theta=options.max2theta)   
       if options.full in  ['no', 'No', 'NO']:
          col_name = {'2theta': xrd.pxrd[:,0], \
                      'd_hkl':  xrd.pxrd[:,1], \
                      'h': xrd.pxrd[:,2], \
                      'k': xrd.pxrd[:,3], \
                      'l': xrd.pxrd[:,4], \
                      'Intensity':xrd.pxrd[:,5]}
       else:
          rank1 = xrd.xrd_intensity > options.minimum_I
          col_name = {'2theta':    np.degrees(xrd.theta2[rank1]), \
                      'd_hkl':     xrd.d_hkl[rank1],\
                      'h':         xrd.hkl_list[rank1,0], \
                      'k':         xrd.hkl_list[rank1,1], \
                      'l':         xrd.hkl_list[rank1,2], \
                      'Intensity': xrd.xrd_intensity[rank1] }

       df = pd.DataFrame(col_name)
       print(tabulate(df, headers='keys')) #, tablefmt='psql'))

       if options.plot == 'yes':
          xrd.plot_pxrd(filename=options.structure+'.png', minimum_I = options.minimum_I)
 
    #for name in ['alpha','gamma','delta']:
    #    fname = 'POSCAR-P3N5-'+name
    #    test = crystal('POSCAR',filename=fname)
    #    xrd = XRD(test, wavelength=0.4959, max2theta=20)   
    #    xrd.plot_pxrd(show_hkl=True, filename=name+'.png', minimum_I = 0.01)
