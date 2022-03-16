# get distance to nearest neighbour

from enum import Enum
import numpy as np

class Structure(Enum):
    bcc = 1
    fcc = 2
    diamond = 3
    default = 0

class Crystal:

    a = 0
    structure = Structure.default
    name = ""

    def __init__(self,_lat_const,_struct,_name ="Unknown"):
        self.a = _lat_const
        self.structure = _struct
        self.name = _name

    def nearest_neighbour_d(self):
        if self.structure == Structure.fcc:
            return self.a / np.sqrt(2)
        elif self.structure == Structure.bcc:
            return self.a * np.sqrt(3)/2
        elif self.structure == Structure.diamond:
            return self.a * np.sqrt(3) / 4

        else: return -1;



# --- START ---
a_given = 5.43095; #in Angstrom
str_given = Structure.diamond
name_given = "Si"

ni = Crystal(a_given,str_given,name_given)
print(ni.nearest_neighbour_d()) # WORKS